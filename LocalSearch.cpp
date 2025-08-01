#include "LocalSearch.h"
#include <algorithm>

void LocalSearch::runSearchSameDay() {
  int nbMoves = 1;
  int nbPhases = 0;
  while (nbMoves > 0 && nbPhases < 1000) {
    nbMoves = 0;
    updateMoves();
    for (unsigned int day = 2; day <= params->nbDays; day++) {
      nbMoves += mutationSameDay(day);
    }
    nbPhases++;
  }
}

void LocalSearch::melangeParcours()
{
  int j;
  for (unsigned int k = 0; k <= params->nbDays; k++) {
    for (int i = 0; i < (int)ordreParcours[k].size() - 1; i++) {
      j = i + params->rng->genrand64_int64() % ((int)ordreParcours[k].size() - i);
      std::swap(ordreParcours[k][i], ordreParcours[k][j]);
      
    }
  }
}

// updates the moves for each node which will be tried in mutationSameDay
void LocalSearch::updateMoves() {
  unsigned int client, client2;
  int size;

  for (unsigned int k = 1; k <= params->nbDays; k++) {
    // pour chaque client present dans ce jour
    for (unsigned int i = 0; i < (int)ordreParcours[k].size(); i++) {
      client = ordreParcours[k][i];
      clients[k][client]->moves.clear();
      size = params->cli[client].sommetsVoisins.size();

      for (unsigned int a1 = 0; a1 < size; a1++) {
        client2 = params->cli[client].sommetsVoisins[a1];
        if (client2 >= params->nbDepots && clients[k][client2]->estPresent)
          clients[k][client]->moves.push_back(client2);
      }
    }
  }

  // params->shuffleProches();
  melangeParcours();
}

int LocalSearch::mutationSameDay(unsigned int day) {
  dayCour = day;
  unsigned int size = ordreParcours[day].size();
  unsigned int size2;
  rechercheTerminee = false;
  bool moveEffectue = false;
  int nbMoves = 0;
  firstLoop = true;

  while (!rechercheTerminee) {
    rechercheTerminee = true;
    moveEffectue = false;
    for (unsigned int posU = 0; posU < size; posU++) {
      posU -= moveEffectue; // on retourne sur le dernier noeud si on a modifie
      nbMoves += moveEffectue;
      moveEffectue = false;
      noeudU = clients[day][ordreParcours[day][posU]];

      noeudUPred = noeudU->pred;
      x = noeudU->suiv;
      noeudXSuiv = x->suiv;
      xSuivCour = x->suiv->idx;
      routeU = noeudU->route;
      noeudUCour = noeudU->idx;
      noeudUPredCour = noeudUPred->idx;
      xCour = x->idx;

      size2 = (unsigned int) noeudU->moves.size();
      for (unsigned int posV = 0; posV < size2 && moveEffectue == 0; posV++) {
        noeudV = clients[day][noeudU->moves[posV]];
        if (!noeudV->route->nodeAndRouteTested[noeudU->idx] ||
            !noeudU->route->nodeAndRouteTested[noeudU->idx] || firstLoop)
        {
          noeudVPred = noeudV->pred;
          y = noeudV->suiv;
          noeudYSuiv = y->suiv;
          ySuivCour = y->suiv->idx;
          routeV = noeudV->route;
          noeudVCour = noeudV->idx;
          noeudVPredCour = noeudVPred->idx;
          yCour = y->idx;

          if (!moveEffectue)
            moveEffectue = mutation1();
          if (!moveEffectue)
            moveEffectue = mutation2();
          if (!moveEffectue)
            moveEffectue = mutation3();

          // les mutations 4 et 6 (switch) , sont sym�triques
          if (noeudU->idx <= noeudV->idx) {
            if (!moveEffectue)
              moveEffectue = mutation4();
            if (!moveEffectue)
              moveEffectue = mutation6();
          }
          if (!moveEffectue)
            moveEffectue = mutation5();

          // mutations 2-opt
          if (!moveEffectue)
            moveEffectue = mutation7();
          if (!moveEffectue)
            moveEffectue = mutation8();
          if (!moveEffectue)
            moveEffectue = mutation9();

          if (moveEffectue) {
            routeU->reinitSingleDayMoves();
            routeV->reinitSingleDayMoves();
          }
        }
      }

  // c'est un depot on tente l'insertion derriere le depot de ce jour
      // si il ya correlation
      if (params->isCorrelated1[noeudU->idx][depots[day][0]->idx] &&
          moveEffectue != 1)
        for (int route = 0; route < (int)depots[day].size(); route++)
        {
          noeudV = depots[day][route];
          if (!noeudV->route->nodeAndRouteTested[noeudU->idx] ||
              !noeudU->route->nodeAndRouteTested[noeudU->idx] || firstLoop)
          {
            noeudVPred = noeudV->pred;
            y = noeudV->suiv;
            noeudYSuiv = y->suiv;
            ySuivCour = y->suiv->idx;
            routeV = noeudV->route;
            noeudVCour = noeudV->idx;
            noeudVPredCour = noeudVPred->idx;
            yCour = y->idx;

            if (!moveEffectue)
              moveEffectue = mutation1();
            if (!moveEffectue)
              moveEffectue = mutation2();
            if (!moveEffectue)
              moveEffectue = mutation3();

            if (!noeudV->suiv->estUnDepot)
            {
              if (!moveEffectue)
                moveEffectue = mutation8();
              if (!moveEffectue)
                moveEffectue = mutation9();
            }

            if (moveEffectue)
            {
              routeU->reinitSingleDayMoves();
              routeV->reinitSingleDayMoves();
            }
          }
        }
      // TODO -- check that memories are working
    }
    firstLoop = false;
  }
  return nbMoves;
}

// pour un noeud, marque que tous les mouvements impliquant ce noeud ont �t�
// test�s pour chaque route du jour day
void LocalSearch::nodeTestedForEachRoute(int cli, int day)
{
  for (int route = 0; route < (int)depots[day].size(); route++)
    routes[day][route]->nodeAndRouteTested[cli] = true;
}

// enleve un client de l'ordre de parcours
void LocalSearch::removeOP(int day, int client)
{
  int it = 0;
  while (ordreParcours[day][it] != client)
  {
    it++;
  }
  ordreParcours[day][it] =
      ordreParcours[day][(int)ordreParcours[day].size() - 1];
  ordreParcours[day].pop_back();
}

// ajoute un client dans l'ordre de parcours
void LocalSearch::addOP(int day, int client)
{
  int it, temp2;
  if (ordreParcours[day].size() != 0)
  {
    it = (int)params->rng->genrand64_int64() % ordreParcours[day].size();
    temp2 = ordreParcours[day][it];
    ordreParcours[day][it] = client;
    ordreParcours[day].push_back(temp2);
  }
  else
    ordreParcours[day].push_back(client);
}

// Evaluates the current objective function of the whole solution
void LocalSearch::printInventoryLevels(std::ostream& file,bool add, std::vector<double> deliveries, double &totalCost)
{
  double inventoryClientCosts = 0.;
  double inventorySupplyCosts = 0.;
  double stockClientCosts = 0;
  double stockClientAmount=0;
  double routeCosts = 0.;
  double loadCosts = 0.;

  // Summing distance and load penalty

  for (unsigned int r = 0; r < params->nombreVehicules[1]; r++) {
    routeCosts += routes[1][r]->temps; // temps: total travel time
    if(!add)  file << "route["<<r<<"]: travel time = " << routes[1][r]->temps << "; capacity = " << routes[1][r]->capacity  << "; charge = " << std::accumulate(deliveries.begin(), deliveries.end(), 0.0) << "; depot = " << routes[1][r]->depot->idx << endl;
    routes[1][r]->printRouteData(file);
  }
  

  // Printing customer inventory and computing customer inventory cost

  double inventoryClient;
  for (unsigned int i = params->nbDepots; i < params->nbDepots + params->nbClients; i++) {
    inventoryClient = params->cli[i].startingInventory;
    if(!add) file  << "CUSTOMER " << i << " bounds (" << params->cli[i].minInventory
          << "," << params->cli[i].maxInventory << ") " << " " << params->meanDemands[i - params->nbDepots] << " " << params->stdDemands[i - params->nbDepots] << " " << params->cli[i].testDemand[params->jVal] << " " ; 

    // print the level in the morning
    if(!add) file << "[morning: " << inventoryClient;
    // print the level after receiving inventory
    inventoryClient += deliveries[i - params->nbDepots];
    if(!add) file  << ", replenishment: " << deliveries[i - params->nbDepots];
    // print the level after consumption
    double stock = std::max<double>(0,params->cli[i].testDemand[params->jVal]-inventoryClient);
    inventoryClient = std::max<double>(0,inventoryClient-params->cli[i].testDemand[params->jVal]);
    
    if(!add) file  << ", evening: " << inventoryClient << "] ";

    inventoryClientCosts += inventoryClient * params->cli[i].inventoryCost ;
    stockClientCosts += stock*params->cli[i].stockoutCost;
    stockClientAmount += stock;

    if(!add) file  << endl;
  }

  file  << endl;
  double inventorySupply = 0;
  if(!add) file  << "SUPPLIER INVENTORY ";

  inventorySupply += params->availableSupply[1];
  // print the level in the morning
  if(!add) file << "[" << inventorySupply << ",";
  for (int i = params->nbDepots; i < params->nbDepots + params->nbClients;
        i++)
    inventorySupply -= deliveries[i - params->nbDepots];
  // print the level after delivery
  if(!add) file  << inventorySupply << "] ";
  inventorySupplyCosts += inventorySupply * params->inventoryCostSupplier;
  
  if(!add) file  << endl;

  file  << "ROUTE COST: " << routeCosts << endl;
  file << "SUPPLY INVENTORY COST: " << inventorySupplyCosts << endl;
  file << "CLIENT INVENTORY COST: " << inventoryClientCosts << endl;
  file << "CLIENT STOCKOUT COST: " << stockClientCosts<<endl;
  file << "CLIENT STOCKOUT Amount: " << stockClientAmount<<endl;
  file  << "COST SUMMARY : OVERALL "
       << routeCosts + loadCosts + inventorySupplyCosts + inventoryClientCosts + stockClientCosts
       << endl;
  totalCost += routeCosts + loadCosts + inventorySupplyCosts + inventoryClientCosts + stockClientCosts;
}

// supprime le noeud
void LocalSearch::removeNoeud(Noeud *U)
{
  // mettre a jour les noeuds
  U->pred->suiv = U->suiv;
 
  U->suiv->pred = U->pred;

  U->route->updateRouteData();

  // on g�re les autres structures de donn�es
  removeOP(U->jour, U->idx);
  U->estPresent = false;

  // signifier que les insertions sur cette route ne sont plus bonnes
  U->route->initiateInsertions();

  // signifier que pour ce jour les insertions de noeuds ne sont plus bonnes
  for (int i = params->nbDepots; i < params->nbDepots + params->nbClients; i++)
    clients[U->jour][i]->coutInsertion = 1.e30;

}

void LocalSearch::addNoeud(Noeud *U)
{
  U->placeInsertion->suiv->pred = U;
  U->pred = U->placeInsertion;
  U->suiv = U->placeInsertion->suiv;
  U->placeInsertion->suiv = U;

  // et mettre a jour les routes
  U->route = U->placeInsertion->route;
  U->route->updateRouteData();

  // on g�re les autres structures de donn�es
  addOP(U->jour, U->idx);
  U->estPresent = true;

  // signifier que les insertions sur cette route ne sont plus bonnes
  U->route->initiateInsertions();

  // signifier que pour ce jour les insertions de noeuds ne sont plus bonnes
  for (int i = params->nbDepots; i < params->nbDepots + params->nbClients; i++)
    clients[U->jour][i]->coutInsertion = 1.e30;
}

// calcule pour un jour donn� et un client donn� (repr�sent� par un noeud)
// les couts d'insertion dans les differentes routes constituant ce jour
void LocalSearch::computeCoutInsertion(Noeud *client)
{
  Route *myRoute;
  client->allInsertions.clear();
  // for each route of this day
  for (unsigned int r = 0; r < routes[client->jour].size(); r++){
    // later on we can simply retrieve
    // calculate the best insertion point as well as its load

    myRoute = routes[client->jour][r];
    myRoute->evalInsertClient(client);
    client->allInsertions.push_back(myRoute->bestInsertion[client->idx]);
  }

  // eliminate dominated insertions
  client->removeDominatedInsertions(params->penalityCapa[idxScenario]);
}

double LocalSearch::evaluateCurrentCost_stockout(unsigned int client, bool test) {
  Noeud *noeudClient;
  double myCost = 0.;
  double I = params->cli[client].startingInventory;
  // Sum up the detour cost, inventory cost, and eventual excess of capacity
  if (test) std::cout << "client "  << client << std::endl;
  for (unsigned int k = 1; k <= params->nbDays; k++) {
    if (test) std::cout << "day "<< k << std::endl;
    noeudClient = clients[k][client];
    if (noeudClient->estPresent){
      // adding the inventory cost
        myCost +=
          params->cli[client].inventoryCost * 
          std::max<double> (0., I + deliveryPerDay[k][client] - params->cli[client].dailyDemand[idxScenario][k]);
          if (test && true) std::cout << "case 1 " << I << " " << deliveryPerDay[k][client] << " " << params->cli[client].dailyDemand[idxScenario][k] << " " << params->cli[client].maxInventory << std::endl;
      //stockout
        myCost +=
          params->cli[client].stockoutCost * std::max<double> (0., params->cli[client].dailyDemand[idxScenario][k] - I - deliveryPerDay[k][client]);

      //-supplier *q[]
        myCost -=  params->inventoryCostSupplier *
            (double)(params->nbDays + 1 - k) * deliveryPerDay[k][client];

      // the detour cost
        myCost +=
            params->timeCost[noeudClient->idx][noeudClient->suiv->idx] +
            params->timeCost[noeudClient->pred->idx][noeudClient->idx] -
            params->timeCost[noeudClient->pred->idx][noeudClient->suiv->idx];

      // and the possible excess capacity, the privous penalty are calculated already.
        double x1 = noeudClient->route->charge -  noeudClient->route->capacity;
        if(eq(x1,0)) x1 = 0;
        double x2=noeudClient->route->charge -
                  noeudClient->route->capacity - deliveryPerDay[k][client];
        if(eq(x2,0)) x2 = 0;
        myCost += params->penalityCapa[idxScenario] *(std::max<double>(0., x1) - std::max<double>(0., x2));
        myCost += 1000000*std::max<double> (0., I + deliveryPerDay[k][client]- params->cli[client].maxInventory);

        I = std::max<double> (0., I + deliveryPerDay[k][client] - params->cli[client].dailyDemand[idxScenario][k]);
      }
      else{     
        myCost += params->cli[client].inventoryCost *  std::max<double>(0., I - params->cli[client].dailyDemand[idxScenario][k]);
        myCost += params->cli[client].stockoutCost * std::max<double>  (0., params->cli[client].dailyDemand[idxScenario][k] - I);

        if (test && I + deliveryPerDay[k][client] > params->cli[client].maxInventory) std::cout << "case 2 " << I << " " << deliveryPerDay[k][client] << " " << params->cli[client].dailyDemand[idxScenario][k]  << " " << params->cli[client].maxInventory << std::endl;
        I = std::max<double> (0., I - params->cli[client].dailyDemand[idxScenario][k]);
        
      }
  }
  return myCost;
}

// constructeur
LocalSearch::LocalSearch(void) {}

// constructeur 2
LocalSearch::LocalSearch(Individu *_individu, Params *_params, int _idxScenario)
    : individu(_individu), params(_params), idxScenario(_idxScenario)
{
  vector<Noeud*> tempNoeud; 
  vector<Route*> tempRoute;

  Noeud *myDepot;
  Noeud *myDepotFin;
  Route *myRoute;
  
  clients.push_back(tempNoeud);
  depots.push_back(tempNoeud);
  depotsFin.push_back(tempNoeud);
  routes.push_back(tempRoute);
  
  for (unsigned int day = 1; day <= params->nbDays; day++) {
    clients.push_back(tempNoeud);
    depots.push_back(tempNoeud);
    depotsFin.push_back(tempNoeud);
    routes.push_back(tempRoute);
    // dimensionnement du champ noeuds a la bonne taille
    for (unsigned int depot = 0; depot < params->nbDepots; depot++) {
      clients[day].push_back(NULL);
    }
    for (unsigned int client = params->nbDepots; client < params->nbClients + params->nbDepots; client++) {
      clients[day].push_back(new Noeud(false, client, day, false, NULL, NULL, NULL));
    }
        
    // dimensionnement du champ depots et routes a la bonne taille       
    for (unsigned int i = 0; i < params->nombreVehicules[day]; i++) {
      myDepot = new Noeud(true, params->ordreVehicules[day][i].depotNumber, day, false, NULL, NULL, NULL);
      myDepotFin = new Noeud(true, params->ordreVehicules[day][i].depotNumber, day, false, NULL, NULL, NULL);
      myRoute = new Route(params, this, i, day, myDepot, 0, 0, params->ordreVehicules[day][i].capacity);
      myDepot->route = myRoute;
      myDepotFin->route = myRoute;
      routes[day].push_back(myRoute);
      depots[day].push_back(myDepot);
      depotsFin[day].push_back(myDepotFin);
    }
  }
              
  // initialisation de la structure ordreParcours 
  ordreParcours = std::vector<std::vector<int>>(params->nbDays + 1);

  for (unsigned int client = params->nbDepots; client < params->nbDepots + params->nbClients; client++)
    ordreParcours[0].push_back(client);
}


// destructeur
LocalSearch::~LocalSearch(void) {
  if (!clients.empty())
    for (unsigned int i = 0; i < clients.size(); i++)
      if (!clients[i].empty())
        for (unsigned int j = 0; j < clients[i].size(); j++)
          delete clients[i][j];

  if (!routes.empty())
    for (unsigned int i = 0; i < routes.size(); i++)
      if (!routes[i].empty())
        for (unsigned int j = 0; j < routes[i].size(); j++)
          delete routes[i][j];

  if (!depots.empty())
    for (unsigned int i = 0; i < depots.size(); i++)
      if (!depots[i].empty())
        for (unsigned int j = 0; j < depots[i].size(); j++)
          delete depots[i][j];
}
