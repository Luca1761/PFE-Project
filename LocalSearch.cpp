#include "LocalSearch.h"
#include <algorithm>
// lance la recherche locale

void LocalSearch::runSearchTotal(bool isRepPhase)
{
  this->isRepPhase = isRepPhase;
  updateMoves();
  for (int day = 1; day <= params->nbDays; day++)
    mutationSameDay(day);

  mutationDifferentDay();
  updateMoves();
  for (int day = 1; day <= params->nbDays; day++)
    mutationSameDay(day);
}


void LocalSearch::melangeParcours()
{
  int j, temp;
  for (int k = 0; k <= params->nbDays; k++)
  {
    for (int i = 0; i < (int)ordreParcours[k].size() - 1; i++)
    {
      j = i +
          params->rng->genrand64_int64() % ((int)ordreParcours[k].size() - i);
      temp = ordreParcours[k][i];
      ordreParcours[k][i] = ordreParcours[k][j];
      ordreParcours[k][j] = temp;
      
    }
  }

  for (int i = 0; i < (int)ordreJours.size() - 1; i++)
  {
    j = i + params->rng->genrand64_int64() % ((int)ordreJours.size() - i);
    temp = ordreJours[i];
    ordreJours[i] = ordreJours[j];
    ordreJours[j] = temp;
  }
}

// updates the moves for each node which will be tried in mutationSameDay
void LocalSearch::updateMoves()
{
  int client, client2;
  int size;

  for (int k = 1; k <= params->nbDays; k++)
  {
    // pour chaque client present dans ce jour
    for (int i = 0; i < (int)ordreParcours[k].size(); i++)
    {
      client = ordreParcours[k][i];
      clients[k][client]->moves.clear();
      size = params->cli[client].sommetsVoisins.size();

      for (int a1 = 0; a1 < size; a1++)
      {
        client2 = params->cli[client].sommetsVoisins[a1];
        if (client2 >= params->nbDepots && clients[k][client2]->estPresent)
          clients[k][client]->moves.push_back(client2);
      }
    }
  }

  params->shuffleProches();
  melangeParcours();
}

int LocalSearch::mutationSameDay(int day)
{
  dayCour = day;
  int size = (int)ordreParcours[day].size();
  int size2;
  rechercheTerminee = false;
  int moveEffectue = 0;
  int nbMoves = 0;
  firstLoop = true;

  while (!rechercheTerminee)
  {
    rechercheTerminee = true;
    moveEffectue = 0;
    for (int posU = 0; posU < size; posU++)
    {
      posU -= moveEffectue; // on retourne sur le dernier noeud si on a modifi�
      nbMoves += moveEffectue;
      moveEffectue = 0;
      noeudU = clients[day][ordreParcours[day][posU]];

      noeudUPred = noeudU->pred;
      x = noeudU->suiv;
      noeudXSuiv = x->suiv;
      xSuivCour = x->suiv->idx;
      routeU = noeudU->route;
      noeudUCour = noeudU->idx;
      noeudUPredCour = noeudUPred->idx;
      xCour = x->idx;

      size2 = (int)noeudU->moves.size();
      for (int posV = 0; posV < size2 && moveEffectue == 0; posV++)
      {
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

          if (moveEffectue != 1)
            moveEffectue = mutation1();
          if (moveEffectue != 1)
            moveEffectue = mutation2();
          if (moveEffectue != 1)
            moveEffectue = mutation3();

          // les mutations 4 et 6 (switch) , sont sym�triques
          if (noeudU->idx <= noeudV->idx)
          {
            if (moveEffectue != 1)
              moveEffectue = mutation4();
            if (moveEffectue != 1)
              moveEffectue = mutation6();
          }
          if (moveEffectue != 1)
            moveEffectue = mutation5();

          // mutations 2-opt
          if (moveEffectue != 1)
            moveEffectue = mutation7();
          if (moveEffectue != 1)
            moveEffectue = mutation8();
          if (moveEffectue != 1)
            moveEffectue = mutation9();

          if (moveEffectue == 1)
          {
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

            if (moveEffectue != 1)
              moveEffectue = mutation1();
            if (moveEffectue != 1)
              moveEffectue = mutation2();
            if (moveEffectue != 1)
              moveEffectue = mutation3();

            if (!noeudV->suiv->estUnDepot)
            {
              if (moveEffectue != 1)
                moveEffectue = mutation8();
              if (moveEffectue != 1)
                moveEffectue = mutation9();
            }

            if (moveEffectue == 1)
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

// trying to change the delivery plan (lot sizing for a given customer)
int LocalSearch::mutationDifferentDay()
{
  rechercheTerminee = false;
  int nbMoves = 0;
  while(!rechercheTerminee){
    rechercheTerminee = true;
    
    for (int posU = 0; posU < params->nbClients; posU++){
      nbMoves += mutation11(ordreParcours[0][posU]);
    }
  }
  return nbMoves;
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

// change the choices of visit periods and quantity for "client"
int LocalSearch::mutation11(int client)
{
  Noeud *noeudTravail;
  double currentCost;
  // First, make sure all insertion costs are computed
  for (int k = 1; k <= params->ancienNbDays; k++){
    noeudTravail = clients[k][client]; //node* day k client
    computeCoutInsertion(noeudTravail); // detour,place (dominated) for each route
  }
  // Compute the current lot sizing solution cost (from the model point of view)
  //before optimizatio  currentCost = evaluateCurrentCost(client);
  currentCost = evaluateCurrentCost_stockout(client);
  /* Generate the structures of the subproblem */
  
  vector<vector<Insertion>> insertions = vector<vector<Insertion>>(params->nbDays);
  vector<double> quantities = vector<double>(params->nbDays);
  vector<int> breakpoints = vector<int>(params->nbDays);
  double objective;
  for (int k = 1; k <= params->nbDays; k++)
  {
    insertions[k - 1] = clients[k][client]->allInsertions;
  }
  
  unique_ptr<LotSizingSolver> lotsizingSolver(
      make_unique<LotSizingSolver>(params, insertions, client));
    
  bool ok = true;
  ok = lotsizingSolver->solveStockoutBackward();
  
  objective = lotsizingSolver->objective;
  quantities = lotsizingSolver->quantities;
  if(lt(currentCost,objective-0.01)) return 0;


  /* APPLYING THE MOVEMENT */
  // Later on we will verify whether it's an improving move or not to trigger a
  // good termination.

  // First, removing all occurences of the node.
  for (int k = 1; k <= params->ancienNbDays; k++)
  {
    noeudTravail = clients[k][client];
    if (noeudTravail->estPresent){
      
      removeNoeud(noeudTravail);
    }
    demandPerDay[k][client] = 0.;

  }
  // Then looking at the solution of the model and inserting in the good place
  for (int k = 1; k <= params->ancienNbDays; k++)
  {
    if (quantities[k - 1] > 0.0001 || (lotsizingSolver->breakpoints[k - 1]&&gt(0,lotsizingSolver->breakpoints[k - 1]->detour) )) // don't forget that in the model the index      // goes from 0 to t-1
    {
      demandPerDay[k][client] = round(quantities[k - 1]);
      
      clients[k][client]->placeInsertion = lotsizingSolver->breakpoints[k - 1]->place;
 
      addNoeud(clients[k][client]);
    }
  }

  double realCost = evaluateCurrentCost_stockout(client);
  if (fabs(realCost- objective)>0.01) {
    std::cout << "The solution doesn't give the expected cost" << std::endl;
    std::cout << "Cost: " << realCost << "; Expected cost: " << objective << std::endl;
    throw string("Cost error");
    return 0;
  }
  if (currentCost-objective >=0.01 )// An improving move has been found,
                                        // the search is not finished.
  {
    rechercheTerminee = false;
    return 1;
  }
  else
    return 0;
}

double LocalSearch::evaluateCurrentCost(int client)
{
  
  Noeud *noeudClient;
  double myCost = 0.;
  // Sum up the detour cost, inventory cost, and eventual excess of capacity
  for (int k = 1; k <= params->ancienNbDays; k++)
  {
    noeudClient = clients[k][client];
    if (noeudClient->estPresent)
    {
      // adding the inventory cost
      myCost +=
          (params->cli[client].inventoryCost - params->inventoryCostSupplier) *
          (double)(params->ancienNbDays + 1 - k) * demandPerDay[k][client];

      // the detour cost
      myCost +=
          params->timeCost[noeudClient->idx][noeudClient->suiv->idx] +
          params->timeCost[noeudClient->pred->idx][noeudClient->idx] -
          params->timeCost[noeudClient->pred->idx][noeudClient->suiv->idx];

      // and the possible excess capacity
      myCost += params->penalityCapa *
                (std::max<double>(0., noeudClient->route->charge -
                                          noeudClient->route->capacity) -
                 std::max<double>(0., noeudClient->route->charge -
                                          noeudClient->route->capacity -
                                          demandPerDay[k][client]));
    }
  }
  return myCost;
}

// Evaluates the current objective function of the whole solution
void LocalSearch::printInventoryLevels(std::ostream& file,bool add)
{
  double inventoryClientCosts = 0.;
  double inventorySupplyCosts = 0.;
  double stockClientCosts = 0;
  double stockClientAmount=0;
  double routeCosts = 0.;
  double loadCosts = 0.;

  // Summing distance and load penalty
  for (int k = 1; k <= params->ancienNbDays; k++)
  {
    for (int r = 0; r < params->nombreVehicules[k]; r++)
    {
      routeCosts += routes[k][r]->temps; // temps: total travel time
      
      if(!add)  file <<"day["<<k<<"] route["<<r<<"]: travel time = "<<routes[k][r]->temps<<endl;
      routes[k][r]->printRouteData(file);
      loadCosts +=
          params->penalityCapa *
          std::max<double>(routes[k][r]->charge - routes[k][r]->capacity,
                           0.);
    }
  }

  // Printing customer inventory and computing customer inventory cost
  if(params->isstockout){

    double inventoryClient;
    for (int i = params->nbDepots; i < params->nbDepots + params->nbClients;
         i++)
    {
      inventoryClient = params->cli[i].startingInventory;
      if(!add) file  << "CUSTOMER " << i << " bounds (" << params->cli[i].minInventory
           << "," << params->cli[i].maxInventory << ") ";
      for (int k = 1; k <= params->nbDays; k++)
      {
        // print the level in the morning
        if(!add) file << "[morning: " << inventoryClient;
        // print the level after receiving inventory
        inventoryClient += demandPerDay[k][i];
        if(!add) file  << " ,replinishment: " << demandPerDay[k][i];
        // print the level after consumption
        double stock = std::max<double>(0,params->cli[i].dailyDemand[k]-inventoryClient);
        inventoryClient = std::max<double>(0,inventoryClient-params->cli[i].dailyDemand[k]);
        
        if(!add) file  << ", everning: " << inventoryClient << "] ";

        inventoryClientCosts += inventoryClient * params->cli[i].inventoryCost ;
        stockClientCosts += stock*params->cli[i].stockoutCost;
        stockClientAmount += stock;
      }
      if(!add) file  << endl;
    }
  }
  else{
    double inventoryClient;
    for (int i = params->nbDepots; i < params->nbDepots + params->nbClients; i++)
    {
      inventoryClient = params->cli[i].startingInventory;
      if(!add) file  << "CUSTOMER " << i << " bounds (" << params->cli[i].minInventory
           << "," << params->cli[i].maxInventory << ") ";
      for (int k = 1; k <= params->nbDays; k++)
      {
        // print the level in the morning
        if(!add) file  << "[" << inventoryClient;
        // print the level after receiving inventory
        inventoryClient += demandPerDay[k][i];
        if(!add) file  << "," << inventoryClient;
        // print the level after consumption
        inventoryClient -= params->cli[i].dailyDemand[k];
        if(!add) file << "," << inventoryClient << "] ";

        inventoryClientCosts += inventoryClient * params->cli[i].inventoryCost;
      }
      if(!add) file  << endl;
    }
  }
  

  double inventorySupply = 0;
  if(!add) file  << "SUPPLIER    ";
  for (int k = 1; k <= params->nbDays; k++)
  {
    inventorySupply += params->availableSupply[k];
    // print the level in the morning
    if(!add) file  << "[" << inventorySupply << ",";
    for (int i = params->nbDepots; i < params->nbDepots + params->nbClients;
         i++)
      inventorySupply -= demandPerDay[k][i];
    // print the level after delivery
    if(!add) file  << inventorySupply << "] ";
    inventorySupplyCosts += inventorySupply * params->inventoryCostSupplier;
  }
  if(!add) file  << endl;

  file  << "ROUTE: " << routeCosts << endl;
  file << "LOAD: " << loadCosts << "SUPPLY: " << inventorySupplyCosts << endl;
  file << "CLIENT INVENTORY: " << inventoryClientCosts << endl;
  file << "CLIENT STOCKOUT: " << stockClientCosts<<endl;
  file << "CLIENT STOCKOUT Amount: " << stockClientAmount<<endl;
  file  << "COST SUMMARY : OVERALL "
       << routeCosts + loadCosts + inventorySupplyCosts + inventoryClientCosts+stockClientCosts
       << endl;
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
  for (int r = 0; r < (int)routes[client->jour].size(); r++){
    // later on we can simply retrieve
    // calculate the best insertion point as well as its load

    myRoute = routes[client->jour][r];
    myRoute->evalInsertClient(client);
    client->allInsertions.push_back(myRoute->bestInsertion[client->idx]);
  }

  // eliminate dominated insertions
  client->removeDominatedInsertions(params->penalityCapa);
}

double LocalSearch::evaluateCurrentCost_stockout(int client)
{
  Noeud *noeudClient;
  double myCost = 0.;
  double I = params->cli[client].startingInventory;
  // Sum up the detour cost, inventory cost, and eventual excess of capacity
  for (int k = 1; k <= params->ancienNbDays; k++)
  {
    noeudClient = clients[k][client];
    if (noeudClient->estPresent){
      // adding the inventory cost
        myCost +=
          params->cli[client].inventoryCost * 
          std::max<double> (0., I + demandPerDay[k][client] - params->cli[client].dailyDemand[k]);
      //stockout
        myCost +=
          params->cli[client].stockoutCost * std::max<double> (0., params->cli[client].dailyDemand[k] - I - demandPerDay[k][client]);

      //-supplier *q[]
        myCost -=  params->inventoryCostSupplier *
            (double)(params->ancienNbDays + 1 - k) * demandPerDay[k][client];

      // the detour cost
        myCost +=
            params->timeCost[noeudClient->idx][noeudClient->suiv->idx] +
            params->timeCost[noeudClient->pred->idx][noeudClient->idx] -
            params->timeCost[noeudClient->pred->idx][noeudClient->suiv->idx];

      // and the possible excess capacity, the privous penalty are calculated already.
        double x1 = noeudClient->route->charge -  noeudClient->route->capacity;
        if(eq(x1,0)) x1 = 0;
        double x2=noeudClient->route->charge -
                  noeudClient->route->capacity - demandPerDay[k][client];
        if(eq(x2,0)) x2 = 0;
        myCost += params->penalityCapa *(std::max<double>(0., x1) - std::max<double>(0., x2));
        myCost += 1000000*std::max<double> (0., I + demandPerDay[k][client]- params->cli[client].maxInventory);

        I = std::max<double> (0., I + demandPerDay[k][client] - params->cli[client].dailyDemand[k]);
      }
      else{     
        myCost += params->cli[client].inventoryCost *  std::max<double>(0., I - params->cli[client].dailyDemand[k]);
        myCost += params->cli[client].stockoutCost * std::max<double>  (0., params->cli[client].dailyDemand[k] - I);

        I = std::max<double> (0., I - params->cli[client].dailyDemand[k]);
        
      }
  }
  return myCost;
}

// constructeur
LocalSearch::LocalSearch(void) {}

// constructeur 2
LocalSearch::LocalSearch(Params *params, Individu *individu)
    : params(params), individu(individu)
{
  vector<Noeud *> tempNoeud; 
  vector<Route *> tempRoute;

  vector<int> temp2;
  Noeud *myDepot;
  Noeud *myDepotFin;
  Route *myRoute;

  clients.push_back(tempNoeud);
  depots.push_back(tempNoeud);
  depotsFin.push_back(tempNoeud);
  routes.push_back(tempRoute);

  for (int kk = 1; kk <= params->nbDays; kk++)
  {
    clients.push_back(tempNoeud);
    depots.push_back(tempNoeud);
    depotsFin.push_back(tempNoeud);
    routes.push_back(tempRoute);
    // dimensionnement du champ noeuds a la bonne taille
    for (int i = 0; i < params->nbDepots; i++)
      clients[kk].push_back(NULL);
    for (int i = params->nbDepots; i < params->nbClients + params->nbDepots;
         i++)
      clients[kk].push_back(
          new Noeud(false, i, kk, false, NULL, NULL, NULL, 0));

    // dimensionnement du champ depots et routes � la bonne taille

    for (int i = 0; i < params->nombreVehicules[kk]; i++)
    {
      myDepot = new Noeud(true, params->ordreVehicules[kk][i].depotNumber, kk,
                          false, NULL, NULL, NULL, 0);
      myDepotFin = new Noeud(true, params->ordreVehicules[kk][i].depotNumber,
                             kk, false, NULL, NULL, NULL, 0);
      myRoute = new Route(
          i, kk, myDepot, 0, 0, params->ordreVehicules[kk][i].maxRouteTime,
          params->ordreVehicules[kk][i].capacity, params, this);
      myDepot->route = myRoute;
      myDepotFin->route = myRoute;
      routes[kk].push_back(myRoute);
      depots[kk].push_back(myDepot);
      depotsFin[kk].push_back(myDepotFin);
    }
  }

  // initialisation de la structure ordreParcours 
  for (int day = 0; day <= params->nbDays; day++)
    ordreParcours.push_back(temp2);

  for (int i = params->nbDepots; i < params->nbDepots + params->nbClients; i++)
    ordreParcours[0].push_back(i);

  // initialisation de la structure ordreJours
  for (int day = 1; day <= params->nbDays; day++)
    ordreJours.push_back(day);
}


// destructeur
LocalSearch::~LocalSearch(void)
{
  if (!clients.empty())
    for (int i = 0; i < (int)clients.size(); i++)
      if (!clients[i].empty())
        for (int j = 0; j < (int)clients[i].size(); j++)
          delete clients[i][j];

  if (!routes.empty())
    for (int i = 0; i < (int)routes.size(); i++)
      if (!routes[i].empty())
        for (int j = 0; j < (int)routes[i].size(); j++)
          delete routes[i][j];

  if (!depots.empty())
    for (int i = 0; i < (int)depots.size(); i++)
      if (!depots[i].empty())
        for (int j = 0; j < (int)depots[i].size(); j++)
          delete depots[i][j];
}
