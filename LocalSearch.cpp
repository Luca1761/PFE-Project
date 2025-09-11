#include "LocalSearch.h"
#include <algorithm>

void LocalSearch::runSearchSameDay() {
  int nbMoves = 1;
  int nbPhases = 0;
  while (nbMoves > 0 && nbPhases < 1000) { // limit on nbPhases to fasten the algorithm
    nbMoves = 0;
    updateMoves();
    for (unsigned int day = 2; day <= params->nbDays; day++) {
      nbMoves += mutationSameDay(day);
    }
    nbPhases++;
  }
}

void LocalSearch::shuffleOrder() {
  unsigned int j;
  for (unsigned int k = 0; k <= params->nbDays; k++) {
    for (unsigned int i = 1; i < clientOrder[k].size(); i++) {
      j = i - 1 + (unsigned int) (params->rng->genrand64_int64() % (clientOrder[k].size() - i + 1));
      std::swap(clientOrder[k][i], clientOrder[k][j]);
    }
  }
}

void LocalSearch::updateMoves() {
  unsigned int client, client2;
  unsigned int size;

  for (unsigned int k = 1; k <= params->nbDays; k++) {
    for (unsigned int i = 0; i < clientOrder[k].size(); i++) {
      client = clientOrder[k][i];
      clients[k][client]->moves.clear();
      size = (unsigned int) params->cli[client].neighbors.size();

      for (unsigned int a1 = 0; a1 < size; a1++) {
        client2 = params->cli[client].neighbors[a1];
        if (client2 >= params->nbDepots && clients[k][client2]->isPresent) // we add clients neighbors that are present
          clients[k][client]->moves.push_back(client2);
      }
    }
  }

  params->shuffleNeighbors();
  shuffleOrder();
}

int LocalSearch::mutationSameDay(unsigned int day) {
  currDay = day;
  stopResearch = false;
  firstLoop = true;
  unsigned int size = (unsigned int) clientOrder[day].size();
  unsigned int size2;
  bool moveDone = false;
  int nbMoves = 0;

  while (!stopResearch) {
    stopResearch = true;
    moveDone = false;
    for (unsigned int posU = 0; posU < size; posU++) {
      posU -= moveDone; // if we did a move, we go back to last node
      nbMoves += moveDone;
      moveDone = false;
      nodeU = clients[day][clientOrder[day][posU]];

      nodeUPrev = nodeU->prev;
      x = nodeU->next;
      nodeXNext = x->next;
      idxXNext = x->next->idx;
      routeU = nodeU->route;
      idxNodeU = nodeU->idx;
      idxNodeUPrev = nodeUPrev->idx;
      idxX = x->idx;

      size2 = (unsigned int) nodeU->moves.size();
      for (unsigned int posV = 0; posV < size2 && !moveDone; posV++) {
        nodeV = clients[day][nodeU->moves[posV]];
        if (!nodeV->route->nodeAndRouteTested[nodeU->idx] ||
            !nodeU->route->nodeAndRouteTested[nodeU->idx] || firstLoop)
        {
          nodeVPrev = nodeV->prev;
          y = nodeV->next;
          nodeYNext = y->next;
          idxYNext = y->next->idx;
          routeV = nodeV->route;
          idxNodeV = nodeV->idx;
          idxNodeVPrev = nodeVPrev->idx;
          idxY = y->idx;
          
          // LNS
          if (!moveDone) moveDone = mutation1();
          if (!moveDone) moveDone = mutation2();
          if (!moveDone) moveDone = mutation3();

          // mutations 4 and 6 (switch) are symetrical
          if (nodeU->idx <= nodeV->idx) {
            if (!moveDone) moveDone = mutation4();
            if (!moveDone) moveDone = mutation6();
          }
          if (!moveDone) moveDone = mutation5();

          // mutations 2-OPT
          if (!moveDone) moveDone = mutation7();
          if (!moveDone) moveDone = mutation8();
          if (!moveDone) moveDone = mutation9();

          if (moveDone) {
            routeU->reinitSingleDayMoves();
            routeV->reinitSingleDayMoves();
          }
        }
      }

      // we try to insert just after the depot on this day (if it's correlated)
      if (params->isCorrelated[nodeU->idx][depots[day][0]->idx] && !moveDone)
        for (unsigned int depot = 0; depot < depots[day].size(); depot++) {
          nodeV = depots[day][depot];
          if (!nodeV->route->nodeAndRouteTested[nodeU->idx] ||
              !nodeU->route->nodeAndRouteTested[nodeU->idx] || firstLoop)
          {
            nodeVPrev = nodeV->prev;
            y = nodeV->next;
            nodeYNext = y->next;
            idxYNext = y->next->idx;
            routeV = nodeV->route;
            idxNodeV = nodeV->idx;
            idxNodeVPrev = nodeVPrev->idx;
            idxY = y->idx;

            if (!moveDone) moveDone = mutation1();
            if (!moveDone) moveDone = mutation2();
            if (!moveDone) moveDone = mutation3();

            //mutations 4, 5, 6 and 7 cannot be done with a depot

            if (!nodeV->next->isADepot) {
              if (!moveDone) moveDone = mutation8();
              if (!moveDone) moveDone = mutation9();
            }

            if (moveDone) {
              routeU->reinitSingleDayMoves();
              routeV->reinitSingleDayMoves();
            }
          }
        }
    }
    firstLoop = false;
  }
  return nbMoves;
}

void LocalSearch::removeCO(unsigned int day, unsigned int client) {
  if (clientOrder[day].empty()) throw std::string("ERROR SIZE");
  unsigned int it = 0;
  while (clientOrder[day][it] != client) {
    it++;
  }
  clientOrder[day][it] = clientOrder[day][clientOrder[day].size() - 1];
  clientOrder[day].pop_back();
}

void LocalSearch::addCO(unsigned int day, unsigned int client) {
  unsigned int it;
  if (!clientOrder[day].empty()) {
    it = (unsigned int) (params->rng->genrand64_int64() % clientOrder[day].size());
    clientOrder[day].push_back(clientOrder[day][it]);
    clientOrder[day][it] = client;
  }
  else
    clientOrder[day].push_back(client);
}

void LocalSearch::printInventoryLevels(std::ostream& file, std::vector<double> deliveries, double &totalCost) {
  double inventoryClientCosts = 0.;
  double inventorySupplyCosts = 0.;
  double stockoutClientCosts = 0;
  double stockoutClientAmount = 0;
  double routeCosts = 0.;

  // Summing distance and load penalty
  for (unsigned int r = 0; r < params->vehicleNumber[1]; r++) {
    routeCosts += routes[1][r]->time; //  total travel time
    double loadRoute = 0.0;

    Node* node = routes[1][r]->depot ;
    while (!node->next->isADepot) {
      node = node->next;
      loadRoute += deliveries[node->idx - params->nbDepots]; // compute the load of the route with average deliveries
    }

    if (loadRoute > routes[1][r]->capacity) {
      std::cout << "INVALID CHARGE for route " + to_string(r) << std::endl;
      throw std::string("INVALID CHARGE for route " + to_string(r));
    }
    file << "route["<<r<<"]: travel time = " << routes[1][r]->time << "; capacity = " << routes[1][r]->capacity << "; charge = " << loadRoute << "; depot = " << routes[1][r]->depot->idx << endl;
    routes[1][r]->printRoute(file);
  }
  file << endl;
  

  // Printing customer inventory and computing customer inventory cost
  double inventoryClient;
  for (unsigned int i = params->nbDepots; i < params->nbDepots + params->nbClients; i++) {
    inventoryClient = params->cli[i].startingInventory;
    file  << "CUSTOMER " << i << " bounds (" << "0"
          << "," << params->cli[i].maxInventory << ")" << " | Av: " 
          << params->meanDemands[i - params->nbDepots] << " | Std: " 
          << params->stdDemands[i - params->nbDepots] << " | True: " 
          << params->cli[i].trueDemand[params->jVal] << " " ; 

    // print the level in the morning
    file << "[morning: " << inventoryClient;
    // print the level after receiving inventory
    inventoryClient += deliveries[i - params->nbDepots];
    file  << ", replenishment: " << deliveries[i - params->nbDepots];
    // print the level after consumption
    double stockout = std::max<double>(0, params->cli[i].trueDemand[params->jVal] - inventoryClient);
    inventoryClient = std::max<double>(0, inventoryClient - params->cli[i].trueDemand[params->jVal]);
    
    file  << ", evening: " << inventoryClient << "] ";

    inventoryClientCosts += inventoryClient * params->cli[i].inventoryCost ;
    stockoutClientCosts += stockout * params->cli[i].stockoutCost;
    stockoutClientAmount += stockout;

    file  << endl;
  }

  file  << endl;

  double inventorySupply = 0;
  file  << "SUPPLIER INVENTORY ";
  inventorySupply += params->availableSupply[1];
  // print the level in the morning
  file << "[morning: " << inventorySupply;
  for (unsigned int i = params->nbDepots; i < params->nbDepots + params->nbClients; i++)
    inventorySupply -= deliveries[i - params->nbDepots];
  // print the level after delivery
  file << ", evening: " <<  inventorySupply << "] ";
  inventorySupplyCosts += inventorySupply * params->inventoryCostSupplier;
  
  file  << endl;

  file  << "ROUTE COST: " << routeCosts << endl;
  file << "SUPPLY INVENTORY COST: " << inventorySupplyCosts << endl;
  file << "CLIENT INVENTORY COST: " << inventoryClientCosts << endl;
  file << "CLIENT STOCKOUT COST: " << stockoutClientCosts << endl;
  file << "CLIENT STOCKOUT Amount: " << stockoutClientAmount << endl;
  file  << "COST SUMMARY : THIS DAY "
       << routeCosts + inventorySupplyCosts + inventoryClientCosts + stockoutClientCosts
       << endl;
  totalCost += routeCosts + inventorySupplyCosts + inventoryClientCosts + stockoutClientCosts;
}

void LocalSearch::removeNode(Node *U) {
  // update nodes
  U->prev->next = U->next;
  U->next->prev = U->prev;

  // update related route
  U->route->updateRouteData();

  // remove from client order and its presence
  removeCO(U->day, U->idx);
  U->isPresent = false;

  // insertions on this route are not good anymore
  U->route->initiateInsertions();
}

void LocalSearch::addNode(Node *U) {
  // update nodes
  U->placeInsertion->next->prev = U;
  U->prev = U->placeInsertion;
  U->next = U->placeInsertion->next;
  U->placeInsertion->next = U;

  // update route
  U->route = U->placeInsertion->route;
  U->route->updateRouteData();

  // add to client order and set its presence to true
  addCO(U->day, U->idx);
  U->isPresent = true;

  // insertions on this route are not good anymore
  U->route->initiateInsertions();
}

bool mySort (Insertion i, Insertion j) { 
	if (i.detour < j.detour) return true ;
	else if (i.detour > j.detour) return false ;
	else return (i.load > j.load) ;
}

void LocalSearch::computeInsertionCost(Node *client) {
  Route *myRoute;
  client->allInsertions.clear();
  // for each route of this day
  for (unsigned int r = 0; r < routes[client->day].size(); r++){
    // calculate the best insertion point as well as its load
    myRoute = routes[client->day][r];
    myRoute->evalInsertClient(client);
    client->allInsertions.push_back(myRoute->bestInsertion[client->idx]);
  }
  std::sort(client->allInsertions.begin(), client->allInsertions.end(), mySort);
}

double LocalSearch::evaluateCurrentClientCost(unsigned int client) {
  Node *clientNode;
  double cost = 0.;
  double I = params->cli[client].startingInventory;
  // Sum up the detour cost, inventory cost, and eventual excess of capacity
  for (unsigned int k = 1; k <= params->nbDays; k++) {
    clientNode = clients[k][client];
    if (clientNode->isPresent){
      // adding the inventory cost
      cost +=
        params->cli[client].inventoryCost * 
        std::max<double> (0., I + deliveryPerDay[k][client] - params->cli[client].dailyDemand[idxScenario][k]);
      
      //stockout
      cost +=
        params->cli[client].stockoutCost * std::max<double> (0., params->cli[client].dailyDemand[idxScenario][k] - I - deliveryPerDay[k][client]);

      //-supplier *q[]
      cost -=  params->inventoryCostSupplier *
          (double)(params->nbDays + 1 - k) * deliveryPerDay[k][client];

      // the detour cost
      cost +=
          params->timeCost[clientNode->prev->idx][clientNode->idx] +
          params->timeCost[clientNode->idx][clientNode->next->idx] -
          params->timeCost[clientNode->prev->idx][clientNode->next->idx];

      // and the possible excess capacity, the privous penalty are calculated already.
      double x1 = clientNode->route->load - clientNode->route->capacity;
      if(eq(x1,0)) x1 = 0;
      double x2 = clientNode->route->load - deliveryPerDay[k][client] - clientNode->route->capacity;
      if(eq(x2, 0)) x2 = 0;
      cost += params->penalityCapa[idxScenario] *(std::max<double>(0., x1) - std::max<double>(0., x2));

      double maxDeliverable = (params->endDayInventories) ? params->cli[client].dailyDemand[idxScenario][k] + params->cli[client].maxInventory : params->cli[client].maxInventory;
      cost += 1000000 * std::max<double> (0., I + deliveryPerDay[k][client] - maxDeliverable);

      I = std::max<double> (0., I + deliveryPerDay[k][client] - params->cli[client].dailyDemand[idxScenario][k]);
    } else{     
      cost += params->cli[client].inventoryCost *  std::max<double>(0., I - params->cli[client].dailyDemand[idxScenario][k]);
      cost += params->cli[client].stockoutCost * std::max<double>  (0., params->cli[client].dailyDemand[idxScenario][k] - I);

      I = std::max<double> (0., I - params->cli[client].dailyDemand[idxScenario][k]);
      
    }
  }
  return cost;
}

LocalSearch::LocalSearch(void) {}

LocalSearch::LocalSearch(Individual *_indiv, Params *_params, unsigned int _idxScenario)
    : indiv(_indiv), params(_params), idxScenario(_idxScenario)
{

  Node *myDepot;
  Node *myEndDepot;
  Route *myRoute;

  clients = vector<vector<Node*>>(params->nbDays + 1);
  depots = vector<vector<Node*>>(params->nbDays + 1);
  endDepots = vector<vector<Node*>>(params->nbDays + 1);
  routes = vector<vector<Route*>>(params->nbDays + 1);
  
  for (unsigned int day = 1; day <= params->nbDays; day++) {
    for (unsigned int depot = 0; depot < params->nbDepots; depot++) {
      clients[day].push_back(NULL);
    }
    for (unsigned int client = params->nbDepots; client < params->nbClients + params->nbDepots; client++) {
      clients[day].push_back(new Node(false, client, day, false, NULL, NULL, NULL));
    }
   
    for (unsigned int r = 0; r < params->vehicleNumber[day]; r++) {
      myDepot = new Node(true, params->vehicleOrder[day][r].depotNumber, day, false, NULL, NULL, NULL);
      myEndDepot = new Node(true, params->vehicleOrder[day][r].depotNumber, day, false, NULL, NULL, NULL);
      myRoute = new Route(params, this, r, day, myDepot, 0, 0, params->vehicleOrder[day][r].capacity);
      myDepot->route = myRoute;
      myEndDepot->route = myRoute;
      routes[day].push_back(myRoute);
      depots[day].push_back(myDepot);
      endDepots[day].push_back(myEndDepot);
    }
  }
              
  // initialize clientOrder
  clientOrder = std::vector<std::vector<unsigned int>>(params->nbDays + 1);

  for (unsigned int client = params->nbDepots; client < params->nbDepots + params->nbClients; client++)
    clientOrder[0].push_back(client);
}

LocalSearch::~LocalSearch(void) {
  // delete pointers if they exist
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
