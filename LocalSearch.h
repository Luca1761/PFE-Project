/*                       Algorithme - HGSADC                         */
/*                    Propriete de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

// local search class
// this class is linked to an individual from the population

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <vector>
#include "Individual.h"
#include "LotSizingSolver.h"
#include "Node.h"
#include "Route.h"

using namespace std;

class Individual;

class LocalSearch
{
private:
Individual *indiv;

public:
  Params *params;

  unsigned int idxScenario;

  bool stopResearch;

  // order to explore clients of each day
  vector<vector<unsigned int>> clientOrder;

  // add a client in previous vector
  void addCO(unsigned int day, unsigned int client);

  // remove a client
  void removeCO(unsigned int day, unsigned int client);
  
  // shuffle the order
  void shuffleOrder();
  
  // updates the moves for each node which will be tried in mutationSameDay
  void updateMoves();
  
  // verify if we are in the first loop of local search
  bool firstLoop;

  Node *nodeU;
  Node *nodeV;
  Node *x; // node after U (not called nodeUNext because we also need nodeXNext)
  Node *y; // node after V (not called nodeVNext because we also need nodeYNext)
  Node *nodeUPrev;
  Node *nodeVPrev;
  Node *nodeXNext;
  Node *nodeYNext;
  Route *routeU;
  Route *routeV;
  unsigned int idxNodeU, idxNodeUPrev, idxX, idxXNext, idxYNext, idxNodeV,
      idxNodeVPrev, idxY;
  unsigned int currDay;

  /* nbDays * nbClients sized vector, clients[day][i] gives every information about client i on this day
  (even if this client is not present on this day)
  */

  vector<vector<Node*>> clients;

  // nodes related to used depots
  vector<vector<Node*>> depots;

  // nodes related to routes ends
  vector<vector<Node*>> endDepots;

  // vector with routes information
  vector<vector<Route*>> routes;

  // deliveryPerDay[i][j] -> The load to be delivered to each customer [j] on day [i]
  vector<vector<double>> deliveryPerDay;

  // full procedure doing every possible mutation (pattern changes and swaps)
  // returns the number of move it made
  int mutationSameDay(unsigned int day);

  // launch the previous procedure on every day
  void runSearchSameDay();

  ////////////////////////////  MUTATION PART  ////////////////////////////

  // If nodeU is a client node, remove nodeU then insert it after nodeV
  bool mutation1();

  // If nodeU and x are clients, remove them then insert (nodeU,x) after
  // nodeV
  bool mutation2();

  // If nodeU and x are clients, remove them then insert (x,nodeU) after
  // nodeV
  bool mutation3();

  /* SWAP */

  // If nodeU and nodeV are clients, swap them
  bool mutation4();

  // If nodeU, x and nodeV are clients, swap (nodeU,x) and nodeV
  bool mutation5();

  // If (nodeU,x) and (nodeV,y) are cliens, swap (nodeU,x) and
  // (nodeV,y)
  bool mutation6();

  /* 2-OPT and 2-OPT* */

  // If T(nodeU) = T(nodeV) , replace (nodeU,x) and (nodeV,y) by
  // (nodeU,nodeV) and (x,y)
  bool mutation7();

  // If T(nodeU) != T(nodeV) , replace (nodeU,x) and (nodeV,y) by
  // (nodeU,nodeV) and (x,y)
  bool mutation8();

  // If T(nodeU) != T(nodeV) , replace (nodeU,x) and (nodeV,y) by
  // (nodeU,y) and (nodeV,x)
  bool mutation9();

  //////////////////////////////////////////////////////////////////////

  // Evaluates the current objective function from the model
  double evaluateCurrentClientCost(unsigned int client);

  // Prints some useful information on the current solution (and write those in a given file)
  // Also evaluates the current objective function of the whole solution with true demands
  void printInventoryLevels(std::ostream& file, std::vector<double> deliveries, double &totalCost);

  /* Routines to update the solution */

  // insert node U just after node V in its route
  void insertNode(Node *U, Node *V);

  // swap node U and V in their respective routes
  void swapNode(Node *U, Node *V);

  // remove a given node
  void removeNode(Node *U);

  // add a given node in the place indicated by Node->placeRoute
  void addNode(Node *U);

  // for a given client/day, computes its insertion cost in the different possible routes on this day
  void computeInsertionCost(Node *client);

  // default constructor
  LocalSearch();
  
  // constructor (creates all the required nodes for depots and clients)
  LocalSearch(Individual* _indiv, Params* _params, unsigned int _idxScenario);

  // destructor
  ~LocalSearch(void);
};

#endif
