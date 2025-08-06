/*                       Algorithme - HGSADC                         */
/*                    Propriete de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef ROUTE_H
#define ROUTE_H

#include "Params.h"
#include "Noeud.h"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
using namespace std ;

class Node ;
class LocalSearch ;

class Route
{

private:

// instance parameters
Params * params;

// local search related to this route
LocalSearch * myLS ;

public:

// route index
unsigned int idx;

// day associated to the route
unsigned int day;

// related depot
Node * depot;

// total traveled distance (time) on the road
double time;

// total load on the road
double load;

// maximum capacity of the route
double capacity;

// compute the possible excedent of charge
inline double excedentCharge(double _load) {
	return std::max<double>(0.0, _load - capacity);
}

// update load and time of the route
void updateRouteData();

// print the route
void printRoute(std::ostream& file);

// for each client, stores the best insertion of this one in this route
vector<Insertion> bestInsertion;

// for each client, stores if every move related to this node and this route have been tested without success
vector<bool> nodeAndRouteTested;

// for a given client, find the best position to a less-cost insertion in the route
void evalInsertClient(Node * U);

// initialize insertions with default values
void initiateInsertions();

// re-initialize moves related to thiw route
void reinitSingleDayMoves();

// default constructor
Route(void);

// constructor
Route(Params* _params, LocalSearch* _myLS, unsigned int _idx, unsigned int _day, Node * _depot, double _time, double _load, double _capacity);

// destructor
~Route(void);
};

#endif
