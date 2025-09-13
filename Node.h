/*                       Algorithme - HGSADC                         */
/*                    Propriete de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef NODE_H
#define NODE_H

#include "Route.h"
#include<iostream>
using namespace std;

class Route ;
class Node
{
    
public :
  
// status of the node (customer or supplier)
bool isADepot;

// index of the represented customer/supplier 
unsigned int idx;

// place in its potential route
int place;

// node day
unsigned int day;

// true if this customer is delivered by a route on this day
bool isPresent;

// next customer/depot on its potential route
Node* next;

// next customer/depot on its potential route
Node* prev;

// its potential route 
Route* route;

// total load of the road until this node (including itself)
double previousLoad;

// list of possible best insertions in the different routes
vector<Insertion> allInsertions;

// final choosen place to insert this node
Node* placeInsertion;

// possible moves
vector<unsigned int> moves;

void removeDominatedInsertions (double penalityCapa);

// default constructor	
Node(void);

// true constructor
Node(bool _isADepot, unsigned int _idx, unsigned int _day, bool _isPresent, Node * _next , Node * _prev, Route * _route);

// destructor
~Node(void);
};

#endif
