#include "Node.h"

Node::Node(void){}

Node::Node(bool _isADepot, unsigned int _idx, unsigned  int _day, bool _isPresent, Node * _next , Node * _prev, Route * _route) 
: isADepot(_isADepot),idx(_idx),day(_day), isPresent(_isPresent), next(_next), prev(_prev), route(_route)
{
	placeInsertion = NULL ;
	place = -1 ;
}

Node::~Node(void){}
