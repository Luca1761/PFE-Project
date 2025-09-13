#include "Node.h"

Node::Node(void){}

Node::Node(bool _isADepot, unsigned int _idx, unsigned  int _day, bool _isPresent, Node * _next , Node * _prev, Route * _route) 
: isADepot(_isADepot),idx(_idx),day(_day), isPresent(_isPresent), next(_next), prev(_prev), route(_route)
{
	placeInsertion = NULL ;
	place = -1 ;
}

void Node::removeDominatedInsertions (double penalityCapa) {
	// Then make a new structure that keeps the non-dominated elements.
	vector <Insertion> newVector ;
	newVector.push_back(allInsertions[0]);
	Insertion courInsertion = allInsertions[0];

	for (auto &insertion : allInsertions) {
		if (insertion.load > courInsertion.load + 0.0001 && 
			courInsertion.detour + penalityCapa * (insertion.load - courInsertion.load) > insertion.detour + 0.0001 )
		{
			newVector.push_back(insertion);
			courInsertion = insertion;
		}
	}
	// and replace the old structure by the new one
	allInsertions = newVector;
}

Node::~Node(void){}
