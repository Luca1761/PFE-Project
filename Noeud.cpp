#include "Noeud.h"

Node::Node(void){}
Node::Node(bool _estUnDepot, unsigned int _idx, unsigned  int _jour, bool _estPresent, Node * _suiv , Node * _pred, Route * _route) 
: estUnDepot(_estUnDepot),idx(_idx),jour(_jour), estPresent(_estPresent),suiv(_suiv), pred(_pred), route(_route)
{
coutInsertion = 1.e30 ;
placeInsertion = NULL ;
place = -1 ;
}

bool mySort (Insertion i, Insertion j) 
{ 
	if (i.detour < j.detour) return true ;
	else if (i.detour > j.detour) return false ;
	else return (i.load > j.load) ;
}
void Node::removeDominatedInsertions(double penalityCapa)
{
	if (false) std::cout << penalityCapa << std::endl;
	// First order the elements by increasing detour
	std::sort (allInsertions.begin(), allInsertions.end(), mySort);

	// Then make a new structure that keeps the non-dominated elements.
	// vector <Insertion> newVector ;
	// newVector.push_back(allInsertions[0]);
	// Insertion courInsertion = allInsertions[0];

	// for (int i = 1 ; i < (int)allInsertions.size() ; i++)
	// {
	// 	if (allInsertions[i].load > courInsertion.load + 0.0001 && 
	// 		courInsertion.detour + penalityCapa * (allInsertions[i].load - courInsertion.load) > allInsertions[i].detour + 0.0001 )
	// 	{
	// 		newVector.push_back(allInsertions[i]);
	// 		courInsertion = allInsertions[i];
	// 	}
	// }

	// and replace the old structure by the new one
	// if (allInsertions.size() != newVector.size()) {
	// 	std::cout << "it happens" << std::endl;
	// 	std::cout << allInsertions.size() << " " << newVector.size() << std::endl;
	// }
	// allInsertions = newVector;
}

Node::~Node(void){}
