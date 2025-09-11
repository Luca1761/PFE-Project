#include "LocalSearch.h"

void LocalSearch::insertNode(Node* U, Node* V) {
	// update nodes (deleter U and insert it after V)
	U->prev->next = U->next ;
	U->next->prev = U->prev ;
	V->next->prev = U ;
	U->prev = V ;
	U->next = V->next ;
	V->next = U ;

	// update routes
	U->route = routeV ;
	routeU->updateRouteData();
	routeV->updateRouteData();
}

void LocalSearch::swapNode(Node* U, Node* V) {
	// update nodes
	Node* VPrev = V->prev ;
	Node* VNext = V->next ;
	Node* UPrev = U->prev ;
	Node* UNExt = U->next ;
	Route* myRouteU = U->route ;
	Route* myRouteV = V->route ;

	UPrev->next = V ;
	UNExt->prev = V ;
	VPrev->next = U ;
	VNext->prev = U ;

	U->prev = VPrev ;
	U->next = VNext ;
	V->prev = UPrev ;
	V->next = UNExt ;

	// update routes
	U->route = myRouteV ;
	V->route = myRouteU ;
	U->route->updateRouteData();
	V->route->updateRouteData();
}

bool LocalSearch::mutation1 () {
	double costSuppU = params->timeCost[idxNodeUPrev][idxX] 
	- params->timeCost[idxNodeUPrev][idxNodeU]  
	- params->timeCost[idxNodeU][idxX];

	double costSuppV = params->timeCost[idxNodeV][idxNodeU] 
	+ params->timeCost[idxNodeU][idxY] 
	- params->timeCost[idxNodeV][idxY];

	if (routeU != routeV) {
		costSuppU += routeU->excedentCharge(routeU->load - deliveryPerDay[currDay][idxNodeU])*params->penalityCapa[idxScenario]
		- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;

		costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[currDay][idxNodeU])*params->penalityCapa[idxScenario]
		- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if (costSuppU + costSuppV > -0.0001) return false;
	if (idxNodeU == idxY) return false;

	// update nodes
	insertNode(nodeU,nodeV);

	stopResearch = false ; 
	return true;
}

// If noeudU and x are clients, remove them then insert (noeudU,x) after noeudV
// teste si x n'est pas un depot , et si x different de noeudV, et si noeudU pas deja apres noeudV
bool LocalSearch::mutation2 ()
{
	double costSuppU = params->timeCost[idxNodeUPrev][idxXNext] 
	- params->timeCost[idxNodeUPrev][idxNodeU] 
	- params->timeCost[idxNodeU][idxX] 
	- params->timeCost[idxX][idxXNext];

	double costSuppV = params->timeCost[idxNodeV][idxNodeU] 
	+ params->timeCost[idxNodeU][idxX] 
	+ params->timeCost[idxX][idxY] 
	- params->timeCost[idxNodeV][idxY];

	if (routeU != routeV) {
	costSuppU += routeU->excedentCharge(routeU->load - deliveryPerDay[currDay][idxNodeU] - deliveryPerDay[currDay][idxX])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[currDay][idxNodeU] + deliveryPerDay[currDay][idxX])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( nodeU == y || nodeV == x || x->isADepot ) { return 0 ;}

	// mettre a jour les noeuds
	insertNode(nodeU,nodeV);
	insertNode(x,nodeU);

	stopResearch = false ; 
	return 1 ;
}

// If noeudU and x are clients, remove them then insert (x,noeudU) after noeudV
// teste si x n'est pas un depot , et si x different de noeudV, et si noeudU pas d�ja apres noeudV
bool LocalSearch::mutation3 ()
{
	double costSuppU = params->timeCost[idxNodeUPrev][idxXNext] 
	- params->timeCost[idxNodeUPrev][idxNodeU] 
	- params->timeCost[idxNodeU][idxX] 
	- params->timeCost[idxX][idxXNext];

	double costSuppV = params->timeCost[idxNodeV][idxX] 
	+ params->timeCost[idxX][idxNodeU] 
	+ params->timeCost[idxNodeU][idxY] 
	- params->timeCost[idxNodeV][idxY];

	if (routeU != routeV) {
	costSuppU += routeU->excedentCharge(routeU->load - deliveryPerDay[currDay][idxNodeU] - deliveryPerDay[currDay][idxX])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[currDay][idxNodeU] + deliveryPerDay[currDay][idxX])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) return 0;
	if ( nodeU == y ||  x == nodeV || x->isADepot ) return 0;

	// mettre a jour les noeuds
	insertNode(x,nodeV);
	insertNode(nodeU,x);

	stopResearch = false ; 
	return 1 ;
}
// SWAP
// If noeudU and noeudV are clients, swap noeudU and noeudV
// sauf si noeudU et noeudV se succedent
bool LocalSearch::mutation4 ()
{
	double costSuppU = params->timeCost[idxNodeUPrev][idxNodeV] 
	+ params->timeCost[idxNodeV][idxX]
	- params->timeCost[idxNodeUPrev][idxNodeU] 
	- params->timeCost[idxNodeU][idxX];

	double costSuppV = params->timeCost[idxNodeVPrev][idxNodeU] 
	+ params->timeCost[idxNodeU][idxY]
	- params->timeCost[idxNodeVPrev][idxNodeV] 
	- params->timeCost[idxNodeV][idxY];

	if (routeU != routeV) {
	costSuppU += routeU->excedentCharge(routeU->load + deliveryPerDay[currDay][idxNodeV] - deliveryPerDay[currDay][idxNodeU])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[currDay][idxNodeU] - deliveryPerDay[currDay][idxNodeV])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( idxNodeU == idxNodeVPrev || idxNodeU == idxY) { return 0 ;}

	// mettre a jour les noeuds
	swapNode(nodeU, nodeV) ;

	stopResearch = false ; 
	return 1 ;
}

// If noeudU, x and noeudV are clients, swap (noeudU,x) and noeudV
bool LocalSearch::mutation5 ()
{
	// on ne fait pas le cas ou x et noeudVCour se suivent
	// car il faut traiter autrement
	// et la mutation 2 entre (noeudUCour,x) et noeudVCour fait le meme travail correctement

	double costSuppU = params->timeCost[idxNodeUPrev][idxNodeV] 
	+ params->timeCost[idxNodeV][idxXNext]
	- params->timeCost[idxNodeUPrev][idxNodeU] 
	- params->timeCost[idxNodeU][idxX] 
	- params->timeCost[idxX][idxXNext];

	double costSuppV = params->timeCost[idxNodeVPrev][idxNodeU] 
	+ params->timeCost[idxX][idxY]
	+ params->timeCost[idxNodeU][idxX]
	- params->timeCost[nodeVPrev->idx][idxNodeV] 
	- params->timeCost[idxNodeV][idxY];
	
	if (routeU != routeV) {
	costSuppU += routeU->excedentCharge(routeU->load + deliveryPerDay[currDay][idxNodeV] - deliveryPerDay[currDay][idxNodeU] - deliveryPerDay[currDay][idxX])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;

	costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[currDay][idxNodeU] + deliveryPerDay[currDay][idxX] - deliveryPerDay[currDay][idxNodeV])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( nodeU == nodeVPrev || x == nodeVPrev || nodeU == y || x->isADepot ) { return 0 ;}

	// mettre a jour les noeuds
	swapNode(nodeU, nodeV) ;
	insertNode(x,nodeU);

	stopResearch = false ; 
	return 1 ;
}
// If (noeudU,x) and (noeudV,y) are clients, swap (noeudU,x) and (noeudV,y)
bool LocalSearch::mutation6 ()
{
	double costSuppU = params->timeCost[idxNodeUPrev][idxNodeV]  
	+ params->timeCost[idxNodeV][idxY]
	+ params->timeCost[idxY][idxXNext]
	- params->timeCost[idxNodeUPrev][idxNodeU] 
	- params->timeCost[idxNodeU][idxX] 
	- params->timeCost[idxX][idxXNext];

	double costSuppV = params->timeCost[idxNodeVPrev][idxNodeU] 
	+ params->timeCost[idxNodeU][idxX]
	+ params->timeCost[idxX][idxYNext]
	- params->timeCost[idxNodeVPrev][idxNodeV] 
	- params->timeCost[idxNodeV][idxY]
	- params->timeCost[idxY][idxYNext];
	
	if (routeU != routeV) {
	costSuppU += routeU->excedentCharge(routeU->load + deliveryPerDay[currDay][idxNodeV] + deliveryPerDay[currDay][idxY] - deliveryPerDay[currDay][idxNodeU] - deliveryPerDay[currDay][idxX])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[currDay][idxNodeU] + deliveryPerDay[currDay][idxX] - deliveryPerDay[currDay][idxNodeV] - deliveryPerDay[currDay][idxY])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( x->isADepot || y->isADepot || y == nodeUPrev || nodeU == y || x == nodeV || nodeV == nodeXNext ) { return 0 ;}

	// mettre a jour les noeuds
	swapNode(nodeU, nodeV) ;
	swapNode(x,y) ;

	stopResearch = false ; 
	return 1 ;
}
// 2-OPT
// If T(noeudU) = T(noeudV) , replace (noeudU,x) and (noeudV,y) by (noeudU,noeudV) and (x,y)
// effectue si noeudU place avant noeudV sur la meme route, et que noeudU n'est pas immediatement avant noeudV
// noeudU x et noeudV sont des sommets , y peut etre un depot.
// on n'a pas le cas ou noeudUCour est un depot.
bool LocalSearch::mutation7 ()
{
	Node * nodeNum = nodeXNext ;
	Node * temp ;

	if  ((routeU->idx != routeV->idx) || nodeU->next == nodeV || nodeU->place > nodeV->place ) {  return 0 ; }

	double cost = params->timeCost[idxNodeU][idxNodeV] + params->timeCost[idxX][idxY]
	- params->timeCost[idxNodeU][idxX] - params->timeCost[idxNodeV][idxY] ;

	if ( cost > -0.0001 ) { return 0 ;}

	// mettre a jour les noeuds
	x->prev = nodeNum ;
	x->next = y ;

	while ( nodeNum != nodeV )
	{
		temp = nodeNum->next ;
		nodeNum->next = nodeNum->prev ;
		nodeNum->prev = temp ;
		nodeNum = temp ;
	}

	nodeV->next = nodeV->prev ;
	nodeV->prev = nodeU ;
	nodeU->next = nodeV ;
	y->prev = x ;

	// et mettre a jour les routes
	routeU->updateRouteData();

	stopResearch = false ; 
	return 1 ;
}

// 2-OPT* (avec inversion du sens de parties de routes)
// If T(noeudU) != T(noeudV) , replace (noeudU,x) and (noeudV,y) by (noeudU,noeudV) and (x,y)
bool LocalSearch::mutation8 ()
{
	// TODO : heterogenous fleet, 2 types de mutations suivant les camions choisis pour chaque segment
	if  ( routeU->idx == routeV->idx || routeU->depot->idx != routeV->depot->idx) { return 0 ; }

	double chargeResteU = routeU->load - nodeU->previousLoad ;
	double chargeResteV = routeV->load - nodeV->previousLoad ;

	double cost = params->timeCost[idxNodeU][idxNodeV] 
	+ params->timeCost[idxX][idxY]
	- params->timeCost[idxNodeU][idxX] 
	- params->timeCost[idxNodeV][idxY]
    + routeU->excedentCharge(nodeU->previousLoad + nodeV->previousLoad)*params->penalityCapa[idxScenario]
	+ routeV->excedentCharge(chargeResteV + chargeResteU)*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;

	if ( cost > -0.0001 ) { return 0 ; } 

	/////////////////////////// ON EFFECTUE LA MUTATION ///////////////////////////////

	Node * depotU = routeU->depot ;
	Node * depotV = routeV->depot ;
	Node * depotUFin = routeU->depot->prev ;
	Node * depotVFin = routeV->depot->prev ;
	Node * depotVSuiv = depotV->next ;

	// on inverse le sens et on change le nom des routes
	Node * temp ;
	Node * xx = x ;
	Node * vv = nodeV ;

	while ( !xx->isADepot )
	{
		temp = xx->next ;
		xx->next = xx->prev ;
		xx->prev = temp ;
		xx->route = routeV ;
		xx = temp ;
	}

	while ( !vv->isADepot )
	{
		temp = vv->prev ;
		vv->prev = vv->next ;
		vv->next = temp ;
		vv->route = routeU ;
		vv = temp ;
	}

	// mettre a jour les noeuds
	nodeU->next = nodeV ;
	nodeV->prev = nodeU ;
	x->next = y ;
	y->prev = x ;

	// mettre � jour les extr�mit�s
	if (x->isADepot)
	{
		depotUFin->next = depotU ;
		depotUFin->prev = depotVSuiv ;
		depotUFin->prev->next = depotUFin ;
		depotV->next = y ;
		y->prev = depotV ;
	}
	else if ( nodeV->isADepot )
	{
		depotV->next = depotUFin->prev ;
		depotV->next->prev = depotV ;
		depotV->prev = depotVFin ;
		depotUFin->prev = nodeU ;
		nodeU->next = depotUFin ;
	}
	else
	{
		depotV->next = depotUFin->prev ;
		depotV->next->prev = depotV ;
		depotUFin->prev = depotVSuiv ;
		depotUFin->prev->next = depotUFin ;
	}

	// et mettre a jour les routes
	routeU->updateRouteData();
	routeV->updateRouteData();

	stopResearch = false ; 
	return 1 ;
}

// 2-OPT* (sans inversion du sens de parties de routes)
// If T(noeudU) != T(noeudV) , replace (noeudU,x) and (noeudV,y) by (noeudU,y) and (noeudV,x)
bool LocalSearch::mutation9 ()
{
	if  (routeU->idx == routeV->idx || routeU->depot->idx != routeV->depot->idx) { return 0 ; }

	Node * count ;

	double chargeResteU = routeU->load - nodeU->previousLoad ;
	double chargeResteV = routeV->load - nodeV->previousLoad ;

	double cost = params->timeCost[idxNodeU][idxY] 
	+ params->timeCost[idxNodeV][idxX]
	- params->timeCost[idxNodeU][idxX] 
	- params->timeCost[idxNodeV][idxY]
	+ (routeU->excedentCharge(nodeU->previousLoad + chargeResteV)
	+ routeV->excedentCharge(nodeV->previousLoad + chargeResteU)
	- routeU->excedentCharge(routeU->load)
	- routeV->excedentCharge(routeV->load))*params->penalityCapa[idxScenario] ;

	if ( cost > -0.0001 ) { return 0 ; } 


	/////////////////////////// ON EFFECTUE LA MUTATION ///////////////////////////////
	// on parcourt les noeuds pour les associer aux bonnes routes

	Node * depotU = routeU->depot ;
	Node * depotV = routeV->depot ;
	Node * depotUFin = depotU->prev ;
	Node * depotVFin = depotV->prev ;
	Node * depotUpred = depotUFin->prev ;

	count = y ;
	while ( !count->isADepot )
	{
		count->route = routeU ;
		count = count->next ;
	}

	count = x ;
	while ( !count->isADepot )
	{
		count->route = routeV ;
		count = count->next ;
	}

	// mettre a jour les noeuds
	nodeU->next = y ;
	y->prev = nodeU ;
	nodeV->next = x ;
	x->prev = nodeV ;

	// mettre � jour les extr�mit�s
	if (x->isADepot)
	{
		depotUFin->prev = depotVFin->prev ;
		depotUFin->prev->next = depotUFin ;
		nodeV->next = depotVFin ;
		depotVFin->prev = nodeV ;
	}
	else
	{
		depotUFin->prev = depotVFin->prev ;
		depotUFin->prev->next = depotUFin ;
		depotVFin->prev = depotUpred ;
		depotVFin->prev->next = depotVFin ;
	}

	routeU->updateRouteData();
	routeV->updateRouteData();

	stopResearch = false ;
	return 1 ;
}
