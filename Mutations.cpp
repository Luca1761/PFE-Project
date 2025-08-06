#include "LocalSearch.h"

// effectue l'insertion du client U apres V
void LocalSearch::insertNoeud(Node * U, Node * V) {
	// mettre a jour les noeuds
	U->prev->next = U->next ;
	U->next->prev = U->prev ;
	V->next->prev = U ;
	U->prev = V ;
	U->next = V->next ;
	V->next = U ;

	U->route = routeV ;
	routeU->updateRouteData();
	routeV->updateRouteData();
}
// effectue le swap du client U avec V
void LocalSearch::swapNoeud(Node * U, Node * V) 
{
	Node * VPred = V->prev ;
	Node * VSuiv = V->next ;
	Node * UPred = U->prev ;
	Node * USuiv = U->next ;
	Route * myRouteU = U->route ;
	Route * myRouteV = V->route ;

	UPred->next = V ;
	USuiv->prev = V ;
	VPred->next = U ;
	VSuiv->prev = U ;

	U->prev = VPred ;
	U->next = VSuiv ;
	V->prev = UPred ;
	V->next = USuiv ;

	U->route = myRouteV ;
	V->route = myRouteU ;
	U->route->updateRouteData();
	V->route->updateRouteData();
}
// INSERT
// If noeudU is a client node, remove noeudU then insert it after noeudV
int LocalSearch::mutation1 ()
{
	double costSuppU = params->timeCost[noeudUPredCour][xCour] 
	- params->timeCost[noeudUPredCour][noeudUCour]  
	- params->timeCost[noeudUCour][xCour];

	double costSuppV = params->timeCost[noeudVCour][noeudUCour] 
	+ params->timeCost[noeudUCour][yCour] 
	- params->timeCost[noeudVCour][yCour];

	// dans le cas ou l'on est dans la meme route , le cout n'est pas calcule correctement en realite
	// tout ce qu'on sait c'est que si il est negatif c'est qu'il est bien reellement negatif
	// pas d'incidence pour l'instant mais attention
	if (routeU != routeV) {
		costSuppU += routeU->excedentCharge(routeU->load - deliveryPerDay[dayCour][noeudUCour])*params->penalityCapa[idxScenario]
		- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;

		costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[dayCour][noeudUCour])*params->penalityCapa[idxScenario]
		- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if (costSuppU + costSuppV > -0.0001) return 0 ;
	if (noeudUCour == yCour) return 0;

	// mettre a jour les noeuds
	insertNoeud(noeudU,noeudV);

	rechercheTerminee = false ; 
	return 1 ;
}

// If noeudU and x are clients, remove them then insert (noeudU,x) after noeudV
// teste si x n'est pas un depot , et si x different de noeudV, et si noeudU pas deja apres noeudV
int LocalSearch::mutation2 ()
{
	double costSuppU = params->timeCost[noeudUPredCour][xSuivCour] 
	- params->timeCost[noeudUPredCour][noeudUCour] 
	- params->timeCost[noeudUCour][xCour] 
	- params->timeCost[xCour][xSuivCour];

	double costSuppV = params->timeCost[noeudVCour][noeudUCour] 
	+ params->timeCost[noeudUCour][xCour] 
	+ params->timeCost[xCour][yCour] 
	- params->timeCost[noeudVCour][yCour];

	if (routeU != routeV) {
	costSuppU += routeU->excedentCharge(routeU->load - deliveryPerDay[dayCour][noeudUCour] - deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[dayCour][noeudUCour] + deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( noeudU == y || noeudV == x || x->isADepot ) { return 0 ;}

	// mettre a jour les noeuds
	insertNoeud(noeudU,noeudV);
	insertNoeud(x,noeudU);

	rechercheTerminee = false ; 
	return 1 ;
}

// If noeudU and x are clients, remove them then insert (x,noeudU) after noeudV
// teste si x n'est pas un depot , et si x different de noeudV, et si noeudU pas d�ja apres noeudV
int LocalSearch::mutation3 ()
{
	double costSuppU = params->timeCost[noeudUPredCour][xSuivCour] 
	- params->timeCost[noeudUPredCour][noeudUCour] 
	- params->timeCost[noeudUCour][xCour] 
	- params->timeCost[xCour][xSuivCour];

	double costSuppV = params->timeCost[noeudVCour][xCour] 
	+ params->timeCost[xCour][noeudUCour] 
	+ params->timeCost[noeudUCour][yCour] 
	- params->timeCost[noeudVCour][yCour];

	if (routeU != routeV) {
	costSuppU += routeU->excedentCharge(routeU->load - deliveryPerDay[dayCour][noeudUCour] - deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[dayCour][noeudUCour] + deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) return 0;
	if ( noeudU == y ||  x == noeudV || x->isADepot ) return 0;

	// mettre a jour les noeuds
	insertNoeud(x,noeudV);
	insertNoeud(noeudU,x);

	rechercheTerminee = false ; 
	return 1 ;
}
// SWAP
// If noeudU and noeudV are clients, swap noeudU and noeudV
// sauf si noeudU et noeudV se succedent
int LocalSearch::mutation4 ()
{
	double costSuppU = params->timeCost[noeudUPredCour][noeudVCour] 
	+ params->timeCost[noeudVCour][xCour]
	- params->timeCost[noeudUPredCour][noeudUCour] 
	- params->timeCost[noeudUCour][xCour];

	double costSuppV = params->timeCost[noeudVPredCour][noeudUCour] 
	+ params->timeCost[noeudUCour][yCour]
	- params->timeCost[noeudVPredCour][noeudVCour] 
	- params->timeCost[noeudVCour][yCour];

	if (routeU != routeV) {
	costSuppU += routeU->excedentCharge(routeU->load + deliveryPerDay[dayCour][noeudVCour] - deliveryPerDay[dayCour][noeudUCour])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[dayCour][noeudUCour] - deliveryPerDay[dayCour][noeudVCour])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( noeudUCour == noeudVPredCour || noeudUCour == yCour) { return 0 ;}

	// mettre a jour les noeuds
	swapNoeud(noeudU, noeudV) ;

	rechercheTerminee = false ; 
	return 1 ;
}

// If noeudU, x and noeudV are clients, swap (noeudU,x) and noeudV
int LocalSearch::mutation5 ()
{
	// on ne fait pas le cas ou x et noeudVCour se suivent
	// car il faut traiter autrement
	// et la mutation 2 entre (noeudUCour,x) et noeudVCour fait le meme travail correctement

	double costSuppU = params->timeCost[noeudUPredCour][noeudVCour] 
	+ params->timeCost[noeudVCour][xSuivCour]
	- params->timeCost[noeudUPredCour][noeudUCour] 
	- params->timeCost[noeudUCour][xCour] 
	- params->timeCost[xCour][xSuivCour];

	double costSuppV = params->timeCost[noeudVPredCour][noeudUCour] 
	+ params->timeCost[xCour][yCour]
	+ params->timeCost[noeudUCour][xCour]
	- params->timeCost[noeudVPred->idx][noeudVCour] 
	- params->timeCost[noeudVCour][yCour];
	
	if (routeU != routeV) {
	costSuppU += routeU->excedentCharge(routeU->load + deliveryPerDay[dayCour][noeudVCour] - deliveryPerDay[dayCour][noeudUCour] - deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;

	costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[dayCour][noeudUCour] + deliveryPerDay[dayCour][xCour] - deliveryPerDay[dayCour][noeudVCour])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( noeudU == noeudVPred || x == noeudVPred || noeudU == y || x->isADepot ) { return 0 ;}

	// mettre a jour les noeuds
	swapNoeud(noeudU, noeudV) ;
	insertNoeud(x,noeudU);

	rechercheTerminee = false ; 
	return 1 ;
}
// If (noeudU,x) and (noeudV,y) are clients, swap (noeudU,x) and (noeudV,y)
int LocalSearch::mutation6 ()
{
	double costSuppU = params->timeCost[noeudUPredCour][noeudVCour]  
	+ params->timeCost[noeudVCour][yCour]
	+ params->timeCost[yCour][xSuivCour]
	- params->timeCost[noeudUPredCour][noeudUCour] 
	- params->timeCost[noeudUCour][xCour] 
	- params->timeCost[xCour][xSuivCour];

	double costSuppV = params->timeCost[noeudVPredCour][noeudUCour] 
	+ params->timeCost[noeudUCour][xCour]
	+ params->timeCost[xCour][ySuivCour]
	- params->timeCost[noeudVPredCour][noeudVCour] 
	- params->timeCost[noeudVCour][yCour]
	- params->timeCost[yCour][ySuivCour];
	
	if (routeU != routeV) {
	costSuppU += routeU->excedentCharge(routeU->load + deliveryPerDay[dayCour][noeudVCour] + deliveryPerDay[dayCour][yCour] - deliveryPerDay[dayCour][noeudUCour] - deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->load)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->load + deliveryPerDay[dayCour][noeudUCour] + deliveryPerDay[dayCour][xCour] - deliveryPerDay[dayCour][noeudVCour] - deliveryPerDay[dayCour][yCour])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->load)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( x->isADepot || y->isADepot || y == noeudUPred || noeudU == y || x == noeudV || noeudV == noeudXSuiv ) { return 0 ;}

	// mettre a jour les noeuds
	swapNoeud(noeudU, noeudV) ;
	swapNoeud(x,y) ;

	rechercheTerminee = false ; 
	return 1 ;
}
// 2-OPT
// If T(noeudU) = T(noeudV) , replace (noeudU,x) and (noeudV,y) by (noeudU,noeudV) and (x,y)
// effectue si noeudU place avant noeudV sur la meme route, et que noeudU n'est pas immediatement avant noeudV
// noeudU x et noeudV sont des sommets , y peut etre un depot.
// on n'a pas le cas ou noeudUCour est un depot.
int LocalSearch::mutation7 ()
{
	Node * nodeNum = noeudXSuiv ;
	Node * temp ;

	if  ((routeU->idx != routeV->idx) || noeudU->next == noeudV || noeudU->place > noeudV->place ) {  return 0 ; }

	double cost = params->timeCost[noeudUCour][noeudVCour] + params->timeCost[xCour][yCour]
	- params->timeCost[noeudUCour][xCour] - params->timeCost[noeudVCour][yCour] ;

	if ( cost > -0.0001 ) { return 0 ;}

	// mettre a jour les noeuds
	x->prev = nodeNum ;
	x->next = y ;

	while ( nodeNum != noeudV )
	{
		temp = nodeNum->next ;
		nodeNum->next = nodeNum->prev ;
		nodeNum->prev = temp ;
		nodeNum = temp ;
	}

	noeudV->next = noeudV->prev ;
	noeudV->prev = noeudU ;
	noeudU->next = noeudV ;
	y->prev = x ;

	// et mettre a jour les routes
	routeU->updateRouteData();

	rechercheTerminee = false ; 
	return 1 ;
}

// 2-OPT* (avec inversion du sens de parties de routes)
// If T(noeudU) != T(noeudV) , replace (noeudU,x) and (noeudV,y) by (noeudU,noeudV) and (x,y)
int LocalSearch::mutation8 ()
{
	// TODO : heterogenous fleet, 2 types de mutations suivant les camions choisis pour chaque segment
	if  ( routeU->idx == routeV->idx || routeU->depot->idx != routeV->depot->idx) { return 0 ; }

	double chargeResteU = routeU->load - noeudU->previousLoad ;
	double chargeResteV = routeV->load - noeudV->previousLoad ;

	double cost = params->timeCost[noeudUCour][noeudVCour] 
	+ params->timeCost[xCour][yCour]
	- params->timeCost[noeudUCour][xCour] 
	- params->timeCost[noeudVCour][yCour]
    + routeU->excedentCharge(noeudU->previousLoad + noeudV->previousLoad)*params->penalityCapa[idxScenario]
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
	Node * vv = noeudV ;

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
	noeudU->next = noeudV ;
	noeudV->prev = noeudU ;
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
	else if ( noeudV->isADepot )
	{
		depotV->next = depotUFin->prev ;
		depotV->next->prev = depotV ;
		depotV->prev = depotVFin ;
		depotUFin->prev = noeudU ;
		noeudU->next = depotUFin ;
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

	rechercheTerminee = false ; 
	return 1 ;
}

// 2-OPT* (sans inversion du sens de parties de routes)
// If T(noeudU) != T(noeudV) , replace (noeudU,x) and (noeudV,y) by (noeudU,y) and (noeudV,x)
int LocalSearch::mutation9 ()
{
	if  (routeU->idx == routeV->idx || routeU->depot->idx != routeV->depot->idx) { return 0 ; }

	Node * count ;

	double chargeResteU = routeU->load - noeudU->previousLoad ;
	double chargeResteV = routeV->load - noeudV->previousLoad ;

	double cost = params->timeCost[noeudUCour][yCour] 
	+ params->timeCost[noeudVCour][xCour]
	- params->timeCost[noeudUCour][xCour] 
	- params->timeCost[noeudVCour][yCour]
	+ (routeU->excedentCharge(noeudU->previousLoad + chargeResteV)
	+ routeV->excedentCharge(noeudV->previousLoad + chargeResteU)
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
	noeudU->next = y ;
	y->prev = noeudU ;
	noeudV->next = x ;
	x->prev = noeudV ;

	// mettre � jour les extr�mit�s
	if (x->isADepot)
	{
		depotUFin->prev = depotVFin->prev ;
		depotUFin->prev->next = depotUFin ;
		noeudV->next = depotVFin ;
		depotVFin->prev = noeudV ;
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

	rechercheTerminee = false ;
	return 1 ;
}
