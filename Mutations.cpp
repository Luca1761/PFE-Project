#include "LocalSearch.h"

// effectue l'insertion du client U apres V
void LocalSearch::insertNoeud(Noeud * U, Noeud * V) {
	// mettre a jour les noeuds
	U->pred->suiv = U->suiv ;
	U->suiv->pred = U->pred ;
	V->suiv->pred = U ;
	U->pred = V ;
	U->suiv = V->suiv ;
	V->suiv = U ;

	U->route = routeV ;
	routeU->updateRouteData();
	routeV->updateRouteData();
}
// effectue le swap du client U avec V
void LocalSearch::swapNoeud(Noeud * U, Noeud * V) 
{
	Noeud * VPred = V->pred ;
	Noeud * VSuiv = V->suiv ;
	Noeud * UPred = U->pred ;
	Noeud * USuiv = U->suiv ;
	Route * myRouteU = U->route ;
	Route * myRouteV = V->route ;

	UPred->suiv = V ;
	USuiv->pred = V ;
	VPred->suiv = U ;
	VSuiv->pred = U ;

	U->pred = VPred ;
	U->suiv = VSuiv ;
	V->pred = UPred ;
	V->suiv = USuiv ;

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
		costSuppU += routeU->excedentCharge(routeU->charge - deliveryPerDay[dayCour][noeudUCour])*params->penalityCapa[idxScenario]
		- routeU->excedentCharge(routeU->charge)*params->penalityCapa[idxScenario] ;

		costSuppV += routeV->excedentCharge(routeV->charge + deliveryPerDay[dayCour][noeudUCour])*params->penalityCapa[idxScenario]
		- routeV->excedentCharge(routeV->charge)*params->penalityCapa[idxScenario] ;
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

	if (routeU != routeV)
	{
	costSuppU += routeU->excedentCharge(routeU->charge - deliveryPerDay[dayCour][noeudUCour] - deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->charge)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->charge + deliveryPerDay[dayCour][noeudUCour] + deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->charge)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( noeudU == y || noeudV == x || x->estUnDepot ) { return 0 ;}

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

	if (routeU != routeV)
	{
	costSuppU += routeU->excedentCharge(routeU->charge - deliveryPerDay[dayCour][noeudUCour] - deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->charge)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->charge + deliveryPerDay[dayCour][noeudUCour] + deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->charge)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) return 0;
	if ( noeudU == y ||  x == noeudV || x->estUnDepot ) return 0;

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

	if (routeU != routeV)
	{
	costSuppU += routeU->excedentCharge(routeU->charge + deliveryPerDay[dayCour][noeudVCour] - deliveryPerDay[dayCour][noeudUCour])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->charge)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->charge + deliveryPerDay[dayCour][noeudUCour] - deliveryPerDay[dayCour][noeudVCour])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->charge)*params->penalityCapa[idxScenario] ;
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
	
	if (routeU != routeV)
	{
	costSuppU += routeU->excedentCharge(routeU->charge + deliveryPerDay[dayCour][noeudVCour] - deliveryPerDay[dayCour][noeudUCour] - deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->charge)*params->penalityCapa[idxScenario] ;

	costSuppV += routeV->excedentCharge(routeV->charge + deliveryPerDay[dayCour][noeudUCour] + deliveryPerDay[dayCour][xCour] - deliveryPerDay[dayCour][noeudVCour])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->charge)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( noeudU == noeudVPred || x == noeudVPred || noeudU == y || x->estUnDepot ) { return 0 ;}

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
	
	if (routeU != routeV)
	{
	costSuppU += routeU->excedentCharge(routeU->charge + deliveryPerDay[dayCour][noeudVCour] + deliveryPerDay[dayCour][yCour] - deliveryPerDay[dayCour][noeudUCour] - deliveryPerDay[dayCour][xCour])*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->charge)*params->penalityCapa[idxScenario] ;
	
	costSuppV += routeV->excedentCharge(routeV->charge + deliveryPerDay[dayCour][noeudUCour] + deliveryPerDay[dayCour][xCour] - deliveryPerDay[dayCour][noeudVCour] - deliveryPerDay[dayCour][yCour])*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->charge)*params->penalityCapa[idxScenario] ;
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( x->estUnDepot || y->estUnDepot || y == noeudUPred || noeudU == y || x == noeudV || noeudV == noeudXSuiv ) { return 0 ;}

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
	Noeud * nodeNum = noeudXSuiv ;
	Noeud * temp ;

	if  ((routeU->idx != routeV->idx) || noeudU->suiv == noeudV || noeudU->place > noeudV->place ) {  return 0 ; }

	double cost = params->timeCost[noeudUCour][noeudVCour] + params->timeCost[xCour][yCour]
	- params->timeCost[noeudUCour][xCour] - params->timeCost[noeudVCour][yCour] ;

	if ( cost > -0.0001 ) { return 0 ;}

	// mettre a jour les noeuds
	x->pred = nodeNum ;
	x->suiv = y ;

	while ( nodeNum != noeudV )
	{
		temp = nodeNum->suiv ;
		nodeNum->suiv = nodeNum->pred ;
		nodeNum->pred = temp ;
		nodeNum = temp ;
	}

	noeudV->suiv = noeudV->pred ;
	noeudV->pred = noeudU ;
	noeudU->suiv = noeudV ;
	y->pred = x ;

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

	double chargeResteU = routeU->charge - noeudU->chargeAvant ;
	double chargeResteV = routeV->charge - noeudV->chargeAvant ;

	double cost = params->timeCost[noeudUCour][noeudVCour] 
	+ params->timeCost[xCour][yCour]
	- params->timeCost[noeudUCour][xCour] 
	- params->timeCost[noeudVCour][yCour]
    + routeU->excedentCharge(noeudU->chargeAvant + noeudV->chargeAvant)*params->penalityCapa[idxScenario]
	+ routeV->excedentCharge(chargeResteV + chargeResteU)*params->penalityCapa[idxScenario]
	- routeU->excedentCharge(routeU->charge)*params->penalityCapa[idxScenario]
	- routeV->excedentCharge(routeV->charge)*params->penalityCapa[idxScenario] ;

	if ( cost > -0.0001 ) { return 0 ; } 

	/////////////////////////// ON EFFECTUE LA MUTATION ///////////////////////////////

	Noeud * depotU = routeU->depot ;
	Noeud * depotV = routeV->depot ;
	Noeud * depotUFin = routeU->depot->pred ;
	Noeud * depotVFin = routeV->depot->pred ;
	Noeud * depotVSuiv = depotV->suiv ;
	Noeud * depotUSuiv = depotU->suiv ;
	Noeud * depotVPred = depotVFin->pred ;

	// on inverse le sens et on change le nom des routes
	Noeud * temp ;
	Noeud * xx = x ;
	Noeud * vv = noeudV ;

	while ( !xx->estUnDepot )
	{
		temp = xx->suiv ;
		xx->suiv = xx->pred ;
		xx->pred = temp ;
		xx->route = routeV ;
		xx = temp ;
	}

	while ( !vv->estUnDepot )
	{
		temp = vv->pred ;
		vv->pred = vv->suiv ;
		vv->suiv = temp ;
		vv->route = routeU ;
		vv = temp ;
	}

	// mettre a jour les noeuds
	noeudU->suiv = noeudV ;
	noeudV->pred = noeudU ;
	x->suiv = y ;
	y->pred = x ;

	// mettre � jour les extr�mit�s
	if (x->estUnDepot)
	{
		depotUFin->suiv = depotU ;
		depotUFin->pred = depotVSuiv ;
		depotUFin->pred->suiv = depotUFin ;
		depotV->suiv = y ;
		y->pred = depotV ;
	}
	else if ( noeudV->estUnDepot )
	{
		depotV->suiv = depotUFin->pred ;
		depotV->suiv->pred = depotV ;
		depotV->pred = depotVFin ;
		depotUFin->pred = noeudU ;
		noeudU->suiv = depotUFin ;
	}
	else
	{
		depotV->suiv = depotUFin->pred ;
		depotV->suiv->pred = depotV ;
		depotUFin->pred = depotVSuiv ;
		depotUFin->pred->suiv = depotUFin ;
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

	Noeud * count ;

	double chargeResteU = routeU->charge - noeudU->chargeAvant ;
	double chargeResteV = routeV->charge - noeudV->chargeAvant ;

	double cost = params->timeCost[noeudUCour][yCour] 
	+ params->timeCost[noeudVCour][xCour]
	- params->timeCost[noeudUCour][xCour] 
	- params->timeCost[noeudVCour][yCour]
	+ (routeU->excedentCharge(noeudU->chargeAvant + chargeResteV)
	+ routeV->excedentCharge(noeudV->chargeAvant + chargeResteU)
	- routeU->excedentCharge(routeU->charge)
	- routeV->excedentCharge(routeV->charge))*params->penalityCapa[idxScenario] ;

	if ( cost > -0.0001 ) { return 0 ; } 


	/////////////////////////// ON EFFECTUE LA MUTATION ///////////////////////////////
	// on parcourt les noeuds pour les associer aux bonnes routes

	Noeud * depotU = routeU->depot ;
	Noeud * depotV = routeV->depot ;
	Noeud * depotUFin = depotU->pred ;
	Noeud * depotVFin = depotV->pred ;
	Noeud * depotUpred = depotUFin->pred ;
	Noeud * depotUSuiv = depotU->suiv ;
	Noeud * depotVPred = depotVFin->pred ;

	count = y ;
	while ( !count->estUnDepot )
	{
		count->route = routeU ;
		count = count->suiv ;
	}

	count = x ;
	while ( !count->estUnDepot )
	{
		count->route = routeV ;
		count = count->suiv ;
	}

	// mettre a jour les noeuds
	noeudU->suiv = y ;
	y->pred = noeudU ;
	noeudV->suiv = x ;
	x->pred = noeudV ;

	// mettre � jour les extr�mit�s
	if (x->estUnDepot)
	{
		depotUFin->pred = depotVFin->pred ;
		depotUFin->pred->suiv = depotUFin ;
		noeudV->suiv = depotVFin ;
		depotVFin->pred = noeudV ;
	}
	else
	{
		depotUFin->pred = depotVFin->pred ;
		depotUFin->pred->suiv = depotUFin ;
		depotVFin->pred = depotUpred ;
		depotVFin->pred->suiv = depotVFin ;
	}

	routeU->updateRouteData();
	routeV->updateRouteData();

	rechercheTerminee = false ;
	return 1 ;
}
