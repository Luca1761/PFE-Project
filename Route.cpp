#include "Route.h"
#include "LocalSearch.h"

Route::Route(void){}
Route::Route(Params* _params, LocalSearch* _myLS, int _idx, int _day, Noeud * _depot, double _temps, double _charge, double _capacity) : 
params(_params), myLS(_myLS), idx(_idx), day(_day), depot(_depot), temps(_temps) , charge(_charge), capacity(_capacity)
{
	bestInsertion = vector<Insertion>(params->nbClients + params->nbDepots);
	nodeAndRouteTested = vector<bool>(params->nbClients + params->nbDepots, false);
}

Route::~Route(void){}

void Route::printRouteData(std::ostream& file)
{
	bool firstIt = true ;
	int place = 0 ;
	double charge = 0 ;
	double earlT = 0 ;

	// on parcourt du debut a la fin
	Noeud * noeud = depot ;
	noeud->place = place ;
	depot->chargeAvant = 0 ;
	depot->est = 0 ;

	while ( !noeud->estUnDepot || firstIt )
	{
		firstIt = false ;
		file <<" node[ "<<noeud->idx <<" ] ->";
		noeud = noeud->suiv ;
		place ++ ;
		noeud->place = place ;
		charge += myLS->deliveryPerDay[day][noeud->idx];
		earlT += params->timeCost[noeud->pred->idx][noeud->idx] ;
		noeud->chargeAvant = charge ;
		noeud->est = earlT ;
	}
	file <<"depot"<<endl;
	noeud->route->temps = earlT ;
	noeud->route->charge = charge ;

}

void Route::updateRouteData () {
	bool firstIt = true ;
	int place = 0 ;
	double charge = 0 ;
	double earlT = 0 ;

	// on parcourt du debut � la fin
	Noeud * noeud = depot ;
	noeud->place = place ;
	depot->chargeAvant = 0 ;
	depot->est = 0 ;

	while (!noeud->estUnDepot || firstIt) {
		firstIt = false ;
		noeud = noeud->suiv ;
		place ++ ;
		noeud->place = place ;
		charge += myLS->deliveryPerDay[day][noeud->idx];
		earlT += params->timeCost[noeud->pred->idx][noeud->idx] ;
		noeud->chargeAvant = charge ;
		noeud->est = earlT ;
	}

	noeud->route->temps = earlT ;
	noeud->route->charge = charge ;

	if (charge < capacity + 0.0001) isFeasible = true ;
	else isFeasible = false;

	initiateInsertions();
}
void Route::evalInsertClientp (Noeud * U) 
{
	Noeud * courNoeud ;
	double cout;
	bestInsertion[U->idx].detour = 1.e30 ;
	bestInsertion[U->idx].place = NULL ;
	bestInsertion[U->idx].load = -1.e30 ;
	
	bool firstIt = true ;
	if ( U->route != this || !U->estPresent )
	{
		bestInsertion[U->idx].load = std::max<double>(0.,capacity - charge);
		courNoeud = depot ;
		while (!courNoeud->estUnDepot || firstIt == true )
		{
			firstIt = false ;
			cout = params->timeCost[courNoeud->idx][U->idx] 
			+ params->timeCost[U->idx][courNoeud->suiv->idx] 
			- params->timeCost[courNoeud->idx][courNoeud->suiv->idx] ;

			if ( cout < bestInsertion[U->idx].detour - 0.0001 )
			{ 
				bestInsertion[U->idx].detour = cout ;
				bestInsertion[U->idx].place = courNoeud ;
			}
			courNoeud = courNoeud->suiv ;
		}
	}
	else // U is already  in our route
	{
		bestInsertion[U->idx].load = std::max<double>(0.,capacity + myLS->deliveryPerDay[day][U->idx] - charge);
		bestInsertion[U->idx].detour = params->timeCost[U->pred->idx][U->idx] - params->timeCost[U->pred->idx][U->suiv->idx]   
										+ params->timeCost[U->idx][U->suiv->idx] ;
		bestInsertion[U->idx].place = U->pred ;
		std::cout << "else detour "<<bestInsertion[U->idx].detour<< " pre "<<U->pred->idx <<" next "<<U->suiv->idx<<endl;

		// however, we'll see if there's a better insertion possible
		// temporarily we'll remove the node from the chain�e list (in O(1))
		U->pred->suiv = U->suiv ;
		U->suiv->pred = U->pred ;
		courNoeud = depot ;

		// et parcourir la route � nouveau
		while (!courNoeud->estUnDepot || firstIt == true )
		{
			firstIt = false ;
			cout = params->timeCost[courNoeud->idx][U->idx] 
			+ params->timeCost[U->idx][courNoeud->suiv->idx]
			- params->timeCost[courNoeud->idx][courNoeud->suiv->idx] ;

			// au final on peut placer d'une meilleure mani�re
			if (cout < bestInsertion[U->idx].detour - 0.0001)
			{ 
				bestInsertion[U->idx].detour = cout ;
				bestInsertion[U->idx].place = courNoeud ;
				std::cout << "better "<<bestInsertion[U->idx].detour;
			}
			courNoeud = courNoeud->suiv ;
		}
		
		// on replace le noeud
		U->pred->suiv = U ;
		U->suiv->pred = U ;
	}
}


// pour un client donn�, trouve la meilleure position d'insertion et sa quantit�
void Route::evalInsertClient (Noeud * U) 
{
	
	Noeud * courNoeud ;
	double cout;
	bestInsertion[U->idx].detour = 1.e30 ;
	bestInsertion[U->idx].place = NULL ;
	bestInsertion[U->idx].load = -1.e30 ;
	
	bool firstIt = true ;
	if ( U->route != this || !U->estPresent )
	{
		bestInsertion[U->idx].load = std::max<double>(0.,capacity - charge);
		courNoeud = depot ;
		while (!courNoeud->estUnDepot || firstIt == true )
		{
			firstIt = false ;
			cout = params->timeCost[courNoeud->idx][U->idx] 
			+ params->timeCost[U->idx][courNoeud->suiv->idx] 
			- params->timeCost[courNoeud->idx][courNoeud->suiv->idx] ;
			
			
			if ( cout < bestInsertion[U->idx].detour - 0.0001 )
			{ 
				bestInsertion[U->idx].detour = cout ;
				bestInsertion[U->idx].place = courNoeud ;
			}
			courNoeud = courNoeud->suiv ;
		}
	}
	else // U is already  in our route
	{
		bestInsertion[U->idx].load = std::max<double>(0.,capacity + myLS->deliveryPerDay[day][U->idx] - charge);
		bestInsertion[U->idx].detour = params->timeCost[U->pred->idx][U->idx] - params->timeCost[U->pred->idx][U->suiv->idx]   
										+ params->timeCost[U->idx][U->suiv->idx] ;
		bestInsertion[U->idx].place = U->pred ;

		// however, we'll see if there's a better insertion possible
		// temporarily we'll remove the node from the chain�e list (in O(1))
		U->pred->suiv = U->suiv ;
		U->suiv->pred = U->pred ;
		courNoeud = depot ;

		// et parcourir la route � nouveau
		while (!courNoeud->estUnDepot || firstIt == true )
		{
			firstIt = false ;
			cout = params->timeCost[courNoeud->idx][U->idx] 
			+ params->timeCost[U->idx][courNoeud->suiv->idx]
			- params->timeCost[courNoeud->idx][courNoeud->suiv->idx] ;

			// au final on peut placer d'une meilleure mani�re
			if (cout < bestInsertion[U->idx].detour - 0.0001)
			{ 
				bestInsertion[U->idx].detour = cout ;
				bestInsertion[U->idx].place = courNoeud ;
			}
			courNoeud = courNoeud->suiv ;
		}
		
		// on replace le noeud
		U->pred->suiv = U ;
		U->suiv->pred = U ;
	}
}

// no insertion is calculated
void Route::initiateInsertions() {
	for (unsigned int i = 0 ; i < params->nbClients + params->nbDepots ; i++) {
		bestInsertion[i].detour = 1.e30;
		bestInsertion[i].load = -1.e30;
		bestInsertion[i].place = NULL;
	}
}

// moves having nodes in this route need to be examined again
void Route::reinitSingleDayMoves()
{
	for (unsigned int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
		nodeAndRouteTested[i] = false ;
}
