/*                       Algorithme - HGSADC                         */
/*                    Propri�t� de Thibaut VIDAL                     */
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

class Noeud ;
class LocalSearch ;

class Route
{

private:

// acces aux donnees de l'instance
Params * params;

// access  to other features
LocalSearch * myLS ;

public:

// numero de la route
int idx ;

// day associated to the route
int day ;

// depot associe a la route
Noeud * depot ;

// distance total de parcours sur la route
double temps ;

// chargement total sur la route
double charge ;

// chargement maximum de la route
double capacity ;

// valide ou non
bool isFeasible ;

inline double excedentCharge(double _charge) {
	return std::max<double>(0, _charge - capacity);
}

// met a jour les charges partielles de la route associ�e au noeud U
void updateRouteData () ;
void printRouteData (std::ostream& file) ;
// pour chaque noeud, stocke le cout de l'insertion dans la route
vector<Insertion> bestInsertion ;

// pour chaque noeud, booleen indiquant si tous les mouvements impliquant ce noeud
// et cette route ont ete testes sans succes
vector<bool> nodeAndRouteTested ;

// pour un client donne, trouve la meilleure position d'insertion
// eventuellement fait le calcul pour tous les clients d'un coup, ou pour plusieurs
void evalInsertClient (Noeud * U) ;
void evalInsertClientp (Noeud * U) ;
// no insertion is calculated
void initiateInsertions();

// moves having nodes in this route need to be examined again:
void reinitSingleDayMoves();

Route(void);

Route(int idx, int day, Noeud * depot, double temps, double charge, double capacity, Params * params, LocalSearch * myLS);

~Route(void);
};

#endif
