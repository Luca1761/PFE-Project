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
Params * params ;

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

// duree maximum de la route
double maxRouteTime ;

// chargement maximum de la route
double capacity ;

// valide ou non
bool isFeasible ;

inline double excedentCharge(double charge)
{
	return std::max<double>(0,charge-capacity);
}

inline double excedentLength(double length)
{
	return std::max<double>(0,length-maxRouteTime);
}

// coordonnees du centroide de la route
// ainsi que l'angle pris par rapport au segment (0,0) (0,1)
// utilise pour refaire le giant tour apres la LS
// pas d'utilite lors de la recherche locale
double centroidX ;
double centroidY ;
double centroidAngle ;

// calcule les coordonnees du centroide
void updateCentroidCoord ();

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

Route(int idx, int day, Noeud * depot, double temps, double charge, double maxRouteTime, double capacity, Params * params, LocalSearch * myLS);

~Route(void);
};

#endif
