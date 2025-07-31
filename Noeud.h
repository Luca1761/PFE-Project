/*                       Algorithme - HGSADC                         */
/*                    Propri�t� de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef NOEUD_H
#define NOEUD_H

#include "Route.h"
#include<iostream>
using namespace std;

class Route ;

struct Move {
	
int destination ;

double cout ;

};

class Noeud
{

public :

  
bool estUnDepot ;
// est un depot ou un client

int idx ;
// indice du depot ou du client represente

int place ;
// place dans la route


unsigned int jour ;
// indice du jour en question

bool estPresent ;
// presence de le ce client � ce jour ci


Noeud * suiv ;
// depot ou client suivant dans la route


Noeud * pred ;


Route * route ;
// depot ou client precedent dans la route
// route associee 

double chargeAvant ;
// charge de la portion de route situee avant lui (lui compris)

//  List of possible insertions in different routes
vector <Insertion> allInsertions ;

void removeDominatedInsertions (double penalityCapa);
// Removing dominated insertions

double coutInsertion ;
// cout insertion dans ce jour si il devait �tre ins�r�
Noeud * placeInsertion ;
// noeud o� serait ins�r�

vector < int > moves ;
// mouvements possibles
Noeud(void);
// constructeur 1	

Noeud(bool _estUnDepot, int _idx, int _jour, bool _estPresent, Noeud * _suiv , Noeud * _pred, Route * _route);
// constructeur 2

~Noeud(void);
// destructeur
};

#endif
