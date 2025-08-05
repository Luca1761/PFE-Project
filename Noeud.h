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

class Node
{

public :

  
bool estUnDepot ;
// est un depot ou un client

unsigned int idx ;
// indice du depot ou du client represente

int place ;
// place dans la route


unsigned int jour ;
// indice du jour en question

bool estPresent ;
// presence de le ce client � ce jour ci


Node * suiv ;
// depot ou client suivant dans la route


Node * pred ;


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
Node * placeInsertion ;
// noeud o� serait ins�r�

vector <unsigned int> moves ;
// mouvements possibles
Node(void);
// constructeur 1	

Node(bool _estUnDepot, unsigned int _idx, unsigned int _jour, bool _estPresent, Node * _suiv , Node * _pred, Route * _route);
// constructeur 2

~Node(void);
// destructeur
};

#endif
