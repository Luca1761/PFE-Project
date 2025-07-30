/*                       Algorithme - HGSADC                         */
/*                    Propri�t� de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

// structure de donnees adaptee a la recherche locale
// une structure de ce type associee a chaque individu de la population

// nouvelle representation des structures

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <vector>
#include "Individu.h"
#include "LotSizingSolver.h"
#include "Noeud.h"
#include "Route.h"

using namespace std;

class Individu;

struct paireJours
{
  int j1;
  int j2;
};

class LocalSearch
{
private:
Individu *individu;

public:
  Params *params;

  int idxScenario;

  bool rechercheTerminee;
  // vecteur donnant l'ordre de parcours des sommets pour chaque jour, ne
  // contenant pas les sommets
  // qui n'existent pas pour le jour donn�
  // afin de diversifier la recherche
  vector<vector<int>> ordreParcours;

  // ajoute un client dans l'ordre de parcours
  void addOP(int day, int client);

  void removeOP(int day, int client);
  
  void melangeParcours();
  
  void updateMoves();
  
  bool firstLoop;

  Noeud *noeudU;
  Noeud *noeudUPred;
  Noeud *x;
  Noeud *noeudXSuiv;
  Noeud *noeudV;
  Noeud *noeudVPred;
  Noeud *y;
  Noeud *noeudYSuiv;
  Route *routeU;
  Route *routeV;
  int noeudUCour, noeudUPredCour, xCour, xSuivCour, ySuivCour, noeudVCour,
      noeudVPredCour, yCour;
  int dayCour;

  /* vecteur de taille nbClients , l'element client(day)(i) contient des donnees
  // relatives
  // a l'emplacement de la visite du client i+1 dans les routes
  clients
│
├── Day 1
│   ├── Client 1 -> Noeud Pointer
│   ├── Client 2 -> Noeud Pointer
│   └── ...
│
├── Day 2
│   ├── Client 1 -> Noeud Pointer
│   ├── Client 2 -> Noeud Pointer
│   └── ...
│
└── ...
*/
  vector<vector<Noeud*>> clients;

  // noeuds associes aux depots utilises
  vector<vector<Noeud*>> depots;

  // noeuds associes aux terminaisons des routes (doublon des depots)
  vector<vector<Noeud*>> depotsFin;

  // vecteur repertoriant des donnees sur les routes routes
  vector<vector<Route*>> routes;

  // deliveryPerDay[i][j] -> The load to be delivered to each customer [j] on day [i]
  vector<vector<double>> deliveryPerDay;

  // effectue une parcours complet de toutes les mutations possibles
  // retourne le nombre de mouvements effectues
  int mutationSameDay(unsigned int day);

  // pour un client, marque que tous les mouvements impliquant ce noeud ont �t�
  // test�s pour chaque route du jour day
  void nodeTestedForEachRoute(int cli, int day);

  // effectue un parcours complet de tous les changement de pattern et swap
  // intra-jours possibles
  // retourne le nombre de mouvements effectu�s

  void runSearchSameDay();
  // Neighborhoods

  /* RELOCATE */

  // If posUreal is a client node, remove posUreal then insert it after posVreal
  int mutation1();

  // If posUreal and x are clients, remove them then insert (posUreal,x) after
  // posVreal
  int mutation2();

  // If posUreal and x are clients, remove them then insert (x,posUreal) after
  // posVreal
  int mutation3();

  /* SWAP */

  // If posUreal and posVreal are clients, swap posUreal and posVreal
  int mutation4();

  // If posUreal, x and posVreal are clients, swap (posUreal,x) and posVreal
  int mutation5();

  // If (posUreal,x) and (posVreal,y) are cliens, swap (posUreal,x) and
  // (posVreal,y)
  int mutation6();

  /* 2-OPT and 2-OPT* */

  // If T(posUreal) = T(posVreal) , replace (posUreal,x) and (posVreal,y) by
  // (posUreal,posVreal) and (x,y)
  int mutation7();

  // If T(posUreal) != T(posVreal) , replace (posUreal,x) and (posVreal,y) by
  // (posUreal,posVreal) and (x,y)
  int mutation8();

  // If T(posUreal) != T(posVreal) , replace (posUreal,x) and (posVreal,y) by
  // (posUreal,y) and (posVreal,x)
  int mutation9();

  // Evaluates the current objective function from the model
  double evaluateCurrentCost_stockout(int client);

  // Prints some useful information on the current solution
  void printInventoryLevels(std::ostream& file,bool add);

  /* Routines to update the solution */

  // effectue l'insertion du client U apres V
  void insertNoeud(Noeud *U, Noeud *V);

  // effectue le swap du client U avec V
  void swapNoeud(Noeud *U, Noeud *V);

  // supprime le noeud
  void removeNoeud(Noeud *U);

  // ajoute un noeud � l'endroit indique dans Noeud->placeRoute
  void addNoeud(Noeud *U);

  // calcule pour un jour donn� et un client donn� (repr�sent� par un noeud)
  // les couts d'insertion dans les differentes routes constituant ce jour
  void computeCoutInsertion(Noeud *client);

  LocalSearch();

  // constructeur, cree les structures de noeuds
  // n'initialise pas pas la pile ni les routes
  LocalSearch(Individu* _individu, Params* _params, int _idxScenario);

  ~LocalSearch(void);
};

#endif
