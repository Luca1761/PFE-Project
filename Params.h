/*                       Algorithme - HGSADC                         */
/*                    Propri�t� de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef PARAMS_H
#define PARAMS_H

#include <math.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <random>
#include "Client.h"
#include "Rng.h"
#include "Vehicle.h"
using namespace std;

const double MAXCOST = 1.e30;

class Vehicle;
class Noeud;

// needed structure for a few places in the code (easily accessible from here)
struct Insertion
{
       double detour;
       //the remain load of this route
       double load;

       Noeud *place;

       Insertion()
       {
              detour = 1.e30;
              load = -1.e30;
              place = NULL;
       }
       Insertion(double _detour, double _load, Noeud *_place)
           : detour(_detour), load(_load), place(_place) {}
       void print()
       {
              cout << "(detour: " << detour << " possible_load:" << load << ") ";
              cout << endl;
       }
};

class Params {
 public:
  // generateur pseudo-aleatoire
  Rng* rng;

  // graine du generateur
  int seed;

  unsigned int jVal;

  unsigned int pHorizon;

  normal_distribution<double> dist;

  // debut de l'algo
  clock_t debut;

  // PARAMETRES DE L'ALGORITHME GENETIQUE //

  // constante d'espacement entre les fitness des individus
  double delta;

  // limite du split
  vector<double> borneSplit;

  // crit�re de proximit� des individus (RI)
  int prox;

  // crit�re de proximit� des individus (RI -- constante)
  int proxCst;

  // nombre d'individus pris en compte dans la mesure de distance
  int nbCountDistMeasure;

  // distance min
  double distMin;

  // nombre d'individus elite
  int el;

  // nombre d'individus dans la population
  int mu;

  // nombre d'offspring dans une generation
  int lambda;

  unsigned int nbScenarios;

  // probabilite de recherche locale totale pour la reparation (PVRP)
  double pRep;

  // coefficient de penalite associe a une violation de capacite
  vector<double> penalityCapa;

  // limite basse sur les indiv valides
  double minValides;

  // limite haute sur les indiv valides
  double maxValides;

  // fraction de la population conserv�e lors de la diversification
  double rho;

   std::vector<std::vector<double>> oldDemands;

  std::vector<double> meanDemands;

  std::vector<double> stdDemands;

  // PARAMETRES DE L'INSTANCE //

  // rounding convention
  bool isRoundingInteger;
  bool isRoundingTwoDigits;

  // Constant value in the objective
  double objectiveConstant_stockout;
  void computeConstant_stockout();

  // nombre de sommets clients
  unsigned int nbClients;

  // nombre de jours
  unsigned int nbDays;

  // nombre de vehicules par d�pot
  unsigned int nbVehiculesPerDep;

  bool traces;

  // nombre de depots (MDVRP)
  // correspond � l'indice du premier client dans le tableau C
  unsigned int nbDepots;

  // sequence des vehicules utilisables chaque jour avec les contraintes et
  // depots associes
  vector<vector<Vehicle>> ordreVehicules;

  // nombre de vehicules utilisables par jour
  vector<unsigned int> nombreVehicules;

  // vecteur des depots et clients 客户，depot向量
  //vector<Client> cli;
  vector<Client> cli;

  // temps de trajet , calcules lors du parsing
  vector<vector<double>> timeCost;

  // critere de correlation
  vector<vector<bool> > isCorrelated1;

  // SPECIFIC DATA FOR THE INVENTORY ROUTING PROBLEM //

  // availableSupply[t] gives the new additional supply at day t.
  // availableSupply[1] gives the initial supply + production for day 1
  vector<double> availableSupply;
  vector<double> allSupply;

  // inventory cost per day at the supplier
  double inventoryCostSupplier;

  // ROUTINES DE PARSING //

  // flux d'entree du parser
  ifstream fichier;

  // initializes the parameters of the method
  void setMethodParams();

  // effectue le prelevement des donnees du fichier
  void preleveDonnees(string nomInstance);

  Client getClientDSIRP();

  void setJval(unsigned int _jVal) {
       jVal = _jVal;
       nbDays = pHorizon - jVal + 1;
  }

  void updateToDay(unsigned int j, std::vector<double> deliveries);

  // computes the distance matrix
  void computeDistancierFromCoords();

  // calcule les autres structures du programme
  void calculeStructures();

  // modifie al�atoirement les tableaux de proximit� des clients
  void shuffleProches();

  void adjustDemands();
  // constructeur de Params qui remplit les structures en les pr�levant dans le
  // fichier
  Params(string nomInstance, int seedRNG, unsigned int nbScenario, unsigned int nbVeh, bool trace);


  // destructeur de Params
  ~Params(void);
};
#endif
