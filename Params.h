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
// #include "Noeud.h"
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
       Insertion(double detour, double load, Noeud *place)
           : detour(detour), load(load), place(place) {}
       void print()
       {
              cout << "(detour: " << detour << " possible_load:" << load << ") ";
              cout << endl;
       }
};

class Params {
 public:
  // g�n�rateur pseudo-aleatoire
  Rng* rng;

  // graine du g�n�rateur
  int seed;

  normal_distribution<double> dist;

  // adresse de l'instance
  string pathToInstance;

  // adresse de la solution
  string pathToSolution;

  // adresse de la BKS
  string pathToBKS;

  // debut de l'algo
  clock_t debut;

  // PARAMETRES DE L'ALGORITHME GENETIQUE //

  // constante d'espacement entre les fitness des individus
  double delta;

  // limite du split
  double borneSplit;

  // crit�re de proximit� des individus (RI)
  int prox;

  // crit�re de proximit� des individus (RI -- constante)
  int proxCst;

  // crit�re de proximit� des individus (PI -- constante)
  int prox2Cst;

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

  int nbScenarios;

  // probabilite de recherche locale totale pour la reparation (PVRP)
  double pRep;

  // coefficient de penalite associe a une violation de capacite
  double penalityCapa;

  // limite basse sur les indiv valides
  double minValides;

  // limite haute sur les indiv valides
  double maxValides;

  // fraction de la population conserv�e lors de la diversification
  double rho;

  // PARAMETRES DE L'INSTANCE //

  // rounding convention
  bool isRoundingInteger;
  bool isRoundingTwoDigits;

  // Constant value in the objective
  double objectiveConstant_stockout;
  void computeConstant_stockout();

  // pr�sence d'un probl�me IRP ;
  bool isstockout;

  // nombre de sommets clients
  int nbClients;

  // nombre de jours
  int nbDays;

  // ancien nombre de jours
  int ancienNbDays;

  // nombre de vehicules par d�pot
  int nbVehiculesPerDep;

  // nombre de depots (MDVRP)
  // correspond � l'indice du premier client dans le tableau C
  int nbDepots;

  // sequence des vehicules utilisables chaque jour avec les contraintes et
  // depots associes
  vector<vector<Vehicle>> ordreVehicules;

  // nombre de vehicules utilisables par jour
  vector<int> nombreVehicules;

  // vecteur des depots et clients 客户，depot向量
  //vector<Client> cli;
  vector<Client> cli;

  // temps de trajet , calcules lors du parsing
  vector<vector<double>> timeCost;

  // critere de correlation
  vector<vector<bool> > isCorrelated1;

  // critere de correlation
  vector<vector<bool> > isCorrelated2;

  // SPECIFIC DATA FOR THE INVENTORY ROUTING PROBLEM //

  // availableSupply[t] gives the new additional supply at day t.
  // availableSupply[1] gives the initial supply + production for day 1
  vector<double> availableSupply;

  // inventory cost per day at the supplier
  double inventoryCostSupplier;

  // TRANSFORMATIONS D'INSTANCES //

  // table de correspondance : le client i dans le nouveau pb correspond �
  // correspondanceTable[i] dans l'ancien
  vector<int> correspondanceTable;

  // table de correspondance : le client i dans le nouveau pb correspond aux
  // elements de correspondanceTable[i] (dans l'ordre) dans l'ancien
  // utile lorsque des d�compositions de probl�me avec shrinking sont
  // envisag�es.
  vector<vector<int> > correspondanceTableExtended;

  // table de correspondance : le client i dans l'ancien pb correspond �
  // correspondanceTable2[i] dans le nouveau
  vector<int> correspondanceTable2;

  // ROUTINES DE PARSING //

  // flux d'entree du parser
  ifstream fichier;

  // initializes the parameters of the method
  void setMethodParams();

  // effectue le prelevement des donnees du fichier
  void preleveDonnees(string nomInstance,int rou, bool stockout);

  // sous routine du prelevement de donnees
  Client getClient(int i,int rou);

  // computes the distance matrix
  void computeDistancierFromCoords();

  // calcule les autres structures du programme
  void calculeStructures();

  // modifie al�atoirement les tableaux de proximit� des clients
  void shuffleProches();

  void adjustDemands(double randomValue);
  // constructeur de Params qui remplit les structures en les pr�levant dans le
  // fichier
  Params(string nomInstance, string nomSolution, int nbVeh,
         string nomBKS, int seedRNG);
  Params(string nomInstance, string nomSolution, int nbVeh, int seedRNG, int rou,bool stockout, 
          double randomValue);

  // Transformation de probl�me, le nouveau fichier params cr�� correspond � un
  // sous-probl�me:
  // et est pr�t � �tre r�solu ind�pendamment
  // si decom = 2 -> depots fix�s, extraction du PVRP associ� au d�pot (MDPVRP
  // -> PVRP et MDVRP -> VRP).
  // si decom = 1 -> patterns fix�s, extraction du VRP associ� au jour "jour",
  // (PVRP->VRP), (SDVRP->VRP)
  // si decom = 0 -> on extrait un probl�me de VRP, qui contient l'ensemble de
  // clients debutSeq ... finSeq (debut et finseq sont des valeurs et non des
  // indices). (VRP->VRP)
  Params(Params* params, int decom, int* serieVisites, Vehicle** serieVehicles,
         int* affectDepots, int* affectPatterns, int depot, int jour,
         int nbVisites, int nbVeh);
  

  void decomposeRoutes(Params* params, int* serieVisites,
                       Vehicle** serieVehicles, int nbVisites, int nbVeh);


  // destructeur de Params
  ~Params(void);
};
#endif
