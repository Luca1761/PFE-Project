/*                       Algorithme - HGSADC                         */
/*                    Propri�t� de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef INDIVIDU_H
#define INDIVIDU_H

/*
Classe individu : chaque individu est represente par son chromT traduisant la sequence de parcours des sommets, sans inserer la depot nulle part. Le calcul du chemin VRP reel et de son fitness se fait avec la fonction split. 

Attention, certains champs ne sont calcules qu'apres une execution de split, ne pas tenter d'y acceder avant.

Des getters auraient pu etre faits afin de restreindre l'acces mais dans un souci d'efficacite algorithmique on se contentera de la conscience du codeur. 
*/

#include <vector>
#include <list>
#include <iostream> 
#include "Noeud.h"
#include "Params.h"
#include "LocalSearch.h"
using namespace std ;

class LocalSearch ;

struct coutSol {
  // valeur du fitness comportant les penalites, si il a ete calcule
  double evaluation ;

  // valeur du fitness non penalise
  double fitness ;

  // violations de capacite
  double capacityViol ;

  double day1Cost;
};

class Individu ;

struct proxData {
  // individu en question
  Individu * indiv ;

  // sa distance
  double dist ;
};

class Individu
{

 private:

 // Acces aux parametres de l'instance et du genetique
 vector<Params*> paramsList;

 public:

  int nbScenario;

  // fitness etendu
  double fitnessEtendu ;

  // rang diversite
  float divRank ;

  // rang fitness 
  float fitRank ;

  // evaluation de la solution
  struct coutSol coutSol;
  struct coutSol_scenario {
    vector<double> evaluation;
    vector<double> fitness;
    vector<double> capacityViol;
  };
  coutSol_scenario coutSol_scenario;

  // The giant tour of each individual 
  // chromT [i][j] -> jour i, client j dans la succession des clients du jour i
  vector<vector<int>> chromT ;

  // chromL [i][j] -> The load to be delivered to each customer [j] on day [i]
  vector<vector<double>> chromL ;

  // Auxiliary data structure to run the Split
  // potentiels[i+1] -> distance pour atteindre le sommet i
  // de la sequence
  // potentiels[0] = 0 et non   
  // potentiels[1] = distance du sommet 0
  vector<vector<double>> potentiels ;

  // pour chaque jour le tableau de [nbCamions] [predecesseur]
  // potentiels[i+1] -> predecesseur de i
  vector<vector<vector<int>>> pred ;

  // says if the individual is a feasible solution
  bool estValide ;

  // distance measure
  double distance(Individu * indiv2);

  // individus class�s par proximit� dans la population, pour les politiques de remplacement
  list<proxData> plusProches;

  // ajoute un element proche dans les structures de proximit�
  void addProche(Individu * indiv) ;

  // enleve un element dans les structures de proximit� 
  void removeProche(Individu * indiv) ;

  // distance moyenne avec les n individus les plus proches
  double distPlusProche(int n) ;

  // structure de donnee associee a l'individu au sein de la recherche locale
  // seul le rejeton de Genetic.cpp possede cette structure
  // sinon elle n'est pas initialisee

  vector<LocalSearch *> localSearchList;

  // fonction Split pour tous les jours
  // essaye deja le split simple
  // si la solution ne respecte pas le nombre de camions : essaye le split a flotte limitee
  void generalSplit_scenario();

  // fonction split ne respectant pas forcement le nombre de vehicules
  // retourne 1 si succes, 0 sinon
  int splitSimple_scenario(int k, Params *paramsTemp, int scenario) ;
  int splitSimpleDay1();

  // fonction split pour problemes a flotte limitee
  void splitLF(int k) ;
  void splitLF_scenario(int k, Params *paramsTemp, int scenario);

  // fonction qui se charge d'evaluer exactement les violations
  // et de remplir tous les champs d'evaluation de solution
  void measureSol() ;
  void measureSol_scenario();

  // initialisation du vecteur potentiels
  void initPot_scenario(int k, int scenario) ;

  // mise a jour de l'objet localSearch, 
  // Attention, Split doit avoir ete calcule avant
  void updateLS_scenario() ;

  // Inverse procedure, after local search to return to a giant tour solution representation and thus fill the chromT table.
  void updateIndiv_scenario() ;

  // Computes the maximum amount of load that can be delivered to client on a day k without exceeding the 
  // customer maximum inventory
  double maxFeasibleDeliveryQuantity(int day, int client); 

  Individu(vector<Params*> pl);

  //destructeur
  ~Individu();
};
#endif
