/*                       Algorithme - HGSADC                         */
/*                    Propriete de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

/*
Individual: each individual is represented by its chromT. ChromT contains the big tour (without depot) with clients to visit, for each day and each scenario.
ChromL contains the quantity to deliver to these clients

We compute the real VRP (with the depot and different routes from different vehicles) later with split functions, giving the real cost of the solution.

Be careful, many different fields of this class are not accessible before the general_split (particularly the cost)
*/

#include <vector>
#include <list>
#include <iostream> 
#include "Node.h"
#include "Params.h"
#include "LocalSearch.h"
#include <thread>
using namespace std ;

class LocalSearch ;

struct coutSol {
  // fitness value, including penalities (from inventory or capacity overcharge)
  double evaluation ;

  // fitness value without penalities
  double fitness ;

  // penalities
  double capacityViol ;
};

// literally the same but for each scenario (coutSol is then compute as an average of these)
struct coutSol_scenario {
  vector<double> evaluation;
  vector<double> fitness;
  vector<double> capacityViol;
};

class Individual ;

struct proxData {
  // Individual
  Individual * indiv ;

  // distance to the individual 
  double dist ;
};

class Individual
{

 private:

 // parameters of instance and genetic algorithm
 Params* params;

 public:

  // number of scenarios ( = params->nbScenario but we use it a lot so define it again)
  unsigned int nbScenario;

  // fitness etendu
  double fitnessEtendu;

  // rang diversite
  float divRank;

  // rang fitness 
  float fitRank;

  // average gloval solution cost
  struct coutSol coutSol; 

  // solution cost, scenario per scenario
  coutSol_scenario coutSol_scenario;

  // The giant tour of each individual 
  // chromT[1][j] -> day 1, j-th customer of the corresponding tour (shared by every scenario)
  // chromT[t + s * (nbDays - 1)][j], t=2...nbDays -> day t, scenario s, j-th customer of the corresponding tour
  vector<vector<unsigned int>> chromT ;

  // chromL[t + s * nbDays][j], t=1...nbDays -> day t, scenario s, the load to be delivered to customer whose index is j (not the j-th) on this day and scenario
  vector<vector<double>> chromL ;

  // Auxiliary data structure to run the Split
  // potentials[i + 1] ->  distance to join i-th customer of the tour
  // potentials[0] = 0  
  // potentials[1] = distance to 0
  vector<vector<double>> potentials ;

  // for each day, the array of size [nbVehicles] [predecessor]
  // pred[c][i+1] -> predecessor of i-th customer on route of vehicle c
  vector<vector<vector<unsigned int>>> pred ;

  // says if the individual is a feasible solution
  bool isFeasible ;

  // distance measure between two individuals
  double distance(Individual * indiv2);

  // individual sorted by proximity in their population, for replacement policies
  list<proxData> plusProches;

  // add an individual in proximity structures
  void addProche(Individual * indiv) ;

  // remove an individual in proximity structures
  void removeProche(Individual * indiv) ;

  // average distance with the n nearest Individuals
  double distPlusProche(int n);

  // one local search per scenario, only rejeton from Genetic has this structure, otherwise its not initialized
  // this structure will help to apply local search to this individual
  vector<LocalSearch*> localSearchList;

  // global split function, to split every giant tour from every scenario and day
  // we try a simple split and if the solution we get has an invalid number of vehicles, we try the limited fleed split
  void generalSplit_scenario();

  // fonction split ne respectant pas forcement le nombre de vehicules
  // retourne 1 si succes, 0 sinon
  // simple split to apply to each tour with day >= 2 for each scenario
  int splitSimple_scenario(unsigned int k, unsigned int scenario) ;

  // same but for day 1 (with average cost on each scenario)
  bool splitSimpleDay1();

  // split function with limited fleet to apply to each tour with day >= 2 for each scenario
  bool splitLF_scenario(unsigned int k, unsigned int scenario);

   // same but for day 1 (with average cost on each scenario)
  bool splitLF_scenario_day1();

  // function that evaluate violations and cost of the solution
  // also fills every evaluation field of the solution
  void measureSol_scenario();

  // measure the cost of the solution with true demand (and fills the vector delivers with the quantity we choose to deliver)
  double measureTrueCost(std::vector<double> &delivers);

  // re-initialize potentials structure
  void initPot_scenario(unsigned int k, unsigned int scenario) ;

  // update localSearch list, (fills localSearch structures with chromT and chromL)
  // BE CAREFUL: Split needs to be compute before
  void updateLS_scenario() ;

  // start the local search procedure
  void runLocalSearch_scenario();

  // start the backward dynamic programming operator
  void backwardDPOperator();

  // start the backward DP operator for a given client
  int backwardDPSingleClient(unsigned int client, bool &rechercheTerminee);

  // start all the other local search operator many times on day 1
  void runSearchDay1();

  // start an iteration of local search
  int mutationDay1();

  // Inverse procedure, after local search to return to a giant tour solution representation and thus fill the chromT table.
  void updateIndiv_scenario() ;

  // every type of day 1 mutations
  int mutation1_indiv();
  int mutation2_indiv();
  int mutation3_indiv();
  int mutation4_indiv();
  int mutation5_indiv();
  int mutation6_indiv();
  int mutation7_indiv();
  int mutation8_indiv();
  int mutation9_indiv();

  // random constructor
  Individual(Params* _params);

  // destructor
  ~Individual();
};
#endif
