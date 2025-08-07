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

class Vehicle;
class Node;

// needed structure for a few places in the code (easily accessible from here)
struct Insertion {
       // the detour cost of this insertion
       double detour;

       //the remain load of the route where we want to insert
       double load;

       // the place of this insertion in this route
       Node *place;

       Insertion() {
              detour = 1.e30;
              load = -1.e30;
              place = NULL;
       }
       Insertion(double _detour, double _load, Node *_place)
           : detour(_detour), load(_load), place(_place) {}
       void print() {
              cout << "(detour: " << detour << " possible_load:" << load << ") ";
              cout << endl;
       }
};

class Params {
 public:
  // pseudo-random generator
  Rng* rng;

  // generator's seed
  int seed;

  // number of cores available for multithreading
  unsigned int nbCores;

  // number of scenarios
  unsigned int nbScenarios;

  // current day in the rolling horizon
  unsigned int jVal;

  // number of days
  unsigned int nbDays;

  void setJval(unsigned int _jVal) {
       jVal = _jVal;
       nbDays = pHorizon - jVal + 1;
  }

  // in the rolling horizon, update every variable according to what happened the previous day
  void updateToDay(unsigned int j, std::vector<double> deliveries);
  
  // true if inventories are taken at the end of the day
  bool endDayInventories;
  
  // verbosity of the algorithm
  bool traces;

  /////////////////////////////////////////// GENETIC ALGORITHM PARAMETERS ///////////////////////////////////////////

  // individual fitness difference
  double delta;

  // split limit (how far to search beyond the capacity limit, one by scenario)
  vector<double> splitBounds;

  // first proximity criterion for inidividuals (granularity parameter - expressed as a percentage of the problem size)
  unsigned int prox;

  // second proximity criterion for inidividual (granularity parameter - expressed as a fixed maximum)
  unsigned int proxCst;

  // number of individuals to take into account in distance measurement
  int nbCountDistMeasure;

  // distance in terms of objective function under which the solutions are considered to be the same
  double distMin;

  // number of elite individuals
  int el;

  // number of individuals in a population
  unsigned int mu;

  // offspring number of individuals in a population
  unsigned int lambda;

  // probability to repair
  double pRep;

  // penality coefficient for a capacity violation (one per scenario)
  vector<double> penalityCapa;

  // maximal number of feasible solutions
  double minFeasibles;

  // maximal number of feasible solutions
  double maxFeasibles;

  // population fraction we keep during diversification
  double rho;

  /////////////////////////////////////////// INSTANCE PARAMETERS ///////////////////////////////////////////
  
  // rounding convention
  bool isRoundingInteger;
  bool isRoundingTwoDigits;
  
  // constant value in the objective
  double objectiveConstant;
  void computeConstant();
  
  // number of clients
  unsigned int nbClients;

  // rolling horizon
  unsigned int pHorizon;
  
  // number of vehicles per depot
  unsigned int nbVehiculesPerDep;
  
  // historical demand data (50 + jVal days in our instances) for each client
  std::vector<std::vector<double>> oldDemands;

  // mean demand from historical data for each client
  std::vector<double> meanDemands;

  // std demand from historical data for each client
  std::vector<double> stdDemands;

  // compute the demands for every client (according to mean and std or to true demand)
  void adjustDemands();

  // if true, the true demand is used (only one scenario)
  bool deterministic;

  // if true, the true demand of day 1 is used
  bool trueDemandDay1;

  // number of depots (1 in this work) | index 0 
  unsigned int nbDepots;

  // sequence des vehicules utilisables chaque jour avec les contraintes et
  // depots associes
  vector<vector<Vehicle>> vehicleOrder;

  // nombre de vehicules utilisables par jour
  vector<unsigned int> vehicleNumber;

  // depot and customers vector
  vector<Client> cli;

  // travel times, computed during parsing
  vector<vector<double>> timeCost;

  // correlation criterion
  vector<vector<bool>> isCorrelated;

  // availableSupply[t] gives the new additional supply at day t.
  // availableSupply[1] gives the initial supply + production for day 1 (starting at jVal)
  vector<double> availableSupply; 
  vector<double> allSupply; // additional supply on the whole rolling horizon

  // inventory cost per day at the supplier
  double inventoryCostSupplier;

  /////////////////////////////////////////// PARSING ///////////////////////////////////////////

  // file to parse
  ifstream file;

  // initializes the parameters of the method
  void setMethodParams();

  // collect all the data from the file
  void collectData();

  // get the next client from the file
  Client getNextClient();

  // computes the distance matrix
  void computeDistancesFromCoords();

  // compute the other structures for the algorithm
  void computeStructures();

  // randomly shuffle proximity arrays
  void shuffleProches();
  
  // constructor
  Params(string nomInstance, int seedRNG, unsigned int nbCores, unsigned int nbScenario, unsigned int nbVeh, bool trace, bool trueDemand);

  // destructor
  ~Params(void);
};
#endif
