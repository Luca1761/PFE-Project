/*                       Algorithme - HGSADC                         */
/*                    Propriete de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef POPULATION_H
#define POPULATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <math.h>
#include <time.h>
#include "Node.h"
#include "Individual.h"

using namespace std ;

struct SubPopulation
{
   // individuals from sub-population
	vector<Individual*> individuals ;

	// number of individuals
	unsigned int nbIndiv ;
};

class Population
{
   private:

   // parameters of instance and genetic algorithm
   Params* params;

   // number of scenarios
   unsigned int nbScenario;

   // indicates if the last XX produced individuals were feasible (capacity point of view)
   list<bool> capaValidityList ;

   // auxiliary data structure to apply the local search
   Individual* trainer;

   // education of an individual (with split and localsearch)
   void train(Individual* indiv);

   // boolean function checking if the fitness already exists
   bool fitExist(SubPopulation* pop, Individual* indiv);

   // place an individual in the correct subpopulation and returns its place (Individuals are increasingly sorted in the subpopulation)
   unsigned int placeIndividual(SubPopulation* pop, Individual* indiv);

   public:

   // clock time when the best individual was found
   clock_t timeBest ;

   // total clock time
   clock_t totalTime ;

   // computes the extended fitness of subpopulation individuals
   void evalExtFit(SubPopulation* pop);

   // same than placeIndividual but we also remove potential individuals from the population
   // if it's full
   unsigned int addIndividual(Individual* indiv);

   // remove an individual from the population, depending on replacement policy
   void removeIndividual(SubPopulation* pop, unsigned int p);

   // update the nearest individuals of a population (because of new individual indiv)
   void updateProximity(SubPopulation* pop, Individual* indiv);

   // "restart process" where we replace a part of the population by random individuals
   void diversify();

   // copy an individual into another
   // WARNING!! copy only chromT and fitness attributes
   void copyIndividual(Individual* destination , Individual* source);
   
   // feasible individuals from the population
   SubPopulation* valid;

   // infeasible individuals from the population
   SubPopulation* invalid;

   // "getter" of 1 individual with binary tournament
   Individual* getIndividualBinT();

   // "getter" of best feasible individual
   // return NULL if no feasible
   Individual* getBestFeasIndividual();

   // "getter" of best infeasible individual
   // return NULL if no infeasible
   Individual* getBestInfeasIndividual();

   // compromise between fitness and diversity
   unsigned int selectCompromise(SubPopulation * souspop) ;

   // after penalties changes, computes the new evaluation of individuals and make a new sort of individuals
   // with their new cost
   void validatePen();

   // export the best solution into a solution file
   void ExportPop(string nomFichier, std::vector<double> deliveries, double &totalCost) ;

   // return the part of feasible individuals
   double validChargePart() ;

   // computes population diversity
   double getDiversity(SubPopulation * pop);

   // return the mean cost of feasible solutions
   double getAverageFeasible();

   // return the mean cost of infeasible solutions
   double getAverageInfeasible();

   // display the current state of population
   void displayState(unsigned int NbIter);

   // update the number of valid individuals
   void updateNbValid (Individual * indiv);

   // measure the cost of our best solution applied to true demands
   // also updates the quantities to deliver on day 1 by computing average quantities on
   // all scenarios
   void measureAndUpdateQuantities(std::vector<double> &deliveries, double &totalCost);

   // constructor
   Population(Params* params);

   // destructor
   ~Population();
};

#endif
