/*                       Algorithme - HGSADC                         */
/*                    Propriete de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef GENETIC_H
#define GENETIC_H

// class for genetic algorithm

#include "Population.h"
#include "Params.h"
#include "Individual.h"
#include "time.h"
#include <stdlib.h>
#include <stdio.h> 
#include <vector>
#include <list>
#include <math.h>
#include <random>
using namespace std ;

class Genetic
{
private:
    // global parameters
	Params* params;

	// population of solutions to work on
	Population* population ;	

	// Number of non-improving iterations before stopping
	unsigned int nbIterNonProd ;

	// Number of iterations before stopping
	unsigned int nbIter ;
	
	// Number of scenarios
	unsigned int nbScenario;

	// First individual to be used as input for the crossover
	Individual* child; 

	// Second individual to be used as input for the crossover
	Individual* child2;

	// launch the mutation
	void mutate();

	// repair the solution if it's not feasible
	void repair();

	// crossover
	void crossPOX(); // for one scenario (deterministic for example)
	void crossPOX_scenario(); // for multiple scenarios
	
	// manage penalties
	void managePenalties();

public:
	// launch the genetic algorithm using parameters params on population
	// the loop stops when we reach maxIterProd iterations or maxIterNonProd of non productive iterations
	// (or if we reach maxTime)
    void evolve(unsigned int maxIter, unsigned int maxIterNonProd, unsigned int maxTime) ;

	// constructor
	Genetic(Params* _params, Population * _population);

	// destructor
	~Genetic(void);
};

#endif
