/*                       Algorithme - HGSADC                         */
/*                    Propri�t� de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef GENETIC_H
#define GENETIC_H

// definition de la structure d'un algorithme genetique
// un algorithme genetique contient ses parametres (Params)
// ainsi que l'adresse de la population sur laquelle il travaille

#include "Population.h"
#include "Params.h"
#include "Individu.h"
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
    //global parameters
	Params* params;

	//population of solutions
	Population* population ;	

	// Number of non-improving iterations before stopping
	unsigned int nbIterNonProd ;

	// Number of iterations before stopping
	unsigned int nbIter ;
	
	// Number of scenarios
	unsigned int nbScenario;

	// First individual to be used as input for the crossover
	Individu* rejeton ; 

	// Second individual to be used as input for the crossover
	Individu* rejeton2 ;

	// effectue la mutation
	void muter_scenario();

	// eventuellement effectue une reparation de la solution
	void reparer_scenario();

	// procedure de crossover
	// retourne -1 si l'individu ainsi cree n'est pas valide
	void crossPOX2();
	int crossPOX_scenario();
	
	// gestion des penalites
	void gererPenalites_scenario();

public:

    // lancement de l'algorithme sur les parametres "params" et la population "population"
	// la boucle s'arrete lorsque l'on a effectue maxIterProd iterations productives
	// ou maxIterNonProd iterations non productives 
    void evolve(unsigned int maxIter, unsigned int maxIterNonProd, unsigned int maxTime) ;

	// constructeur du solveur genetique
	Genetic(Params* _params, Population * _population);

	// destructeur
	~Genetic(void);
};

#endif
