/*                       Algorithme - HGSADC                         */
/*                    Propri�t� de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */
/*  Utilisation non autoris�e sans permission explicite des auteurs  */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <random>
#include "Genetic.h"
#include "LotSizingSolver.h"
#include "commandline.h"

using namespace std;

#define TRACES

int mainSIRP(int argc, char *argv[])
{
  commandline c(argc, argv);

  // Nombre de ticks horloge que le programme est autorise a durer
  clock_t nb_ticks_allowed; 
  nb_ticks_allowed =  CLOCKS_PER_SEC;

  int nbScenarios = c.get_nb_scenarios();
  mt19937 gen(c.get_seed());  
  normal_distribution<double> dist(0, 0);

  vector<double> randomNumbers = vector<double>(nbScenarios, 0.0);
  for (int i = 0; i < nbScenarios; i++) {
    randomNumbers[i] = dist(gen);
  }

  vector<Params*> paramsList;
  for(int i = 0; i < nbScenarios; i++) {
    double randomValue = randomNumbers[i];
    Params* param = new Params(
      c.get_path_to_instance(), 
      c.get_path_to_solution(),
      c.get_nbVeh(), 
      c.get_seed(),
      c.get_rou(), 
      c.get_stockout(),
      randomValue,   
      dist           
    );
    paramsList.push_back(param);
  }

  // Params *mesParametres = new Params(
  //     c.get_path_to_instance(), c.get_path_to_solution(), c.get_type(),
  //     c.get_nbVeh(), c.get_path_to_BKS(), c.get_seed(),c.get_rou(), c.get_stockout(), c.get_nb_scenarios());
      
  // initial population 
  Population *population = new Population(paramsList);

  // we create the solver
  Genetic solver(paramsList, population, nb_ticks_allowed, true);
  
  // on lance l'evolution   launch evolution
  
  int max_iter = 100000;
  int maxIterNonProd = 10000;
  solver.evolve(max_iter, maxIterNonProd, 1);







  
  population->ExportPop(c.get_path_to_solution(),true);
  
  population->ExportBKS(c.get_path_to_BKS());

  // on desalloue la memoire
  for (int i = 0; i < nbScenarios; i++)
    delete paramsList[i];
  return 0;
  
}

int main(int argc, char *argv[])
{
  mainSIRP(argc, argv);
}
