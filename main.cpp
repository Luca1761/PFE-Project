/*                       Algorithme - HGSADC                         */
/*                    Propriete de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */
/*  Utilisation non autorisee sans permission explicite des auteurs  */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <random>
#include "Genetic.h"
#include "LotSizingSolver.h"
#include "commandline.h"

using namespace std;

#define TRACES

int main(int argc, char *argv[]) {
  commandline command(argc, argv);
  
  // Nombre de ticks horloge que le programme est autorise a durer
  clock_t nb_ticks_allowed; 
  nb_ticks_allowed =  CLOCKS_PER_SEC;
  int max_iter = 505;
  int maxIterNonProd = 505;
  
  std::vector<double> deliveries;
  Params* params = new Params(
    command.get_path_to_instance(), 
    command.get_seed(),
    command.get_nb_scenarios(),
    command.get_nbVeh()
  );
  std::cout << "Read instance: " << command.get_path_to_instance() << std::endl;
  std::cout << "Solution stored in: "<< command.get_path_to_solution() << std::endl;
  
  double totalCost = 0;
  for (unsigned int j = 1; j <= params->pHorizon; j++) { // Solve on rolling horizon
    //update initial inventories/available supply/structures according to what happened previously
    params->updateToDay(j, deliveries);
    
    // initial population 
    Population *population = new Population(params);
    
    // we create the solver
    Genetic solver(params, nb_ticks_allowed, population, true);
    
    //launch evolution
    std::cout << "Solve day " << j << "..." << std::endl;
    solver.evolve(max_iter, maxIterNonProd);

    // population->ExportPop(c.get_path_to_solution(),true);

    population->measureAndUpdateQuantities(deliveries, totalCost, j);
  }
  std::cout << "FINAL COST: " << totalCost << std::endl;
  
  // desallocating of params memory
  delete params;
  return 0;
}