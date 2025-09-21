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


int main(int argc, char *argv[]) {
  commandline command(argc, argv);
  
  std::vector<double> deliveries;
  Params* params = new Params(
    command.get_path_to_instance(), 
    command.get_seed(),
    command.get_nb_cores(),
    command.get_nb_scenarios(),
    command.get_nbVeh(),
    command.get_end_day_inventories(),
    command.get_trace(),
    command.get_deterministic(),
    command.get_true_demand()
  );
  if (params->traces)std::cout << "Read instance: " << command.get_path_to_instance() << std::endl;
  if (params->traces)std::cout << "Solution stored in: "<< command.get_path_to_solution() << std::endl;
  
  double totalCost = 0;
  double totalCostBis = 0;
  for (unsigned int j = 1; j <= params->pHorizon; j++) { // Solve on rolling horizon
    //update initial inventories/available supply/structures according to what happened previously
    params->updateToDay(j, deliveries);
    
    // initial population 
    Population *population = new Population(params);
    
    // we create the solver
    Genetic solver(params, population);

    //launch evolution
    if (params->traces) std::cout << "Solve day " << j << "..." << std::endl;
    solver.evolve(command.get_maxIter(), command.get_maxIterNonProd(), command.get_maxTime());

    population->measureAndUpdateQuantities(deliveries, totalCost);
    population->ExportPop(command.get_path_to_solution(), deliveries, totalCostBis);
    if (params->deterministic) return 0;
  }
  if (params->traces) std::cout << "FINAL COST: " << totalCost << std::endl;
  else std::cout << totalCost << std::endl;
  
  // desallocating of params memory
  delete params;
  return 0;
}