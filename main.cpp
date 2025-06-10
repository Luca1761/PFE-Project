/*                       Algorithme - HGSADC                         */
/*                    Propri�t� de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */
/*  Utilisation non autoris�e sans permission explicite des auteurs  */

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "Genetic.h"
#include "LotSizingSolver.h"
#include "commandline.h"

using namespace std;

#define TRACES

int mainSIRP(int argc, char *argv[])
{
  commandline c(argc, argv);

  if (c.is_valid())
  {
    // Nombre de ticks horloge que le programme est autorise a durer
    clock_t nb_ticks_allowed; 
    nb_ticks_allowed =  CLOCKS_PER_SEC;

    // initialisation des Parametres � partir du fichier d'instance et du chemin
    // de la solution

    Params *mesParametres = new Params(
        c.get_path_to_instance(), c.get_path_to_solution(), c.get_type(),
        c.get_nbVeh(), c.get_path_to_BKS(), c.get_seed(),c.get_rou(), c.get_stockout(), c.get_nb_scenarios());


        
        
    // initial population 
    Population *population = new Population(mesParametres);

    // on cree le solver we create the solver,we create the solver we create the solver
    Genetic solver(mesParametres, population, nb_ticks_allowed, true, true);
    
    // on lance l'evolution   launch evolution
    
    int max_iter = 100000;
    int maxIterNonProd = 10000;
    solver.evolve(max_iter, maxIterNonProd, 1);

    population->ExportPop(c.get_path_to_solution(),true);
    
    population->ExportBKS(c.get_path_to_BKS());
  
    // on desalloue la memoire
    delete mesParametres;
    return 0;
  }
  else
    throw string(
        "ligne de commande non parsable, Usage : genvrp instance [-t cpu-time] "
        "[-sol solution]");
}

int main(int argc, char *argv[])
{
  mainSIRP(argc, argv);
}
