#include "commandline.h"

void commandline::set_instance_name(string to_parse)
{
  instance_name = to_parse;
}

void commandline::set_sortie_name(string to_parse) { sortie_name = to_parse; }

void commandline::set_default_sorties_name(string to_parse)
{
  char caractere1 = '/';
  char caractere2 = '\\';
  int pos = (int) to_parse.find_last_of(caractere1);
  int pos2 = (int) to_parse.find_last_of(caractere2);
  if (pos2 > pos)
    pos = pos2;

  if (pos != -1) {
    unsigned int position = (unsigned int) pos;
    string directory = "dsirp_results/";
    string filename = to_parse.substr(position + 1, to_parse.length() - position - 1-4);

    sortie_name = directory + "STsol-" + filename+ "_veh-" + std::to_string(nbVeh);
  } else {
    sortie_name = to_parse.substr(0, 0) + "STsol-" +
                  to_parse.substr(0, to_parse.length() - 1) +"-seed-"+to_parse.substr((unsigned int) seed)+"-veh-"+to_parse.substr(nbVeh);
  }
}

// constructeur
commandline::commandline(int argc, char *argv[])
{
  if (argc > 24 || argc < 2)
  {
    std::cout << "incorrect command line" << std::endl;
    throw std::string("incorrect command line");
  }
  else
  {
    // default values
    set_instance_name(string(argv[1]));

    seed = 0;
    nb_cores = 1;
    nbVeh = 1; // unknown
    maxTime = 10;
    maxIter = 100;
    maxIterNonProd = 100;
    traces = true;
    bool hasSolName = false;

    // parameters
    for (int i = 2; i < argc; i += 2)
    {
      if (string(argv[i]) == "-sol") {
        set_sortie_name(string(argv[i + 1]));
        hasSolName = true;
      } else if (string(argv[i]) == "-seed")
        seed = atoi(argv[i + 1]);
      else if (string(argv[i]) == "-cores")
        nb_cores = (unsigned int) atoi(argv[i+1]);
      else if (string(argv[i]) == "-traces")
        traces = (atoi(argv[i + 1]) == 1);
      else if (string(argv[i]) == "-time")
        maxTime = (unsigned int) atoi(argv[i + 1]);
      else if (string(argv[i]) == "-iter")
        maxIter = (unsigned int) atoi(argv[i + 1]);
      else if (string(argv[i]) == "-iterNonProd")
        maxIterNonProd = (unsigned int) atoi(argv[i + 1]);
      else if (string(argv[i]) == "-veh") {
        int nbV = atoi(argv[i + 1]);
        if (nbV < 1) {
          std::cout << "Needs at least one vehicle" << std::endl;
          throw std::string("Needs at least one vehicle");
        }
        nbVeh = (unsigned int) nbV;
      } else if (string(argv[i]) == "-scenarios"){
          int nbS = atoi(argv[i + 1]);
          if (nbS < 1) {
            std::cout << "Needs at least one scenario" << std::endl;
            throw std::string("Needs at least one scenario");
          }
          nb_scenarios = (unsigned int) nbS;
      } else {
        std::cout << "Commande non reconnue : " << string(argv[i]) << std::endl;
        throw std::string("Commande non reconnue");
      }
    }
    if (!hasSolName)
      set_default_sorties_name(string(argv[1]));
  }
}

void commandline::set_debug_params(string instance)
{

  set_instance_name(instance);
  set_default_sorties_name(instance);

  seed = 1000;
  nbVeh = 2;  // unknown
}

// destructeur
commandline::~commandline() {}

// renvoie le chemin de l'instance

string commandline::get_path_to_instance() { return instance_name; }

string commandline::get_path_to_solution() { return sortie_name; }

unsigned int commandline::get_nb_scenarios() {return nb_scenarios; }

// renvoie le nombre de v�hicules optimal connu
unsigned int commandline::get_nbVeh() { return nbVeh; }

// renvoie la seed
int commandline::get_seed() { return seed; }

// renvoie la seed
unsigned int commandline::get_nb_cores() { return nb_cores; }

// max iterations
unsigned int commandline::get_maxIter(){ return  maxIter;}

// max non productive iterations
unsigned int commandline::get_maxIterNonProd() {return maxIterNonProd;}

// max time
unsigned int commandline::get_maxTime() {return maxTime;}

// traces
bool commandline::get_trace() {return traces;}