#include "commandline.h"

void commandline::set_instance_name(string to_parse) { instance_name = to_parse; }

void commandline::set_sortie_name(string to_parse) { sortie_name = to_parse; }

void commandline::set_default_sorties_name(string to_parse) {
  char caractere1 = '/';
  char caractere2 = '\\';
  int pos = (int) to_parse.find_last_of(caractere1);
  int pos2 = (int) to_parse.find_last_of(caractere2);
  if (pos2 > pos)
    pos = pos2;

  if (pos != -1) {
    unsigned int position = (unsigned int) pos;
    string directory = (true_demand) ? "irp_results/" : "dsirp_results/";
    string filename = to_parse.substr(position + 1, to_parse.length() - position - 1-4);
    string type = (true_demand) ? "IRPsol-" : "DSIRPsol-";
    if (end_day_inventories) type += "end-";

    sortie_name = directory + type + filename + "_seed-" + std::to_string(seed) + "_veh-" + std::to_string(nbVeh) + "_scenarios-" + std::to_string(nb_scenarios);
  } else {
    sortie_name = to_parse.substr(0, 0) + "DSIRPsol-" +
                  to_parse.substr(0, to_parse.length() - 1) + "-seed-" + to_parse.substr((unsigned int) seed) + "-veh-" + to_parse.substr(nbVeh);
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
    true_demand = false;
    end_day_inventories = false;
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
        nb_cores = (unsigned int) atoi(argv[i + 1]);
      else if (string(argv[i]) == "-endDay")
        end_day_inventories = (atoi(argv[i + 1]) == 1);
      else if (string(argv[i]) == "-trueDemand")
        true_demand = (atoi(argv[i + 1]) == 1);
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
    if (true_demand) nb_scenarios = 1;
    if (end_day_inventories) {
      nb_scenarios = 1;
      true_demand = true;
    }
    if (!hasSolName)
      set_default_sorties_name(string(argv[1]));
  }
}

commandline::~commandline() {}

int commandline::get_seed() { return seed; }

unsigned int commandline::get_nb_cores() { return nb_cores; }

bool commandline::get_true_demand() {return true_demand;}

bool commandline::get_end_day_inventories() {return end_day_inventories;}

unsigned int commandline::get_nb_scenarios() {return nb_scenarios; }

unsigned int commandline::get_nbVeh() { return nbVeh; }

string commandline::get_path_to_instance() { return instance_name; }

string commandline::get_path_to_solution() { return sortie_name; }

unsigned int commandline::get_maxIter(){ return  maxIter;}

unsigned int commandline::get_maxIterNonProd() {return maxIterNonProd;}

unsigned int commandline::get_maxTime() {return maxTime;}

bool commandline::get_trace() {return traces;}
