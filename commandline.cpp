#include "commandline.h"

void commandline::set_instance_name(string to_parse)
{
  instance_name = to_parse;
}

void commandline::set_sortie_name(string to_parse) { sortie_name = to_parse; }

void commandline::set_BKS_name(string to_parse) { BKS_name = to_parse; }

void commandline::set_default_sorties_name(string to_parse)
{
  char caractere1 = '/';
  char caractere2 = '\\';
  int position = (int)to_parse.find_last_of(caractere1);
  int position2 = (int)to_parse.find_last_of(caractere2);
  if (position2 > position)
    position = position2;

  if (position != -1)
  {
    string directory = to_parse.substr(0, position + 1) + "diff/";
    string filename = to_parse.substr(position + 1, to_parse.length() - position - 1-4);

    sortie_name = directory + "STsol-" + filename+ "_veh-" + std::to_string(nbVeh) + "_rou-" + std::to_string(rou);
    BKS_name = directory + "STbks-" + filename+ "_veh-" + std::to_string(nbVeh) + "_rou-" + std::to_string(rou);
  }
  else
  {
    sortie_name = to_parse.substr(0, position + 1) + "STsol-" +
                  to_parse.substr(position + 1, to_parse.length() - 1) +"-seed-"+to_parse.substr(seed)+"-veh-"+to_parse.substr(nbVeh)+"-rou-"+to_parse.substr(rou);
    BKS_name = "STbks-" + to_parse;
  }
}

// constructeur
commandline::commandline(int argc, char *argv[])
{
  bool isTime = false;
  bool isOutput = false;
  bool isBKS = false;
  if (argc > 24 || argc < 2)
  {
    cout << "incorrect command line" << endl;
    throw std::string("incorrect command line");
  }
  else
  {
    // default values
    set_instance_name(string(argv[1]));

    seed = 0;
    nbVeh = -1; // unknown
    rou = -1;
    stockout = false;

    // parameters
    for (int i = 2; i < argc; i += 2)
    {
      if (string(argv[i]) == "-sol")
        set_sortie_name(string(argv[i + 1]));
      else if (string(argv[i]) == "-bks")
        set_BKS_name(string(argv[i + 1]));
      else if (string(argv[i]) == "-seed")
        seed = atoi(argv[i + 1]);
      else if (string(argv[i]) == "-veh")
        nbVeh = atoi(argv[i + 1]);
      else if (string(argv[i]) == "-stock"){
          rou = atoi(argv[i + 1]);
          stockout = true;
      } else if (string(argv[i]) == "-scenarios"){
          nb_scenarios = atoi(argv[i + 1]);
      } else {
        cout << "Commande non reconnue : " << string(argv[i]) << endl;
        throw std::string("Commande non reconnue");
      }
    }
    set_default_sorties_name(string(argv[1]));
  }
}

void commandline::set_debug_params(string instance)
{
  bool isTime = false;
  bool isOutput = false;
  bool isBKS = false;

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

string commandline::get_path_to_BKS() { return BKS_name; }

int commandline::get_stockout() { return stockout; }

int commandline::get_nb_scenarios() {return nb_scenarios; }

int commandline::get_rou() { return rou; }

// renvoie le nombre de v�hicules optimal connu
int commandline::get_nbVeh() { return nbVeh; }

// renvoie la seed
int commandline::get_seed() { return seed; }