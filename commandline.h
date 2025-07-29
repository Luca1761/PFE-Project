#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

// This class parses user's command line
class commandline
{
private:

    // seed
    int seed;

    // number of scenarios
    unsigned int nb_scenarios;

    // nbVehicles, if given (1 in our setup)
    unsigned int nbVeh;

    // instance name
    string instance_name;

    // nom de la sortie
    string sortie_name;

    // fill instance_name variable
    void set_instance_name(string to_parse);

    // fill sortie_name variable
    void set_sortie_name(string to_parse);

    // gives a default value to solution file name (using instance name)
    void set_default_sorties_name(string to_parse);

public:
    // constructor
    commandline(int argc, char *argv[]);

    // destructor
    ~commandline();

    void set_debug_params(string instance);

    // instance path
    string get_path_to_instance();

    // solution path
    string get_path_to_solution();

    // number of scenarios
    unsigned int get_nb_scenarios(); 

    // number of vehicles
    unsigned int get_nbVeh();

    // seed
    int get_seed();

};
#endif
