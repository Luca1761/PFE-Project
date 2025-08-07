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

    // nb cores for multithreading in expes
    unsigned int nb_cores;

    // to solve with true demand (if true, nb_scenarios should be 1)
    bool true_demand;

    // if true, constraint is I_i + q_i - d_i < U_i. Else, it is I_i + q_i < U_i (beginning of the day)
    bool end_day_inventories;

    // number of scenarios
    unsigned int nb_scenarios;

    // nbVehicles, if given (1 in our setup)
    unsigned int nbVeh;

    // instance name
    string instance_name;

    // sol file name
    string sortie_name;

    // max_iteration
    unsigned int maxIter;
    
    // max non productive iteration
    unsigned int maxIterNonProd;

    // maximum time
    unsigned int maxTime;

    // show or don't show the trace
    bool traces;

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

    // seed
    int get_seed();
    
    // nb cores for multithreading
    unsigned int get_nb_cores();

    // get information to set true demand or not
    bool get_true_demand();

    // get the information about inventory time
    bool get_end_day_inventories();

    // number of scenarios
    unsigned int get_nb_scenarios(); 

    // number of vehicles
    unsigned int get_nbVeh();

    // instance path
    string get_path_to_instance();

    // solution path
    string get_path_to_solution();

    // max iterations
    unsigned int get_maxIter();

    // max non productive iterations
    unsigned int get_maxIterNonProd();

    // max time
    unsigned int get_maxTime();

    // get trace
    bool get_trace();
};
#endif
