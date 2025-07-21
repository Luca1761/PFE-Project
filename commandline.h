#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

// Cette classe parse la ligne de commande entre par le user
class commandline
{
private:

    // seed
    int seed;
    int rou;

    int nb_scenarios;
    int var;

    // nbVehicles, if given
    int nbVeh;

    // nom de l'instance
    string instance_name;

    // nom de la sortie
    string sortie_name;

    // nom de la BKS
    string BKS_name;

    // remplit l'attribut instance_name
    void set_instance_name(string to_parse);

    // remplit l'attribut sortie_name
    void set_sortie_name(string to_parse);

    // remplit l'attribut BKS_name
    void set_BKS_name(string to_parse);

    // donne une valeur par defaut a la solution
    // en fonction du nom de l'instance
    void set_default_sorties_name(string to_parse);

public:
    // constructeur
    commandline(int argc, char *argv[]);

    // destructeur
    ~commandline();

    void set_debug_params(string instance);

    // renvoie le chemin de l'instance
    string get_path_to_instance();

    // renvoie le chemin vers la solution
    string get_path_to_solution();

    // renvoie le chemin vers la meilleure solution connue
    string get_path_to_BKS();

    // renvoie le temps cpu allou
    int get_cpu_time();
    int get_rou();
    int get_nb_scenarios(); 
    int get_var();

    // renvoie le nombre de v�hicules optimal connu
    int get_nbVeh();

    // renvoie la seed
    int get_seed();

};
#endif
