#ifndef VEHICLE_H
#define VEHICLE_H

#include "Params.h"

class Params ;

class Vehicle
{

private:

// acces aux donnees de l'instance
Params * params ;

public:

// identification of the related depot
unsigned int depotNumber ;     

// capacity limit
double capacity ;

Vehicle(unsigned int _depotNumber, double _capacity);
~Vehicle(void);
};

#endif
