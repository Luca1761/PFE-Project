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
int depotNumber ;     

// duration limit
double maxRouteTime ;

// capacity limit
double capacity ;

Vehicle(int depotNumber,double maxRouteTime,double capacity);
~Vehicle(void);
};

#endif
