#ifndef VEHICLE_H
#define VEHICLE_H

#include "Params.h"

class Params ;

class Vehicle
{

private:

// instance parameters
Params * params ;

public:

// identification of the related depot
unsigned int depotNumber ;     

// capacity limit
double capacity ;

// constructor
Vehicle(unsigned int _depotNumber, double _capacity);

// destructor
~Vehicle(void);

};

#endif
