/*                       Algorithme - HGSADC                         */
/*                    Propri�t� de Thibaut VIDAL                     */
/*                    thibaut.vidal@cirrelt.ca                       */

#ifndef CLIENT_H
#define CLIENT_H

#include <iostream> 
#include <vector>
#include <math.h>
using namespace std ;

struct couple {
	double x;
	double y;
};

class Client // can be a customer or a supplier
{
public:

	// client index
    int custIdx;

    // client coordinates
    couple coord;

	// DATA STRUCTURES USED FOR THE IRP //
	
	// starting inventory level
	double startingInventory;

	// bound for the inventory
	double maxInventory;

	// daily demand of the customer (scenario based) dailyDemand[scenario][t] -> demand on day t for scenario
	vector<vector<double>> dailyDemand; 

	// true demand of the customer on full horizon (only for test)
	vector<double> trueDemand ;

	// daily inventory cost of the customer
	double inventoryCost ;

	// daily stockout cost of the customer
	double stockoutCost;

	// customer and depot sorted by proximity
	vector<unsigned int> proximityOrder ;

	// client's neighbors, according to proximity criteria
	vector<unsigned int> neighbors ;

	// constructor
	Client(void);

	// destructor
	~Client(void);
};

#endif
