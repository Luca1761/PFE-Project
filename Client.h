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

class Client
{
public:

	// customer number
    int custIdx ;

    // coordonnees des points
    couple coord ;

	// DATA STRUCTURES USED FOR THE IRP //
	
	// starting inventory level
	double startingInventory ;

	// bounds for the inventory (generally 0)
	double minInventory;

	// bounds for the inventory
	double maxInventory ;

	// daily demand of the customer.
	vector<vector<double>> dailyDemand ;

	vector<double> testDemand ;

	// daily inventory cost of the customer
	double inventoryCost ;
	double stockoutCost;
	// ordre des sommets et depots, par proximit�
	vector<int> ordreProximite ;

	// les sommets les plus proches selon le critere de proximite
	vector<int> sommetsVoisins ;

	Client(void);

	~Client(void);
};

#endif
