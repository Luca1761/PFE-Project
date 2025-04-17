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

	// nombre de jours
    int nbJours ;

	// customer number
    int custNum ;

    // coordonn�es des points
    couple coord ;

	// demande associ�e � un sommet
	double demand ;

	// DATA STRUCTURES USED FOR THE IRP //
	
	// starting inventory level
	double startingInventory ;

	// bounds for the inventory
	double minInventory ;

	// bounds for the inventory
	double maxInventory ;

	// daily demand of the customer.
	vector <double> dailyDemand ;

	// daily inventory cost of the customer
	double inventoryCost ;
	double stockoutCost;
	// ordre des sommets et depots, par proximit�
	vector <int> ordreProximite ;

	// les sommets les plus proches selon le critere de proximite
	vector <int> sommetsVoisins ;

	Client(void);

	~Client(void);
};

#endif
