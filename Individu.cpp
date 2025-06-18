#include "Individu.h"

// constructeur d'un Individu comme simple conteneur

Individu::Individu(vector<Params*> pl) : paramsList(pl)
{
	nbScenario = pl.size();
	vector<int> tempVect;
	vector<vector<int>> tempVect2;
	coutSol.evaluation = 0;
	coutSol.fitness = 0;
	coutSol.capacityViol = 0;
	coutSol.lengthViol = 0;
	int j, temp;

	// IRP, this construction has to be fully changed.
	// There is no more concept of pattern
	// So, for now let's simply assume that we serve all customers in one delivery on day 1 (we ignore the product availability at the supplier for now)
	// Later on we will arrange this construction
	chromT = vector<vector<int>>(paramsList[0]->nbDays + 1);
	chromL = vector<vector<double>>(paramsList[0]->nbDays + 1, vector<double>(paramsList[0]->nbClients + paramsList[0]->nbDepots, 0.));
	chromR = vector<vector<int>>(paramsList[0]->nbDays + 1, vector<int>(paramsList[0]->nombreVehicules[1], -1));

	// OPTION 2 -- JUST IN TIME POLICY //

	double startInventory;
	double dailyDemand;
	for (int i = paramsList[0]->nbDepots; i < paramsList[0]->nbClients + paramsList[0]->nbDepots; i++)
	{
		startInventory = paramsList[0]->cli[i].startingInventory;
		
		for (int k = 1; k <= paramsList[0]->nbDays; k++)
		{
			double currentDayClientDemand = paramsList[0]->cli[i].dailyDemand[k];
			double nextDayClientDemand = paramsList[0]->cli[i].dailyDemand[(k + 1) % paramsList[0]->nbDays];
			// startInventory += chromL[k][i];

			if (startInventory >= currentDayClientDemand)
			{
				// enough initial inventory, no need to service
				chromL[k][i] = 0;

				bool isInventoryEnoughForNextDay = startInventory - currentDayClientDemand >= nextDayClientDemand;
				if (k < paramsList[0]->nbDays && !isInventoryEnoughForNextDay) {
					bool shouldDeliveryForNextDay = paramsList[0]->rng->genrand64_real1() <= 0.3;
					if (shouldDeliveryForNextDay) {
						chromL[k][i] = nextDayClientDemand + startInventory - currentDayClientDemand;
						chromT[k].push_back(i);
					}
				}
				startInventory = startInventory + chromL[k][i] - currentDayClientDemand;

			}
			else
			{
				// not enough initial inventory, just in time policy for the initial solution
				dailyDemand = paramsList[0]->cli[i].dailyDemand[k] - startInventory;
				startInventory = 0;
				chromL[k][i] = dailyDemand;
				chromT[k].push_back(i);
			}
		}
	}

	// And shuffle the whole solution
	for (int k = 1; k <= paramsList[0]->nbDays; k++)
	{
		for (int i = 0; i <= (int)chromT[k].size() - 1; i++)
		{
			j = i + paramsList[0]->rng->genrand64_int64() % ((int)chromT[k].size() - i); 
			// swap i,j
			temp = chromT[k][i];
			chromT[k][i] = chromT[k][j];
			chromT[k][j] = temp;
		}
	}

	// initialisation of the other structures

	for (int k = 0; k <= paramsList[0]->nbDays; k++)
	{
		pred.push_back(tempVect2);
		for (int i = 0; i < paramsList[0]->nombreVehicules[k] + 1; i++)
		{
			pred[k].push_back(tempVect);
			pred[k][i].push_back(0);
			for (int j = 0; j < (int)paramsList[0]->nbClients + paramsList[0]->nbDepots + 1; j++)
				pred[k][i].push_back(0);
		}
	}
	vector<double> potTemp;
	for (int i = 0; i <= paramsList[0]->nombreVehicules[1]; i++)
	{
		potentiels.push_back(potTemp);
		for (int j = 0; j < (int)paramsList[0]->nbClients + paramsList[0]->nbDepots + 1; j++)
			potentiels[i].push_back(1.e30);
	}
	potentiels[0][0] = 0;

	for (int i = 0; i < paramsList[0]->nbDepots + paramsList[0]->nbClients; i++)
	{
		suivants.push_back(tempVect);
		precedents.push_back(tempVect);
	}
	localSearchList = vector<LocalSearch *>(nbScenario, nullptr);
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario] = new LocalSearch(paramsList[scenario], this);
	}
}

// destructeur
Individu::~Individu()
{
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario] != NULL) {
			delete localSearchList[scenario];
		}
	}
}

// The Split is done as follows, we test if it's possible to split without setting a fixed limit on the number of vehicles
// If the resulting solution satisfies the number of vehicle, it's perfectly fine, we return it
// Otherwise we call the Split version with limited fleet (which is a bit slower).
void Individu::generalSplit()
{
	coutSol.evaluation = 0;

	// lancement de la procedure split pour chaque jour
	// on essaye d�ja le split simple, si c'est sans succes , le split LF
	for (int k = 1; k <= paramsList[0]->nbDays; k++)
		if (chromT[k].size() != 0)
		{
			int enoughVehicle = splitSimple(k);
			if (enoughVehicle == 0)
			{
				splitLF(k);
			}
		}

	// After Split
	// we call a function that fills all other data structures and computes the cost
	measureSol();

	if (paramsList[0]->borneSplit > 1000)
		throw string("Erreur Split");

	if (coutSol.evaluation > 1.e20)
	{
		// if no feasible solution has been found,
		// then we relax the limit on the max capacity violation adn try again (which is initially set to 4*Q)
		// it's a very exceptional case (mostly for the PVRP, should not happen for the CVRP).
		cout << " Impossible de Split, augmentation de l'acceptation du split " << endl;
		paramsList[0]->borneSplit *= 1.1;
		generalSplit();
	}
}

// fonction split which does not consider a limit on the number of vehicles
// just uses the line "1" of the "potentiels" table.
int Individu::splitSimple(int k)
{
	// on va utiliser la ligne 1 des potentiels et structures pred
	double load, distance, cost;
	int s0, s1, sb, j;
	potentiels[1][0] = 0;
	s0 = paramsList[0]->ordreVehicules[k][0].depotNumber;

	for (int i = 0; i < (int)chromT[k].size(); i++)
	{
		load = 0;
		distance = 0;
		j = i;
		while (j < (int)chromT[k].size() && load <= paramsList[0]->ordreVehicules[k][0].capacity * paramsList[0]->borneSplit)
		{
			s1 = chromT[k][j];
			load += chromL[k][s1];
			if (i == j)
			{
				distance = paramsList[0]->timeCost[s0][s1];
			}
			else
			{
				sb = chromT[k][j - 1];
				distance += paramsList[0]->timeCost[sb][s1];
			}

			// computing the penalized cost
			cost = distance + paramsList[0]->timeCost[s1][s0];
			if (load > paramsList[0]->ordreVehicules[k][0].capacity)
				cost += (load - paramsList[0]->ordreVehicules[k][0].capacity) * paramsList[0]->penalityCapa;

			if (potentiels[1][i] + cost < potentiels[1][j + 1]) // basic Bellman algorithm
			{
				potentiels[1][j + 1] = potentiels[1][i] + cost;
				pred[k][1][j + 1] = i;
			}
			j++;
		}
	}

	// testing if le number of vehicles is correct
	// in addition, the table pred is updated to keep track of everything
	j = (int)chromT[k].size();
	for (int jj = 0; jj < paramsList[0]->nombreVehicules[k]; jj++)
	{
		pred[k][paramsList[0]->nombreVehicules[k] - jj][j] = pred[k][1][j];
		j = pred[k][paramsList[0]->nombreVehicules[k] - jj][j];
	}

	coutSol.evaluation = -1.e30; // just for security, making sure this value is not used (as it does not contain the constants)
	initPot(k);
	if (j == 0)
		return 1;
	else
		return 0;
}

// fonction split pour probl�mes � flotte limit�e
void Individu::splitLF(int k)
{
	double load, distance, cost;
	int sb, s0, s1, i, j;

	// pour chaque camion
	for (int cam = 0; cam < paramsList[0]->nombreVehicules[k]; cam++)
	{
		i = 0;
		s0 = paramsList[0]->ordreVehicules[k][cam].depotNumber;
		while (i < (int)chromT[k].size() && potentiels[cam][i] < 1.e29)
		{
			if (potentiels[cam][i] < potentiels[cam + 1][i])
			{
				potentiels[cam + 1][i] = potentiels[cam][i];
				pred[k][cam + 1][i] = i;
			}
			load = 0;
			distance = 0;
			j = i;

			while (j < (int)chromT[k].size() && load <= paramsList[0]->ordreVehicules[k][cam].capacity * paramsList[0]->borneSplit)
			{
				s1 = chromT[k][j];
				load += chromL[k][s1];
				if (i == j)
				{
					distance = paramsList[0]->timeCost[s0][s1];
				}
				else
				{
					sb = chromT[k][j - 1];
					distance += paramsList[0]->timeCost[sb][s1];
				}

				// computing the penalized cost
				cost = distance + paramsList[0]->timeCost[s1][s0];
				if (load > paramsList[0]->ordreVehicules[k][cam].capacity)
					cost += (load - paramsList[0]->ordreVehicules[k][cam].capacity) * paramsList[0]->penalityCapa;

				if (potentiels[cam][i] + cost < potentiels[cam + 1][j + 1]) // Basic Bellman iteration
				{
					potentiels[cam + 1][j + 1] = potentiels[cam][i] + cost;
					pred[k][cam + 1][j + 1] = i;
				}
				j++;
			}
			i++;
		}
	}

	// on ajoute le fitness du jour donn�
	coutSol.evaluation = -1.e30; // just for security, making sure this value is not used (as it does not contain the constants)

	// on nettoye ce que l'on a d�plac�
	initPot(k);
}

void Individu::measureSol()
{
	int depot;
	int i, j;
	double distance, load;
	coutSol.fitness = 0;
	coutSol.capacityViol = 0;
	coutSol.lengthViol = 0;
	int nbServices = 0;

	for (int kk = 1; kk <= paramsList[0]->nbDays; kk++)
	{
		// on parcourt les sommets grace aux resultats de split pour
		// remplir les structures
		j = (int)chromT[kk].size();

		for (int jj = 0; jj < paramsList[0]->nombreVehicules[kk]; jj++)
		{
			depot = paramsList[0]->ordreVehicules[kk][paramsList[0]->nombreVehicules[kk] - jj - 1].depotNumber;
			distance = 0;
			load = 0;
			i = (int)pred[kk][paramsList[0]->nombreVehicules[kk] - jj][j];
			chromR[kk][paramsList[0]->nombreVehicules[kk] - jj - 1] = i;

			if (j == i)
			{
				distance = 0;
				load = 0;
			}
			// case where there is only one delivery in the route
			else if (j == i + 1)
			{
				distance = paramsList[0]->timeCost[depot][chromT[kk][i]] + paramsList[0]->timeCost[chromT[kk][i]][depot];
				load = chromL[kk][chromT[kk][i]];
				nbServices++;
			}
			else
			{
				distance = paramsList[0]->timeCost[depot][chromT[kk][i]];
				load = 0;

				// infos sommets milieu
				for (int k = i; k <= j - 2; k++)
				{
					distance += paramsList[0]->timeCost[chromT[kk][k]][chromT[kk][k + 1]];
					load += chromL[kk][chromT[kk][k]];
					nbServices++;
				}

				// infos sommet fin
				distance += paramsList[0]->timeCost[chromT[kk][j - 1]][depot];
				load += chromL[kk][chromT[kk][j - 1]];
				nbServices++;
			}

			coutSol.fitness += distance;
			if (load > paramsList[0]->ordreVehicules[kk][paramsList[0]->nombreVehicules[kk] - jj - 1].capacity)
				coutSol.capacityViol += load - paramsList[0]->ordreVehicules[kk][paramsList[0]->nombreVehicules[kk] - jj - 1].capacity;

			j = i;
		}
	}

	// Add to the fitness the constants and the inventory cost
	if (paramsList[0]->isstockout)
	{
		vector<vector<double>> I_end(paramsList[0]->nbDays + 2, vector<double>(paramsList[0]->nbDepots + paramsList[0]->nbClients));
		for (int i = paramsList[0]->nbDepots; i < paramsList[0]->nbDepots + paramsList[0]->nbClients; i++)
		{
			I_end[0][i] = paramsList[0]->cli[i].startingInventory;
		}

		for (int k = 1; k <= paramsList[0]->nbDays; k++)
		{
			for (int cus = paramsList[0]->nbDepots; cus < paramsList[0]->nbDepots + paramsList[0]->nbClients; cus++)
			{
				coutSol.fitness += paramsList[0]->cli[cus].inventoryCost * std::max<double>(0, I_end[k - 1][cus] + chromL[k][cus] - paramsList[0]->cli[cus].dailyDemand[k]);
				coutSol.fitness += paramsList[0]->cli[cus].stockoutCost * std::max<double>(0, paramsList[0]->cli[cus].dailyDemand[k] - I_end[k - 1][cus] - chromL[k][cus]);

				coutSol.fitness -= chromL[k][cus] * (paramsList[0]->ancienNbDays + 1 - k) * paramsList[0]->inventoryCostSupplier;
				I_end[k][cus] = std::max<double>(0, I_end[k - 1][cus] + chromL[k][cus] - paramsList[0]->cli[cus].dailyDemand[k]);
			}
		}
	}
	// And the necessary constants
	coutSol.fitness += paramsList[0]->objectiveConstant_stockout;

	if (coutSol.capacityViol < 0.0001 && coutSol.lengthViol < 0.0001)
		estValide = true;
	else
		estValide = false;

	coutSol.evaluation = paramsList[0]->penalityCapa * coutSol.capacityViol + paramsList[0]->penalityLength * coutSol.lengthViol + coutSol.fitness;
}

// initialisation du vecteur potentiels
void Individu::initPot(int day)
{
	for (int i = 0; i < paramsList[0]->nombreVehicules[day] + 1; i++)
	{
		for (size_t j = 0; j <= chromT[day].size() + 1; j++)
		{
			potentiels[i][j] = 1.e30;
		}
	}
	potentiels[0][0] = 0;
}
// mise a jour de l'objet localSearch, attention, split doit avoir ete calcule avant
void Individu::updateLS()
{
	int depot;
	int i, j;
	Noeud *myDepot;
	Noeud *myDepotFin;
	Noeud *myClient;
	Route *myRoute;

	// We copy the amount of delivery per day
	// (more clean to make sure that LocalSearch can work totally independently of the Individu structure)
	localSearchList[0]->demandPerDay = chromL;

	for (int kk = 1; kk <= paramsList[0]->nbDays; kk++)
	{
		// on r�initialise l'ordreParcours
		localSearchList[0]->ordreParcours[kk].clear();

		// on replace les champs estPresent � false
		for (i = paramsList[0]->nbDepots; i < (int)localSearchList[0]->clients[kk].size(); i++)
		{
			localSearchList[0]->clients[kk][i]->estPresent = false;
		}

		// on parcourt les sommets grace aux resultats de split pour
		// remplir les structures
		j = (int)chromT[kk].size();

		for (int jj = 0; jj < paramsList[0]->nombreVehicules[kk]; jj++)
		{
			depot = paramsList[0]->ordreVehicules[kk][paramsList[0]->nombreVehicules[kk] - jj - 1].depotNumber;
			i = (int)pred[kk][paramsList[0]->nombreVehicules[kk] - jj][j];

			myDepot = localSearchList[0]->depots[kk][paramsList[0]->nombreVehicules[kk] - jj - 1];
			myDepotFin = localSearchList[0]->depotsFin[kk][paramsList[0]->nombreVehicules[kk] - jj - 1];
			myRoute = localSearchList[0]->routes[kk][paramsList[0]->nombreVehicules[kk] - jj - 1];

			myDepot->suiv = myDepotFin;
			myDepot->pred = myDepotFin;
			myDepotFin->suiv = myDepot;
			myDepotFin->pred = myDepot;

			// cas ou on a un seul sommet dans le cycle
			if (j == i + 1)
			{
				myClient = localSearchList[0]->clients[kk][chromT[kk][i]];
				myClient->pred = myDepot;
				myClient->suiv = myDepotFin;
				myClient->route = myRoute;
				myClient->estPresent = true;
				myDepot->suiv = myClient;
				myDepotFin->pred = myClient;
				localSearchList[0]->ordreParcours[kk].push_back(myClient->idx);
			}
			else if (j > i + 1)
			{
				// on a au moins 2 sommets
				// infos sommet debut
				myClient = localSearchList[0]->clients[kk][chromT[kk][i]];
				myClient->pred = myDepot;
				myClient->suiv = localSearchList[0]->clients[kk][chromT[kk][i + 1]];
				myClient->route = myRoute;
				myClient->estPresent = true;
				myDepot->suiv = myClient;
				localSearchList[0]->ordreParcours[kk].push_back(myClient->idx);

				// infos sommet fin
				myClient = localSearchList[0]->clients[kk][chromT[kk][j - 1]];
				myClient->pred = localSearchList[0]->clients[kk][chromT[kk][j - 2]];
				myClient->suiv = myDepotFin;
				myClient->route = myRoute;
				myClient->estPresent = true;
				myDepotFin->pred = myClient;
				localSearchList[0]->ordreParcours[kk].push_back(myClient->idx);

				// infos sommets milieu
				for (int k = (int)i + 1; k <= j - 2; k++)
				{
					myClient = localSearchList[0]->clients[kk][chromT[kk][k]];
					myClient->pred = localSearchList[0]->clients[kk][chromT[kk][k - 1]];
					myClient->suiv = localSearchList[0]->clients[kk][chromT[kk][k + 1]];
					myClient->route = myRoute;
					myClient->estPresent = true;
					localSearchList[0]->ordreParcours[kk].push_back(myClient->idx);
				}
			}
			j = i;
		}
		// pour chaque route on met les charges partielles � jour
		for (i = 0; i < (int)localSearchList[0]->routes[kk].size(); i++)
			localSearchList[0]->routes[kk][i]->updateRouteData();
	}
}
int partition(std::vector<Route *> &arr, int low, int high)
{
	Route *pivot = arr[high];
	int i = (low - 1);

	for (int j = low; j <= high - 1; j++)
	{
		if (arr[j]->centroidAngle <= pivot->centroidAngle)
		{
			i++;
			std::swap(arr[i], arr[j]);
		}
	}
	std::swap(arr[i + 1], arr[high]);
	return (i + 1);
}

// mise � jour du chromT suite aux modification de localSearch
void Individu::updateIndiv()
{
	// Don't forget to copy back the load delivered to each customer on each day
	chromL = localSearchList[0]->demandPerDay;

	vector<Route *> ordreRoutesAngle;
	Route *temp;
	Noeud *node;

	for (int kk = 1; kk <= paramsList[0]->nbDays; kk++)
	{
		ordreRoutesAngle = localSearchList[0]->routes[kk];

		for (int r = 0; r < (int)ordreRoutesAngle.size(); r++)
			ordreRoutesAngle[r]->updateCentroidCoord();

		// on recopie les noeuds dans le chromosome

		chromT[kk].clear();
		for (int r = 0; r < (int)ordreRoutesAngle.size(); r++)
		{
			node = ordreRoutesAngle[r]->depot->suiv;
			while (!node->estUnDepot)
			{
				chromT[kk].push_back(node->idx);
				node = node->suiv;
			}
		}
	}

	generalSplit();
}

// Computes the maximum amount of load that can be delivered to client on a day k without exceeding the
// customer maximum inventory
double Individu::maxFeasibleDeliveryQuantity(int day, int client)
{
	// Printing customer inventory and computing customer inventory cost
	double inventoryClient;
	double minSlack = 1.e30;
	inventoryClient = paramsList[0]->cli[client].startingInventory;
	for (int k = 1; k <= paramsList[0]->nbDays; k++)
	{
		// here level in the morning
		inventoryClient += chromL[k][client];
		//  level after receiving inventory

		// updating the residual capacity if k is greater than the day
		if (k >= day && paramsList[0]->cli[client].maxInventory - inventoryClient < minSlack)
			minSlack = paramsList[0]->cli[client].maxInventory - inventoryClient;

		if (minSlack < -0.0001)
		{
			cout << "ISSUE : negative slack capacit during crossover, this should not happen" << endl;
			throw("ISSUE : negative slack capacit during crossover, this should not happen");
		}

		inventoryClient -= paramsList[0]->cli[client].dailyDemand[k];
		// level after consumption
	}

	if (minSlack < 0.0001) // avoiding rounding issues
		return 0.;
	else
		return minSlack;
}

// distance g�n�rale
// il faut que les suivants de chaque individu aient �t� calcul�s
double Individu::distance(Individu *indiv2)
{
	int note = 0;
	bool isIdentical;

	// Inventory Routing
	// distance based on number of customers which have different service days

	for (int j = paramsList[0]->nbDepots; j < paramsList[0]->nbClients + paramsList[0]->nbDepots; j++)
	{
		isIdentical = true;
		for (int k = 1; k <= paramsList[0]->nbDays; k++)
			if ((chromL[k][j] < 0.0001 && indiv2->chromL[k][j] > 0.0001) || (indiv2->chromL[k][j] < 0.0001 && chromL[k][j] > 0.0001))
				isIdentical = false;
		if (isIdentical == false)
			note++;
	}

	return ((double)note / (double)(paramsList[0]->nbClients));
}

// ajoute un element proche dans les structures de proximit�
void Individu::addProche(Individu *indiv)
{
	list<proxData>::iterator it;
	proxData data;
	data.indiv = indiv;
	data.dist = distance(indiv);

	if (plusProches.empty())
		plusProches.push_back(data);
	else
	{
		it = plusProches.begin();
		while (it != plusProches.end() && it->dist < data.dist)
			++it;
		plusProches.insert(it, data);
	}
}

// enleve un element dans les structures de proximit�
void Individu::removeProche(Individu *indiv)
{
	list<proxData>::iterator last = plusProches.end();
	for (list<proxData>::iterator first = plusProches.begin(); first != last;)
		if (first->indiv == indiv)
			first = plusProches.erase(first);
		else
			++first;
}

// distance moyenne avec les n individus les plus proches
double Individu::distPlusProche(int n)
{
	double result = 0;
	double compte = 0;
	list<proxData>::iterator it = plusProches.begin();

	for (int i = 0; i < n && it != plusProches.end(); i++)
	{
		result += it->dist;
		compte += 1.0;
		++it;
	}
	return result / compte;
}
