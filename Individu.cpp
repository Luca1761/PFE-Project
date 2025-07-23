#include "Individu.h"

// constructeur d'un Individu comme simple conteneur

Individu::Individu(vector<Params*> pl) : paramsList(pl)
{
	nbScenario = paramsList.size();
	coutSol.evaluation = 0;
	coutSol.fitness = 0;
	coutSol.capacityViol = 0;

	//not same lengths because quantities at day 1 are second stage variables while roads are first stage
	chromT = vector<vector<int>>(nbScenario * (paramsList[0]->nbDays - 1) + 1 + 1);
	chromL = vector<vector<double>>(nbScenario * (paramsList[0]->nbDays) + 1, vector<double>(paramsList[0]->nbClients + paramsList[0]->nbDepots, 0.));

	// OPTION 2 -- JUST IN TIME POLICY //
	double dailyDelivery;
	Params *paramsTemp = paramsList[0];
	vector<vector<double>> startInventory;
	// DAY 1
	bool isFirstOption = (paramsList[0]->rng->genrand64_real1() < 0.1);
	for (unsigned int i = paramsTemp->nbDepots; i < paramsTemp->nbClients + paramsTemp->nbDepots; i++) {
		double initialInventory = paramsTemp->cli[i].startingInventory;
		vector<double> scenariosInventory(nbScenario);
		if (isFirstOption) {
			int nb = 0;
			for (int scenario = 0; scenario < nbScenario; scenario++) {
				nb += (initialInventory >= paramsList[scenario]->cli[i].dailyDemand[1]);
			}
			if (nb < nbScenario && paramsTemp->rng->genrand64_real1() < 0.5) {
				dailyDelivery = paramsTemp->cli[i].maxInventory - initialInventory;
				if (dailyDelivery > 0) chromT[1].push_back(i);
				for (int scenario = 0; scenario < nbScenario; scenario++) {
					chromL[1 + scenario * paramsTemp->nbDays][i] = dailyDelivery;
					scenariosInventory[scenario] = paramsTemp->cli[i].maxInventory - paramsList[scenario]->cli[i].dailyDemand[1];
				}
			} else {
				for (int scenario = 0; scenario < nbScenario; scenario++) {
					chromL[1 + scenario * paramsTemp->nbDays][i] = 0.0;
					scenariosInventory[scenario] = std::max(initialInventory - paramsList[scenario]->cli[i].dailyDemand[1], 0.0);
				}
			}
		} else {
			int nb = 0;
			for (int scenario = 0; scenario < nbScenario; scenario++) {
				nb += (initialInventory >= paramsList[scenario]->cli[i].dailyDemand[1]);
			}
			if (nb == nbScenario) {
				int nb1 = 0;
				for (int scenario = 0; scenario < nbScenario; scenario++) {
					double nextDayClientDemand = paramsList[scenario]->cli[i].dailyDemand[2];
					nb1 += (initialInventory - paramsList[scenario]->cli[i].dailyDemand[1]) >= nextDayClientDemand;
				}
				if (nb1 != nbScenario) {
					bool shouldDeliveryForNextDay = paramsList[0]->rng->genrand64_real1() <= 0.7;
					if (shouldDeliveryForNextDay) {
						chromT[1].push_back(i);
						for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
							double nextDayClientDemand = paramsList[scenario]->cli[i].dailyDemand[2];
							chromL[1 + scenario * paramsTemp->nbDays][i] = std::max(nextDayClientDemand - (initialInventory - paramsList[scenario]->cli[i].dailyDemand[1]), 1.0);
						}
					}
				}
				for (int scenario = 0; scenario < nbScenario; scenario++) {
					scenariosInventory[scenario] = initialInventory - paramsList[scenario]->cli[i].dailyDemand[1] + chromL[1 + scenario * paramsTemp->nbDays][i];
				}
			} else {
				chromT[1].push_back(i);
				for (int scenario = 0; scenario < nbScenario; scenario++) {
					chromL[1 + scenario * paramsTemp->nbDays][i] = std::max(paramsList[scenario]->cli[i].dailyDemand[1] - initialInventory, 1.0);
					scenariosInventory[scenario] = initialInventory - paramsList[scenario]->cli[i].dailyDemand[1] + chromL[1 + scenario * paramsTemp->nbDays][i];
				}
			}

		}
		startInventory.push_back(scenariosInventory);
	}
	//TO CHECK
	//OTHER DAYS
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		paramsTemp = paramsList[scenario];
		for (int i = paramsTemp->nbDepots; i < paramsTemp->nbClients + paramsTemp->nbDepots; i++) {
			double tempInventory = startInventory[i - paramsTemp->nbDepots][scenario];
			if (isFirstOption) {
				for (int day = 2; day <= paramsTemp->nbDays; day++) {
					bool doNotDeliver = (tempInventory >= paramsTemp->cli[i].dailyDemand[day] || paramsTemp->rng->genrand64_real1() < 0.5);
					if (doNotDeliver) {
						tempInventory = std::max<double>(0., tempInventory - paramsTemp->cli[i].dailyDemand[day]);
					} else {
						dailyDelivery = paramsTemp->cli[i].maxInventory - tempInventory; 
						tempInventory = paramsTemp->cli[i].maxInventory - paramsTemp->cli[i].dailyDemand[day];
						chromL[day + scenario * paramsTemp->nbDays][i] = dailyDelivery;
						if (dailyDelivery > 0)	chromT[day + scenario * (paramsTemp->nbDays - 1)].push_back(i);
					}
				}
			} else {
				for (int k = 2; k <= paramsTemp->nbDays; k++) {
					double currentDayClientDemand = paramsTemp->cli[i].dailyDemand[k];
					double nextDayClientDemand = paramsTemp->cli[i].dailyDemand[(k + 1) % paramsTemp->nbDays];

					if (tempInventory >= currentDayClientDemand) {
						// enough initial inventory, no need to service
						chromL[k + scenario * paramsTemp->nbDays][i] = 0;

						bool isInventoryEnoughForNextDay = tempInventory - currentDayClientDemand >= nextDayClientDemand;
						if (k < paramsTemp->nbDays && !isInventoryEnoughForNextDay) {
							bool shouldDeliveryForNextDay = paramsTemp->rng->genrand64_real1() <= 0.3;
							if (shouldDeliveryForNextDay) {
								chromL[k + scenario * paramsTemp->nbDays][i] = nextDayClientDemand - (tempInventory - currentDayClientDemand);
								chromT[k + scenario * (paramsTemp->nbDays - 1)].push_back(i);
							}
						}
						tempInventory = tempInventory + chromL[k + scenario * paramsTemp->nbDays][i] - currentDayClientDemand;
					} else {
						// not enough initial inventory, just in time policy for the initial solution
						double dailyDelivery = paramsTemp->cli[i].dailyDemand[k] - tempInventory;
						tempInventory = 0;
						chromL[k + scenario * paramsTemp->nbDays][i] = dailyDelivery;
						chromT[k + scenario * (paramsTemp->nbDays - 1)].push_back(i);
					}
				}
			}
		}
	}
	// And shuffle the whole solution
	for (int day = 1; day <= paramsList[0]->nbDays; day++) {
		for (unsigned int i = 0; i < chromT[day].size(); i++) {
			int j = i + paramsList[0]->rng->genrand64_int64() % ((int)chromT[day].size() - i); 
			// swap i and j elements
			std::swap(chromT[day][i], chromT[day][j]);
		}
	}

	// initialisation of the other structures
	pred = vector<vector<vector<int>>>(nbScenario * (paramsList[0]->nbDays - 1) + 1 + 1, vector<vector<int>>(paramsTemp->nombreVehicules[1] + 1, vector<int>(paramsTemp->nbClients + paramsTemp->nbDepots + 1, 0)));
	potentiels = vector<vector<double>>(paramsList[0]->nombreVehicules[1] + 1, vector<double>((int)paramsList[0]->nbClients + paramsList[0]->nbDepots + 1, 1.e30));
	potentiels[0][0] = 0;

	localSearchList = vector<LocalSearch*>(nbScenario, nullptr);
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
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
void Individu::generalSplit_scenario() {
	// lancement de la procedure split pour chaque jour
	// on essaye deja le split simple, si c'est sans succes , le split LF
	if (chromT[1].size() > 0) {
		// TO CHECK
		if (!splitSimpleDay1()) {
			std::cout << "ouUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUps" << std::endl;
		}
	}
	Params *paramsTemp;
	for (unsigned int k = 2; k < chromT.size(); k++) {
		int consideredScenario = (k-2) / (paramsList[0]->nbDays - 1);
		paramsTemp = paramsList[consideredScenario];
		if (chromT[k].size() > 0 && !splitSimple_scenario(k, paramsTemp, consideredScenario)) {
			splitLF_scenario(k, paramsTemp, consideredScenario);
		}
	}

	// After Split
	// we call a function that fills all other data structures and computes the cost
	measureSol_scenario();

	if (coutSol.evaluation > 1.e20) {
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			if (coutSol_scenario.capacityViol[scenario] > 0.0001)
				paramsList[scenario]->borneSplit *= 1.1;
		}
		generalSplit_scenario();
	}
}

// function split which does not consider a limit on the number of vehicles
// just uses the line "1" of the "potentiels" table.
int Individu::splitSimple_scenario(int k, Params *paramsTemp, int scenario)
{
	// on va utiliser la ligne 1 des potentiels et structures pred
	double load, distance, cost;
	int s0, s1, sb;
	potentiels[1][0] = 0;
	int day = k - scenario * (paramsTemp->nbDays - 1);
	int otherDay = scenario * (paramsList[0]->nbDays) + day;
	s0 = paramsTemp->ordreVehicules[day][0].depotNumber;
	for (int i = 0; i < (int)chromT[k].size(); i++)
	{
		load = 0;
		distance = 0;
		for (int j = i; j < (int)chromT[k].size() && load <= paramsTemp->ordreVehicules[day][0].capacity * paramsTemp->borneSplit; j++)
		{
			s1 = chromT[k][j];
			load += chromL[otherDay][s1];
			sb = (i == j) ? s0 : chromT[k][j - 1];
			distance += paramsTemp->timeCost[sb][s1];

			// computing the penalized cost
			cost = distance + paramsTemp->timeCost[s1][s0];
			if (load > paramsTemp->ordreVehicules[day][0].capacity)
				cost += (load - paramsTemp->ordreVehicules[day][0].capacity) * paramsTemp->penalityCapa;

			if (potentiels[1][i] + cost < potentiels[1][j + 1]) // basic Bellman algorithm
			{
				potentiels[1][j + 1] = potentiels[1][i] + cost;
				pred[k][1][j + 1] = i;
			}
		}
	}

	// testing if le number of vehicles is correct
	// in addition, the table pred is updated to keep track of everything
	int l = (int)chromT[k].size();
	for (int jj = 0; jj < paramsTemp->nombreVehicules[day]; jj++)
	{
		pred[k][paramsTemp->nombreVehicules[day] - jj][l] = pred[k][1][l];
		l = pred[k][1][l];
	}

	coutSol.evaluation = -1.e30; // just for security, making sure this value is not used (as it does not contain the constants)
	initPot_scenario(k, scenario);
	return (l == 0);
}

bool Individu::splitSimpleDay1() {
	// on va utiliser la ligne 1 des potentiels et structures pred
	double distance, averageCost;
	vector<double> loadScenario;
	int s0, s1, sb;
	potentiels[1][0] = 0;
	s0 = paramsList[0]->ordreVehicules[1][0].depotNumber;
	for (int i = 0; i < (int)chromT[1].size(); i++) {
		loadScenario = std::vector<double>(nbScenario, 0.0);
		distance = 0;
		bool isCapacityOk = true;
		for (int j = i; j < (int)chromT[1].size() && isCapacityOk; j++)
		{
			s1 = chromT[1][j];
			sb = (i == j) ? s0 : chromT[1][j - 1];
			distance += paramsList[0]->timeCost[sb][s1];
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {

				loadScenario[scenario] += chromL[scenario * (paramsList[0]->nbDays) + 1][s1];
				if (loadScenario[scenario] > paramsList[scenario]->ordreVehicules[1][0].capacity * paramsList[scenario]->borneSplit) isCapacityOk = false;
				if (loadScenario[scenario] > paramsList[scenario]->ordreVehicules[1][0].capacity)
					averageCost += (loadScenario[scenario] - paramsList[scenario]->ordreVehicules[1][0].capacity) * paramsList[scenario]->penalityCapa / nbScenario;
			}
			if (!isCapacityOk) break;
			// computing the penalized cost
			averageCost += distance + paramsList[0]->timeCost[s1][s0];

			if (potentiels[1][i] + averageCost < potentiels[1][j + 1]) // basic Bellman algorithm
			{
				potentiels[1][j + 1] = potentiels[1][i] + averageCost;
				pred[1][1][j + 1] = i;
			}
		}
	}

	// testing if the number of vehicles is correct
	// in addition, the table pred is updated to keep track of everything
	int l = (int)chromT[1].size();
	for (int jj = 0; jj < paramsList[0]->nombreVehicules[1]; jj++) {
		pred[1][paramsList[0]->nombreVehicules[1] - jj][l] = pred[1][1][l];
		l = pred[1][1][l];
	}

	coutSol.evaluation = -1.e30; // just for security, making sure this value is not used (as it does not contain the constants)
	initPot_scenario(1, 0);
	return (l == 0);
}

// fonction split pour probl�mes � flotte limit�e
void Individu::splitLF_scenario(int k, Params *paramsTemp, int scenario)
{
	double load, distance, cost;
	int sb, s0, s1, i, j;

	// 根据index和scenario序号获得 当前scenario的相对index
	int day = k - scenario * (paramsTemp->nbDays - 1);
	// pour chaque camion
	for (int cam = 0; cam < paramsTemp->nombreVehicules[day]; cam++) {
		i = 0;
		s0 = paramsTemp->ordreVehicules[day][cam].depotNumber;
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

			while (j < (int)chromT[k].size() && load <= paramsTemp->ordreVehicules[day][cam].capacity * paramsTemp->borneSplit)
			{
				s1 = chromT[k][j];
				load += chromL[k][s1];
				if (i == j)
				{
					distance = paramsTemp->timeCost[s0][s1];
				}
				else
				{
					sb = chromT[k][j - 1];
					distance += paramsTemp->timeCost[sb][s1];
				}

				// computing the penalized cost
				cost = distance + paramsTemp->timeCost[s1][s0];
				if (load > paramsTemp->ordreVehicules[day][cam].capacity)
					cost += (load - paramsTemp->ordreVehicules[day][cam].capacity) * paramsTemp->penalityCapa;

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

	// on nettoye ce que l'on a déplacé
	initPot_scenario(k, scenario);
}

//TO CHECK
void Individu::measureSol() {
	int depot;
	int i, j;
	double distance, load;
	double inventoryCost = 0.0;
	double routeCost = 0.0;
	double capaViol = 0.0;
	double fitness = 0.0;

	Params *paramsTemp = paramsList[0];
	vector<double> I_end(paramsTemp->nbDepots + paramsTemp->nbClients);
	for (unsigned int l = paramsTemp->nbDepots; l < paramsTemp->nbDepots + paramsTemp->nbClients; l++){
		I_end[l] = paramsTemp->cli[l].startingInventory;
	}

	inventoryCost += paramsTemp->availableSupply[1] * paramsTemp->inventoryCostSupplier;
	for (unsigned int cus = paramsTemp->nbDepots; cus < paramsTemp->nbDepots + paramsTemp->nbClients; cus++) {
		double toDeliver = 0.0;
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			toDeliver += chromL[1 + scenario * (paramsTemp->nbDays)][cus];
		}
		toDeliver /= nbScenario;
		
		inventoryCost += paramsTemp->cli[cus].inventoryCost * std::max<double>(0, I_end[cus] + toDeliver - paramsTemp->cli[cus].dailyDemand[1]);
		inventoryCost += paramsTemp->cli[cus].stockoutCost * std::max<double>(0, paramsTemp->cli[cus].dailyDemand[1] - I_end[cus] - toDeliver);
		inventoryCost -= toDeliver * paramsTemp->inventoryCostSupplier;
	}

	j = (int)chromT[1].size();
	for (unsigned int jj = 0; jj < paramsTemp->nombreVehicules[1]; jj++) {
		depot = paramsTemp->ordreVehicules[1][paramsTemp->nombreVehicules[1] - jj - 1].depotNumber;
		distance = 0;
		load = 0;
		i = (int)pred[1][paramsTemp->nombreVehicules[1] - jj][j];
		
		if (j == i) {
			distance = 0;
			load = 0;
		} else if (j == i + 1) {
			distance = paramsTemp->timeCost[depot][chromT[1][i]] + paramsTemp->timeCost[chromT[1][i]][depot];
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				load += chromL[1 + scenario * (paramsTemp->nbDays)][chromT[1][i]];
			}
			load /= nbScenario;
		} else {
			distance = paramsTemp->timeCost[depot][chromT[1][i]];
			load = 0;

			// infos sommets milieu
			for (int k = i; k <= j - 2; k++) {
				distance += paramsTemp->timeCost[chromT[1][k]][chromT[1][k + 1]];
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					load += chromL[1 + scenario * (paramsTemp->nbDays)][chromT[1][k]];
				}
			}

			// infos sommet fin
			distance += paramsTemp->timeCost[chromT[1][j - 1]][depot];
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++)
				load += chromL[1 + scenario * (paramsTemp->nbDays)][chromT[1][j - 1]];
			load /= nbScenario;
		}

		routeCost += distance;

		if (load > paramsTemp->ordreVehicules[1][paramsTemp->nombreVehicules[1] - jj - 1].capacity) {
			capaViol += load - paramsTemp->ordreVehicules[1][paramsTemp->nombreVehicules[1] - jj - 1].capacity;
		}
		j = i;	
	}

	fitness = routeCost + inventoryCost;

	std::cout << fitness << " " << capaViol << std::endl;
}

void Individu::measureSol_scenario() {
	int depot;
	int i, j;
	double distance, load;
	vector<double> inventoryCost(nbScenario, 0);
	vector<double> routeCost(nbScenario, 0);
	vector<double> capaViol(nbScenario, 0);
	vector<double> fitness(nbScenario, 0);
	coutSol_scenario.fitness = fitness;
	coutSol_scenario.evaluation = fitness;
	coutSol_scenario.capacityViol = capaViol;

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		Params *paramsTemp = paramsList[scenario];
		vector<vector<double>> I_end(paramsTemp->nbDays + 1, vector<double>(paramsTemp->nbDepots + paramsTemp->nbClients));
		for (unsigned int l = paramsTemp->nbDepots; l < paramsTemp->nbDepots + paramsTemp->nbClients; l++){
			I_end[0][l] = paramsTemp->cli[l].startingInventory;
		}

		int startIndexL = scenario * (paramsTemp->nbDays);

		vector<int> dayIndexL(paramsTemp->nbDays + 1, 0);
		for (int k = 1; k <= paramsTemp->nbDays; k++){
			dayIndexL[k] = startIndexL + k;
		}

		int startIndexT = scenario * (paramsTemp->nbDays - 1);

		vector<int> dayIndexT(paramsTemp->nbDays + 1, 0);
		dayIndexT[1] = 1;
		for (int k = 2; k <= paramsTemp->nbDays; k++){
			dayIndexT[k] = startIndexT + k;
		}

		for (unsigned int k = 1; k <= paramsTemp->nbDays; k++) {
			int day = dayIndexL[k];
			for (unsigned int cus = paramsTemp->nbDepots; cus < paramsTemp->nbDepots + paramsTemp->nbClients; cus++) {
				inventoryCost[scenario] += paramsTemp->cli[cus].inventoryCost * std::max<double>(0, I_end[k-1][cus]+chromL[day][cus]-paramsTemp->cli[cus].dailyDemand[k]);
				inventoryCost[scenario] += paramsTemp->cli[cus].stockoutCost * std::max<double>(0, paramsTemp->cli[cus].dailyDemand[k]-I_end[k-1][cus]-chromL[day][cus]);
				inventoryCost[scenario] -= chromL[day][cus] * (paramsTemp->ancienNbDays + 1 - k) * paramsTemp->inventoryCostSupplier;

				I_end[k][cus] = std::max<double>(0,I_end[k-1][cus] + chromL[day][cus] - paramsTemp->cli[cus].dailyDemand[k]);
			}
		}
		for (unsigned int kk = 1; kk <= paramsTemp->nbDays; kk++) {
			int dayT = dayIndexT[kk];
			int dayL = dayIndexL[kk];

			j = (int)chromT[dayT].size();
			for (unsigned int jj = 0; jj < paramsTemp->nombreVehicules[kk]; jj++) {
				depot = paramsTemp->ordreVehicules[kk][paramsTemp->nombreVehicules[kk] - jj - 1].depotNumber;
				distance = 0;
				load = 0;
				i = (int)pred[kk][paramsTemp->nombreVehicules[kk] - jj][j];

				if (j == i) {
					distance = 0;
					load = 0;
				} else if (j == i + 1) {
					distance = paramsTemp->timeCost[depot][chromT[dayT][i]] + paramsTemp->timeCost[chromT[dayT][i]][depot];
					load = chromL[dayL][chromT[dayT][i]];
				} else {
					distance = paramsTemp->timeCost[depot][chromT[dayT][i]];
					load = 0;

					// infos sommets milieu
					for (int k = i; k <= j - 2; k++) {
						distance += paramsTemp->timeCost[chromT[dayT][k]][chromT[dayT][k + 1]];
						load += chromL[dayL][chromT[dayT][k]];
					}

					// infos sommet fin
					distance += paramsTemp->timeCost[chromT[dayT][j - 1]][depot];
					load += chromL[dayL][chromT[dayT][j - 1]];
				}

				routeCost[scenario] += distance;

				if (load > paramsTemp->ordreVehicules[kk][paramsTemp->nombreVehicules[kk] - jj - 1].capacity) {
					capaViol[scenario] += load - paramsTemp->ordreVehicules[kk][paramsTemp->nombreVehicules[kk] - jj - 1].capacity;
				}
				j = i;	
			}
		}

		fitness[scenario] = routeCost[scenario] + inventoryCost[scenario] + paramsTemp->objectiveConstant_stockout;
	}

	estValide = true;
	coutSol.fitness = 0;
	coutSol.evaluation = 0;
	coutSol.capacityViol = 0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		coutSol_scenario.fitness[scenario] = fitness[scenario];
		coutSol.fitness += fitness[scenario];
		coutSol_scenario.evaluation[scenario] = paramsList[scenario]->penalityCapa * capaViol[scenario] + fitness[scenario];
		coutSol.evaluation += coutSol_scenario.evaluation[scenario];
		coutSol_scenario.capacityViol[scenario] = capaViol[scenario];
		if (capaViol[scenario] > 0.0001) {
			coutSol.capacityViol += capaViol[scenario];
			estValide = false;
		}
	}
	coutSol.evaluation /= nbScenario;
	coutSol.fitness /= nbScenario;
	coutSol.capacityViol /= nbScenario;
}

// initialisation du vecteur potentiels
void Individu::initPot_scenario(int k, int scenario)
{
	int day = k - scenario * (paramsList[0]->nbDays - 1);
	for (int i = 0; i < paramsList[0]->nombreVehicules[day] + 1; i++)
	{
		for (size_t j = 0; j <= chromT[k].size() + 1; j++)
		{
			potentiels[i][j] = 1.e30;
		}
	}
	potentiels[0][0] = 0;
}

void Individu::updateLS_scenario() {
	int i, j;
	vector<Noeud*> myDepot(nbScenario);
	vector<Noeud*> myDepotFin(nbScenario);
	vector<Noeud*> myClient(nbScenario);
	vector<Route*> myRoute(nbScenario);

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		Noeud * tempDepot = myDepot[scenario];
		Noeud * tempDepotFin = myDepotFin[scenario];
		Noeud * tempClient = myClient[scenario];
		Route * tempRoute = myRoute[scenario];

		Params *paramsTemp = paramsList[scenario];
		unsigned int startIndex = scenario * (paramsTemp->nbDays);
		unsigned int startIndexT = scenario * (paramsTemp->nbDays - 1);
		vector<vector<double> > deliveryPerDay(paramsTemp->nbDays+1, vector<double>(paramsTemp->nbClients + paramsTemp->nbDepots, 0));

		for (unsigned int day = 1; day <= paramsTemp->nbDays; day++) {
			for (unsigned int client = 0; client < paramsTemp->nbClients + paramsTemp->nbDepots; client++) {
				deliveryPerDay[day][client] = chromL[startIndex + day][client];
			}
		}
		localSearchList[scenario]->deliveryPerDay = deliveryPerDay;

		for (unsigned int kk = 1; kk <= paramsTemp->nbDays; kk++) {
			localSearchList[scenario]->ordreParcours[kk].clear();
			for (unsigned int l = paramsTemp->nbDepots; l < (int)localSearchList[scenario]->clients[kk].size(); l++) {
				localSearchList[scenario]->clients[kk][l]->estPresent = false;
			}

			int chromIndex = (kk == 1) ? 1 : startIndexT + kk;
			j = (int)chromT[chromIndex].size();

			for (unsigned int jj = 0; jj < paramsTemp->nombreVehicules[kk]; jj++) {
				i = (int)pred[kk][paramsTemp->nombreVehicules[kk] - jj][j];

				tempDepot = localSearchList[scenario]->depots[kk][paramsTemp->nombreVehicules[kk] - jj - 1];
				tempDepotFin = localSearchList[scenario]->depotsFin[kk][paramsTemp->nombreVehicules[kk] - jj - 1];
				tempRoute = localSearchList[scenario]->routes[kk][paramsTemp->nombreVehicules[kk] - jj - 1];

				tempDepot->suiv = tempDepotFin;
				tempDepot->pred = tempDepotFin;
				tempDepotFin->suiv = tempDepot;
				tempDepotFin->pred = tempDepot;

				// cas ou on a un seul sommet dans le cycle
				if (j == i + 1) {
					tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][i]];
					tempClient->pred = tempDepot;
					tempClient->suiv = tempDepotFin;
					tempClient->route = tempRoute;
					tempClient->estPresent = true;
					tempDepot->suiv = tempClient;
					tempDepotFin->pred = tempClient;
					localSearchList[scenario]->ordreParcours[kk].push_back(tempClient->idx);
				}
				else if (j > i + 1) {
					// on a au moins 2 sommets
					// infos sommet debut
					tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][i]];
					tempClient->pred = tempDepot;
					tempClient->suiv = localSearchList[scenario]->clients[kk][chromT[chromIndex][i + 1]];
					tempClient->route = tempRoute;
					tempClient->estPresent = true;
					tempDepot->suiv = tempClient;
					localSearchList[scenario]->ordreParcours[kk].push_back(tempClient->idx);

					// infos sommet fin
					tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][j - 1]];
					tempClient->pred = localSearchList[scenario]->clients[kk][chromT[chromIndex][j - 2]];
					tempClient->suiv = tempDepotFin;
					tempClient->route = tempRoute;
					tempClient->estPresent = true;
					tempDepotFin->pred = tempClient;
					localSearchList[scenario]->ordreParcours[kk].push_back(tempClient->idx);

					// infos sommets milieu
					for (int k = (int)i + 1; k <= j - 2; k++)
					{
						tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][k]];
						tempClient->pred = localSearchList[scenario]->clients[kk][chromT[chromIndex][k - 1]];
						tempClient->suiv = localSearchList[scenario]->clients[kk][chromT[chromIndex][k + 1]];
						tempClient->route = tempRoute;
						tempClient->estPresent = true;
						localSearchList[scenario]->ordreParcours[kk].push_back(tempClient->idx);
					}
				}
				j = i;
			}
			// pour chaque route on met les charges partielles à jour
			for (unsigned i = 0; i < localSearchList[scenario]->routes[kk].size(); i++)
				localSearchList[scenario]->routes[kk][i]->updateRouteData();
		}
	}
}

void Individu::localSearchRunSearch_scenario() {
	const int GROUP_SIZE = 1 + nbScenario / 4;
	
	// Local search moves (mutation1-mutation9)
	vector<thread> threads;
	for (int i = 0; i < nbScenario; i += GROUP_SIZE) {
		int end = min(i + GROUP_SIZE, nbScenario);
		threads.emplace_back([this, i, end]() {
			for (unsigned int scenario = i; scenario < end; scenario++) {
				localSearchList[scenario]->runSearchSameDay();
			}
		});
	}
	for (auto& t : threads) {
		t.join();
	}

	// mutation1-mutation9 for day 1 (average cost reduction)
	runSearchDay1();

	// Our brand new operator of local search using dynamic programming
	muterDifferentScenarioDP();
	
	// We repeat the process after dynamic programming
	vector<thread> threads2;
	for (int i = 0; i < nbScenario; i += GROUP_SIZE) {
		int end = min(i + GROUP_SIZE, nbScenario);
		threads2.emplace_back([this, i, end]() {
			for (unsigned int scenario = i; scenario < end; scenario++) {
				localSearchList[scenario]->runSearchSameDay();
			}
		});
	}
	
	for (auto& t : threads2) {
		t.join();
	}
	runSearchDay1();
}

void Individu::runSearchDay1() {
	int nbMoves = 1;
	int nbPhases = 0;
	while (nbMoves > 0 && nbPhases < 100) {
		nbMoves = 0;
		nbMoves += mutationSameDay1();
		nbPhases++;
	}
}

int Individu::mutationSameDay1() {
	localSearchList[0]->dayCour = 1;
	int size = (int)localSearchList[0]->ordreParcours[1].size();
	int size2;
	localSearchList[0]->rechercheTerminee = false;
	bool moveEffectue = false;
	int nbMoves = 0;
	localSearchList[0]->firstLoop = true;
	
	while (!localSearchList[0]->rechercheTerminee) {
		localSearchList[0]->rechercheTerminee = true;
		moveEffectue = false;
		for (int posU = 0; posU < size; posU++) {
			posU -= moveEffectue; // on retourne sur le dernier noeud si on a modifie
			nbMoves += moveEffectue;
			moveEffectue = false;
			
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				localSearchList[scenario]->noeudU = localSearchList[scenario]->clients[1][localSearchList[0]->ordreParcours[1][posU]];
				localSearchList[scenario]->noeudUPred = localSearchList[scenario]->noeudU->pred;
				localSearchList[scenario]->x = localSearchList[scenario]->noeudU->suiv;
				localSearchList[scenario]->noeudXSuiv = localSearchList[scenario]->x->suiv;
				localSearchList[scenario]->xSuivCour = localSearchList[scenario]->x->suiv->idx;
				localSearchList[scenario]->routeU = localSearchList[scenario]->noeudU->route;
				localSearchList[scenario]->noeudUCour = localSearchList[scenario]->noeudU->idx;
				localSearchList[scenario]->noeudUPredCour = localSearchList[scenario]->noeudUPred->idx;
				localSearchList[scenario]->xCour = localSearchList[scenario]->x->idx;
			}

			size2 = (int)localSearchList[0]->noeudU->moves.size();
			for (int posV = 0; posV < size2 && moveEffectue == 0; posV++) {
				for (int scenario = 0; scenario < nbScenario; scenario++) 
				localSearchList[scenario]->noeudV = localSearchList[scenario]->clients[1][localSearchList[0]->noeudU->moves[posV]];
				if (!localSearchList[0]->noeudV->route->nodeAndRouteTested[localSearchList[0]->noeudU->idx] ||
					!localSearchList[0]->noeudU->route->nodeAndRouteTested[localSearchList[0]->noeudU->idx] || localSearchList[0]->firstLoop)
					{
						for (int scenario = 0; scenario < nbScenario; scenario++) {
							localSearchList[scenario]->noeudVPred = localSearchList[scenario]->noeudV->pred;
							localSearchList[scenario]->y = localSearchList[scenario]->noeudV->suiv;
							localSearchList[scenario]->noeudYSuiv = localSearchList[scenario]->y->suiv;
							localSearchList[scenario]->ySuivCour = localSearchList[scenario]->y->suiv->idx;
							localSearchList[scenario]->routeV = localSearchList[scenario]->noeudV->route;
							localSearchList[scenario]->noeudVCour = localSearchList[scenario]->noeudV->idx;
							localSearchList[scenario]->noeudVPredCour = localSearchList[scenario]->noeudVPred->idx;
							localSearchList[scenario]->yCour = localSearchList[scenario]->y->idx;
						}
						
						if (!moveEffectue)
						moveEffectue = mutation1_indiv();
						if (!moveEffectue)
						moveEffectue = mutation2_indiv();
						if (!moveEffectue)
						moveEffectue = mutation3_indiv();
						
						// les mutations 4 et 6 (switch) , sont sym�triques
						if (localSearchList[0]->noeudU->idx <= localSearchList[0]->noeudV->idx) {
							if (!moveEffectue)
							moveEffectue = mutation4_indiv();
							if (!moveEffectue)
							moveEffectue = mutation6_indiv();
						}
						if (!moveEffectue)
						moveEffectue = mutation5_indiv();
						
						// mutations 2-opt
						if (!moveEffectue)
						moveEffectue = mutation7_indiv();
						if (!moveEffectue)
						moveEffectue = mutation8_indiv();
						if (!moveEffectue)
						moveEffectue = mutation9_indiv();
						
						if (moveEffectue) {
							for (int scenario = 0; scenario < nbScenario; scenario++) {
								localSearchList[scenario]->routeU->reinitSingleDayMoves();
								localSearchList[scenario]->routeV->reinitSingleDayMoves();
							}
						}
					}
			}
			
			// c'est un depot on tente l'insertion derriere le depot de ce jour
			// si il ya correlation
			if (localSearchList[0]->params->isCorrelated1[localSearchList[0]->noeudU->idx][localSearchList[0]->depots[1][0]->idx] &&
				!moveEffectue)
				for (int route = 0; route < (int)localSearchList[0]->depots[1].size(); route++)
				{
					for (int scenario = 0; scenario < nbScenario; scenario++) {
						localSearchList[scenario]->noeudV = localSearchList[scenario]->depots[1][route];
					}
					if (!localSearchList[0]->noeudV->route->nodeAndRouteTested[localSearchList[0]->noeudU->idx] ||
						!localSearchList[0]->noeudU->route->nodeAndRouteTested[localSearchList[0]->noeudU->idx] || localSearchList[0]->firstLoop)
					{
						for (int scenario = 0; scenario < nbScenario; scenario++) {
							localSearchList[scenario]->noeudVPred = localSearchList[scenario]->noeudV->pred;
							localSearchList[scenario]->y = localSearchList[scenario]->noeudV->suiv;
							localSearchList[scenario]->noeudYSuiv = localSearchList[scenario]->y->suiv;
							localSearchList[scenario]->ySuivCour = localSearchList[scenario]->y->suiv->idx;
							localSearchList[scenario]->routeV = localSearchList[scenario]->noeudV->route;
							localSearchList[scenario]->noeudVCour = localSearchList[scenario]->noeudV->idx;
							localSearchList[scenario]->noeudVPredCour = localSearchList[scenario]->noeudVPred->idx;
							localSearchList[scenario]->yCour = localSearchList[scenario]->y->idx;
						}

						if (!moveEffectue)
						moveEffectue = mutation1_indiv();
						if (!moveEffectue)
						moveEffectue = mutation2_indiv();
						if (!moveEffectue)
						moveEffectue = mutation3_indiv();

						if (!localSearchList[0]->noeudV->suiv->estUnDepot)
						{
						if (!moveEffectue)
							moveEffectue = mutation8_indiv();
						if (!moveEffectue)
							moveEffectue = mutation9_indiv();
						}

						if (moveEffectue) {
							for (int scenario = 0; scenario < nbScenario; scenario++) {
								localSearchList[scenario]->routeU->reinitSingleDayMoves();
								localSearchList[scenario]->routeV->reinitSingleDayMoves();
							}
						}
					}
				}
		// TODO -- check that memories are working
		}
		localSearchList[0]->firstLoop = false;
	}
	return nbMoves;
}

void Individu::muterDifferentScenarioDP() {
	vector<int> randomClients;

	for (int client = paramsList[0]->nbDepots; client < paramsList[0]->nbClients + paramsList[0]->nbDepots; client++) {
		randomClients.push_back(client);
	}

	std::mt19937 g(paramsList[0]->seed);
	shuffle(randomClients.begin(), randomClients.end(), g);

	bool rechercheTerminee = false;
   	int nbMoves = 0;
	int nbPhases = 0;
	while(!rechercheTerminee){
		rechercheTerminee = true;		
		for (int client : randomClients) {
			nbMoves += mutationDP(client, rechercheTerminee);
		}
		nbPhases++;
  	}
}

int Individu::mutationDP(int client, bool &rechercheTerminee) {
	Noeud *noeudTravail;
	double currentCost = 0.0;
	// First, make sure all insertion costs are computed
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		for (int k = 1; k <= paramsList[0]->ancienNbDays; k++){
			noeudTravail = localSearchList[scenario]->clients[k][client]; //node* day k client
			localSearchList[scenario]->computeCoutInsertion(noeudTravail); // detour,place (dominated) for each route
		}
	}
	// Compute the current lot sizing solution cost (from the model point of view)
	//before optimization
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		currentCost += localSearchList[scenario]->evaluateCurrentCost_stockout(client);
	}
	currentCost /= nbScenario;
	/* Generate the structures of the subproblem */
	vector<vector<vector<Insertion>>> insertions = vector<vector<vector<Insertion>>>(nbScenario, vector<vector<Insertion>>(paramsList[0]->nbDays));
	
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		for (int k = 1; k <= paramsList[0]->nbDays; k++)
		{
			insertions[scenario][k - 1] = localSearchList[scenario]->clients[k][client]->allInsertions;
		}
	}
	
	unique_ptr<LotSizingSolver> lotsizingSolver(
			make_unique<LotSizingSolver>(paramsList, insertions, client));

		
	lotsizingSolver->solveStockoutBackward();
		
	vector<vector<double>> quantities = vector<vector<double>>(nbScenario, vector<double>(paramsList[0]->nbDays));
	vector<vector<int>> breakpoints = vector<vector<int>>(nbScenario, vector<int>(paramsList[0]->nbDays));
	vector<double> objectiveScenarios = lotsizingSolver->objective;
	double objective = 0.0;
	for (int scenario = 0; scenario < nbScenario; scenario++){
		objective += objectiveScenarios[scenario];
	}
	objective /= (double) nbScenario;
	quantities = lotsizingSolver->quantities;
	
	if(lt(currentCost,objective-0.01)) {
		return 0;
	}

	/* APPLYING THE MOVEMENT */
	// Later on we will verify whether it's an improving move or not to trigger a
	// good termination.

	// First, removing all occurences of the node.
	double test = 0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		for (int k = 1; k <= paramsList[scenario]->ancienNbDays; k++) {
			noeudTravail = localSearchList[scenario]->clients[k][client];
			if (noeudTravail->estPresent){
				localSearchList[scenario]->removeNoeud(noeudTravail);
			}
			localSearchList[scenario]->deliveryPerDay[k][client] = 0.;

		}
		// Then looking at the solution of the model and inserting in the good place
		for (int k = 1; k <= paramsList[scenario]->ancienNbDays; k++)
		{
			if (quantities[scenario][k - 1] > 0.0001 || (lotsizingSolver->breakpoints[scenario][k - 1]&&gt(0,lotsizingSolver->breakpoints[scenario][k - 1]->detour) )) // don't forget that in the model the index      // goes from 0 to t-1
			{
			localSearchList[scenario]->deliveryPerDay[k][client] = round(quantities[scenario][k - 1]);
			
			localSearchList[scenario]->clients[k][client]->placeInsertion = lotsizingSolver->breakpoints[scenario][k - 1]->place;
		
			localSearchList[scenario]->addNoeud(localSearchList[scenario]->clients[k][client]);
			}
		}

		double realCost = localSearchList[scenario]->evaluateCurrentCost_stockout(client);
		test+=realCost;
		if (fabs(realCost- objectiveScenarios[scenario])>0.001) {
			std::cout << "The solution doesn't give the expected cost for scenario " << scenario << std::endl;
			std::cout << "Cost: " << realCost << "; Expected cost: " << objectiveScenarios[scenario] << std::endl;
			for (int scenario1 = 0; scenario1 < nbScenario; scenario1++) {
				std::cout << round(quantities[scenario][0]) << std::endl;
				std::cout << quantities[scenario][0] << std::endl;
			}
			throw string("Cost error");
			return 0;
		}
	}

	if (lt(objective + 0.01, currentCost))// An improving move has been found, the search is not finished.
	{
		rechercheTerminee = false;
		return 1;
	}
	else
		return 0;
}

// mise a jour du chromT et chromL suite aux modification de localSearch
void Individu::updateIndiv_scenario() {
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		Params *paramsTemp = paramsList[scenario];
		int startIndex = scenario * (paramsTemp->nbDays) + 1;
		int endIndex = startIndex + (paramsTemp->nbDays);
		for (int i = startIndex; i < endIndex; i++) {  
			if (localSearchList[scenario] != nullptr && i - startIndex + 1 < localSearchList[scenario]->deliveryPerDay.size()) {
				chromL[i] = localSearchList[scenario]->deliveryPerDay[i - startIndex + 1];
			}
		}
		Noeud *node;
		int startIndexT = scenario * (paramsTemp->nbDays - 1);
		for (unsigned int day = 1; day <= paramsTemp->nbDays; day++) {
			int chromIndex = day == 1 ? 1 : startIndexT + day;
			chromT[chromIndex].clear();
			for (Route* temp : localSearchList[scenario]->routes[day]) {
				node = temp->depot->suiv;
				while (!node->estUnDepot) {
					chromT[chromIndex].push_back(node->idx);
					node = node->suiv;
				}
			}
		}
	}
	generalSplit_scenario();
}

// distance generale
double Individu::distance(Individu *indiv2) {
	// TO CHECK
	double note = 0;
	bool isIdentical;
	vector<int> dayIndexL(paramsList[0]->nbDays + 1, 0);

	// Inventory Routing
	// distance based on number of customers which have different service days
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		int noteScenario = 0;
		int startIndex = scenario * (paramsList[0]->nbDays);
		for (int k = 1; k <= paramsList[0]->nbDays; k++)
			dayIndexL[k] = startIndex + k;
		for (int client = paramsList[0]->nbDepots; client < paramsList[0]->nbClients + paramsList[0]->nbDepots; client++) {
			isIdentical = true;
			for (int k : dayIndexL)
				if ((chromL[k][client] < 0.0001 && indiv2->chromL[k][client] > 0.0001) || (indiv2->chromL[k][client] < 0.0001 && chromL[k][client] > 0.0001))
					isIdentical = false;
			if (isIdentical == false)
				noteScenario++;
		}
		note += ((double)noteScenario / (double)(paramsList[0]->nbClients));
	}

	return ((double)note / (double)(nbScenario));
}

// ajoute un element proche dans les structures de proximite
void Individu::addProche(Individu *indiv) {
	list<proxData>::iterator it;
	proxData data;
	data.indiv = indiv;
	data.dist = distance(indiv);

	if (plusProches.empty()) {
		plusProches.push_back(data);
	} else {
		it = plusProches.begin();
		while (it != plusProches.end() && it->dist < data.dist) {
			++it;
		}
		plusProches.insert(it, data);
	}
}

// enleve un element dans les structures de proximite
void Individu::removeProche(Individu *indiv) {
	list<proxData>::iterator last = plusProches.end();
	for (list<proxData>::iterator it = plusProches.begin(); it != last;)
		if (it->indiv == indiv) it = plusProches.erase(it);
		else ++it;
}

// distance moyenne avec les n individus les plus proches
double Individu::distPlusProche(int n) {
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

int Individu::mutation1_indiv() {
	double costSuppU = paramsList[0]->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->xCour] 
	- paramsList[0]->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->noeudUCour]  
	- paramsList[0]->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->xCour];

	double costSuppV = paramsList[0]->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->noeudUCour] 
	+ paramsList[0]->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->yCour] 
	- paramsList[0]->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->yCour];

	// dans le cas ou l'on est dans la meme route , le cout n'est pas calcule correctement en realite
	// tout ce qu'on sait c'est que si il est negatif c'est qu'il est bien reellement negatif
	// pas d'incidence pour l'instant mais attention
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV) {
			costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour])*paramsList[scenario]->penalityCapa
			- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge)*paramsList[scenario]->penalityCapa) / (double)nbScenario ;

			costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour])*paramsList[scenario]->penalityCapa
			- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge)*paramsList[scenario]->penalityCapa) / (double)nbScenario ;
		}
	}

	if (costSuppU + costSuppV > -0.0001) return 0 ;
	if (localSearchList[0]->noeudUCour == localSearchList[0]->yCour) return 0;

	// mettre a jour les noeuds
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->noeudU,localSearchList[scenario]->noeudV);
	}

	localSearchList[0]->rechercheTerminee = false;
	return 1 ;
}

int Individu::mutation2_indiv() {
	double costSuppU = paramsList[0]->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->xSuivCour] 
	- paramsList[0]->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->noeudUCour] 
	- paramsList[0]->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->xCour] 
	- paramsList[0]->timeCost[localSearchList[0]->xCour][localSearchList[0]->xSuivCour];

	double costSuppV = paramsList[0]->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->noeudUCour] 
	+ paramsList[0]->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->xCour] 
	+ paramsList[0]->timeCost[localSearchList[0]->xCour][localSearchList[0]->yCour] 
	- paramsList[0]->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->yCour];

	for (int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV) {
			costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*paramsList[scenario]->penalityCapa
			- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;
			
			costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*paramsList[scenario]->penalityCapa
			- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( localSearchList[0]->noeudU == localSearchList[0]->y || localSearchList[0]->noeudV == localSearchList[0]->x || localSearchList[0]->x->estUnDepot ) { return 0 ;}

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->noeudU,localSearchList[scenario]->noeudV);
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->x,localSearchList[scenario]->noeudU);
	}

	localSearchList[0]->rechercheTerminee = false ; 
	return 1 ;
}

int Individu::mutation3_indiv() {
	Params* paramsTemp = paramsList[0];
	double costSuppU = paramsTemp->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->xSuivCour] 
	- paramsTemp->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->noeudUCour] 
	- paramsTemp->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->xCour] 
	- paramsTemp->timeCost[localSearchList[0]->xCour][localSearchList[0]->xSuivCour];

	double costSuppV = paramsTemp->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->xCour] 
	+ paramsTemp->timeCost[localSearchList[0]->xCour][localSearchList[0]->noeudUCour] 
	+ paramsTemp->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->yCour] 
	- paramsTemp->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->yCour];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*paramsList[scenario]->penalityCapa
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;
		
		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*paramsList[scenario]->penalityCapa
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) return 0;
	if (localSearchList[0]->noeudU == localSearchList[0]->y ||  localSearchList[0]->x == localSearchList[0]->noeudV || localSearchList[0]->x->estUnDepot ) return 0;

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->x,localSearchList[scenario]->noeudV);
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->noeudU,localSearchList[scenario]->x);
	}

	localSearchList[0]->rechercheTerminee = false ; 
	return 1 ;
}

int Individu::mutation4_indiv() {
	Params* paramsTemp = paramsList[0];
	LocalSearch* local = localSearchList[0];
	double costSuppU = paramsTemp->timeCost[local->noeudUPredCour][local->noeudVCour] 
	+ paramsTemp->timeCost[local->noeudVCour][local->xCour]
	- paramsTemp->timeCost[local->noeudUPredCour][local->noeudUCour] 
	- paramsTemp->timeCost[local->noeudUCour][local->xCour];

	double costSuppV = paramsTemp->timeCost[local->noeudVPredCour][local->noeudUCour] 
	+ paramsTemp->timeCost[local->noeudUCour][local->yCour]
	- paramsTemp->timeCost[local->noeudVPredCour][local->noeudVCour] 
	- paramsTemp->timeCost[local->noeudVCour][local->yCour];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour])*paramsList[scenario]->penalityCapa
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;
		
		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour])*paramsList[scenario]->penalityCapa
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( local->noeudUCour == local->noeudVPredCour || local->noeudUCour == local->yCour) { return 0 ;}

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->swapNoeud(localSearchList[scenario]->noeudU, localSearchList[scenario]->noeudV) ;
	}

	localSearchList[0]->rechercheTerminee = false ; 
	return 1 ;
}

int Individu::mutation5_indiv() {
	Params* paramsTemp = paramsList[0];
	LocalSearch* local = localSearchList[0];
	double costSuppU = paramsTemp->timeCost[local->noeudUPredCour][local->noeudVCour] 
	+ paramsTemp->timeCost[local->noeudVCour][local->xSuivCour]
	- paramsTemp->timeCost[local->noeudUPredCour][local->noeudUCour] 
	- paramsTemp->timeCost[local->noeudUCour][local->xCour] 
	- paramsTemp->timeCost[local->xCour][local->xSuivCour];

	double costSuppV = paramsTemp->timeCost[local->noeudVPredCour][local->noeudUCour] 
	+ paramsTemp->timeCost[local->xCour][local->yCour]
	+ paramsTemp->timeCost[local->noeudUCour][local->xCour]
	- paramsTemp->timeCost[local->noeudVPred->idx][local->noeudVCour] 
	- paramsTemp->timeCost[local->noeudVCour][local->yCour];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*paramsList[scenario]->penalityCapa
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;

		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour])*paramsList[scenario]->penalityCapa
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( local->noeudU == local->noeudVPred || local->x == local->noeudVPred || local->noeudU == local->y || local->x->estUnDepot ) { return 0 ;}

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->swapNoeud(localSearchList[scenario]->noeudU, localSearchList[scenario]->noeudV) ;
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->x, localSearchList[scenario]->noeudU);
	}
	local->rechercheTerminee = false ; 
	return 1 ;
}

int Individu::mutation6_indiv() {
	Params* paramsTemp = paramsList[0];
	LocalSearch* local = localSearchList[0];
	double costSuppU = paramsTemp->timeCost[local->noeudUPredCour][local->noeudVCour]  
	+ paramsTemp->timeCost[local->noeudVCour][local->yCour]
	+ paramsTemp->timeCost[local->yCour][local->xSuivCour]
	- paramsTemp->timeCost[local->noeudUPredCour][local->noeudUCour] 
	- paramsTemp->timeCost[local->noeudUCour][local->xCour] 
	- paramsTemp->timeCost[local->xCour][local->xSuivCour];

	double costSuppV = paramsTemp->timeCost[local->noeudVPredCour][local->noeudUCour] 
	+ paramsTemp->timeCost[local->noeudUCour][local->xCour]
	+ paramsTemp->timeCost[local->xCour][local->ySuivCour]
	- paramsTemp->timeCost[local->noeudVPredCour][local->noeudVCour] 
	- paramsTemp->timeCost[local->noeudVCour][local->yCour]
	- paramsTemp->timeCost[local->yCour][local->ySuivCour];
	
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->yCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*paramsList[scenario]->penalityCapa
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;
		
		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->yCour])*paramsList[scenario]->penalityCapa
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( local->x->estUnDepot || local->y->estUnDepot || local->y == local->noeudUPred || local->noeudU == local->y || local->x == local->noeudV || local->noeudV == local->noeudXSuiv ) { return 0 ;}

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->swapNoeud(localSearchList[scenario]->noeudU, localSearchList[scenario]->noeudV) ;
		localSearchList[scenario]->swapNoeud(localSearchList[scenario]->x,localSearchList[scenario]->y) ;
	}

	local->rechercheTerminee = false ; 
	return 1 ;
	
}

int Individu::mutation7_indiv() {
	Params* paramsTemp = paramsList[0];
	LocalSearch* local = localSearchList[0];
	if  ((local->routeU->idx != local->routeV->idx) || local->noeudU->suiv == local->noeudV || local->noeudU->place > local->noeudV->place ) {  return 0 ; }
	
	double cost = paramsTemp->timeCost[local->noeudUCour][local->noeudVCour] + paramsTemp->timeCost[local->xCour][local->yCour]
	- paramsTemp->timeCost[local->noeudUCour][local->xCour] - paramsTemp->timeCost[local->noeudVCour][local->yCour] ;
	
	if ( cost > -0.0001 ) { return 0 ;}
	
	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		Noeud * nodeNum = localSearchList[scenario]->noeudXSuiv ;
		Noeud * temp ;
		localSearchList[scenario]->x->pred = nodeNum ;
		localSearchList[scenario]->x->suiv = localSearchList[scenario]->y ;

		while ( nodeNum != localSearchList[scenario]->noeudV )
		{
			temp = nodeNum->suiv ;
			nodeNum->suiv = nodeNum->pred ;
			nodeNum->pred = temp ;
			nodeNum = temp ;
		}

		localSearchList[scenario]->noeudV->suiv = localSearchList[scenario]->noeudV->pred ;
		localSearchList[scenario]->noeudV->pred = localSearchList[scenario]->noeudU ;
		localSearchList[scenario]->noeudU->suiv = localSearchList[scenario]->noeudV ;
		localSearchList[scenario]->y->pred = localSearchList[scenario]->x ;

		// et mettre a jour les routes
		localSearchList[scenario]->routeU->updateRouteData();
	}

	local->rechercheTerminee = false ; 
	return 1 ;
	
}

int Individu::mutation8_indiv() {
	Params* paramsTemp = paramsList[0];
	LocalSearch* local = localSearchList[0];
	if  (local->routeU->idx == local->routeV->idx || local->routeU->depot->idx != local->routeV->depot->idx) { return 0 ; }
	double cost = 0.0;

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		LocalSearch* localTemp = localSearchList[scenario];
		double chargeResteU = localTemp->routeU->charge - localTemp->noeudU->chargeAvant ;
		double chargeResteV = localTemp->routeV->charge - localTemp->noeudV->chargeAvant ;
		double tempsU = localTemp->noeudU->est;
		double tempsV = localTemp->noeudV->est;
		double tempsResteU = localTemp->routeU->temps - tempsU - paramsTemp->timeCost[localTemp->noeudUCour][localTemp->xCour] ;
		double tempsResteV = localTemp->routeV->temps - tempsV - paramsTemp->timeCost[localTemp->noeudVCour][localTemp->yCour] ;

		cost += (paramsTemp->timeCost[localTemp->noeudUCour][localTemp->noeudVCour] 
		+ paramsTemp->timeCost[localTemp->xCour][localTemp->yCour]
		- paramsTemp->timeCost[localTemp->noeudUCour][localTemp->xCour] 
		- paramsTemp->timeCost[localTemp->noeudVCour][localTemp->yCour]
		+ localTemp->routeU->excedentCharge(localTemp->noeudU->chargeAvant + localTemp->noeudV->chargeAvant)*paramsList[scenario]->penalityCapa
		+ localTemp->routeV->excedentCharge(chargeResteV + chargeResteU)*paramsList[scenario]->penalityCapa
		- localTemp->routeU->excedentCharge(localTemp->routeU->charge)*paramsList[scenario]->penalityCapa
		- localTemp->routeV->excedentCharge(localTemp->routeV->charge)*paramsList[scenario]->penalityCapa) / (double) nbScenario;
	}

	if ( cost > -0.0001 ) { return 0 ; } 

	/////////////////////////// ON EFFECTUE LA MUTATION ///////////////////////////////
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		LocalSearch* localTemp = localSearchList[scenario];
		Noeud * depotU = localTemp->routeU->depot ;
		Noeud * depotV = localTemp->routeV->depot ;
		Noeud * depotUFin = localTemp->routeU->depot->pred ;
		Noeud * depotVFin = localTemp->routeV->depot->pred ;
		Noeud * depotVSuiv = depotV->suiv ;
		Noeud * depotUSuiv = depotU->suiv ;
		Noeud * depotVPred = depotVFin->pred ;

		// on inverse le sens et on change le nom des routes
		Noeud * temp ;
		Noeud * xx = localTemp->x ;
		Noeud * vv = localTemp->noeudV ;

		while ( !xx->estUnDepot )
		{
			temp = xx->suiv ;
			xx->suiv = xx->pred ;
			xx->pred = temp ;
			xx->route = localTemp->routeV ;
			xx = temp ;
		}

		while ( !vv->estUnDepot )
		{
			temp = vv->pred ;
			vv->pred = vv->suiv ;
			vv->suiv = temp ;
			vv->route = localTemp->routeU ;
			vv = temp ;
		}

		// mettre a jour les noeuds
		localTemp->noeudU->suiv = localTemp->noeudV ;
		localTemp->noeudV->pred = localTemp->noeudU ;
		localTemp->x->suiv = localTemp->y ;
		localTemp->y->pred = localTemp->x ;

		// mettre � jour les extr�mit�s
		if (localTemp->x->estUnDepot)
		{
			depotUFin->suiv = depotU ;
			depotUFin->pred = depotVSuiv ;
			depotUFin->pred->suiv = depotUFin ;
			depotV->suiv = localTemp->y ;
			localTemp->y->pred = depotV ;
		}
		else if ( localTemp->noeudV->estUnDepot )
		{
			depotV->suiv = depotUFin->pred ;
			depotV->suiv->pred = depotV ;
			depotV->pred = depotVFin ;
			depotUFin->pred = localTemp->noeudU ;
			localTemp->noeudU->suiv = depotUFin ;
		}
		else
		{
			depotV->suiv = depotUFin->pred ;
			depotV->suiv->pred = depotV ;
			depotUFin->pred = depotVSuiv ;
			depotUFin->pred->suiv = depotUFin ;
		}

		// et mettre a jour les routes
		localTemp->routeU->updateRouteData();
		localTemp->routeV->updateRouteData();
	}

	local->rechercheTerminee = false ; 
	return 1 ;
	
}

int Individu::mutation9_indiv() {
	Params* paramsTemp = paramsList[0];
	LocalSearch* local = localSearchList[0];
	if  (local->routeU->idx == local->routeV->idx || local->routeU->depot->idx != local->routeV->depot->idx) { return 0 ; }

	double cost = 0.0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		LocalSearch* localTemp = localSearchList[scenario];
		double chargeResteU = localTemp->routeU->charge - localTemp->noeudU->chargeAvant ;
		double chargeResteV = localTemp->routeV->charge - localTemp->noeudV->chargeAvant ;
		double tempsU = localTemp->noeudU->est;
		double tempsV = localTemp->noeudV->est;
		double tempsResteU = localTemp->routeU->temps - tempsU - paramsTemp->timeCost[localTemp->noeudUCour][localTemp->xCour] ;
		double tempsResteV = localTemp->routeV->temps - tempsV - paramsTemp->timeCost[localTemp->noeudVCour][localTemp->yCour] ;

		cost += (paramsTemp->timeCost[localTemp->noeudUCour][localTemp->yCour] 
		+ paramsTemp->timeCost[localTemp->noeudVCour][localTemp->xCour]
		- paramsTemp->timeCost[localTemp->noeudUCour][localTemp->xCour] 
		- paramsTemp->timeCost[localTemp->noeudVCour][localTemp->yCour]
		+ (localTemp->routeU->excedentCharge(localTemp->noeudU->chargeAvant + chargeResteV)
		+ localTemp->routeV->excedentCharge(localTemp->noeudV->chargeAvant + chargeResteU)
		- localTemp->routeU->excedentCharge(localTemp->routeU->charge)
		- localTemp->routeV->excedentCharge(localTemp->routeV->charge))*paramsList[scenario]->penalityCapa) / (double) nbScenario;
	}

	if (cost > -0.0001) {return 0;} 

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		LocalSearch* localTemp = localSearchList[scenario];
		Noeud * count ;
		/////////////////////////// ON EFFECTUE LA MUTATION ///////////////////////////////
		// on parcourt les noeuds pour les associer aux bonnes routes

		Noeud * depotU = localTemp->routeU->depot ;
		Noeud * depotV = localTemp->routeV->depot ;
		Noeud * depotUFin = depotU->pred ;
		Noeud * depotVFin = depotV->pred ;
		Noeud * depotUpred = depotUFin->pred ;
		Noeud * depotUSuiv = depotU->suiv ;
		Noeud * depotVPred = depotVFin->pred ;

		count = localTemp->y ;
		while ( !count->estUnDepot )
		{
			count->route = localTemp->routeU ;
			count = count->suiv ;
		}

		count = localTemp->x ;
		while ( !count->estUnDepot )
		{
			count->route = localTemp->routeV ;
			count = count->suiv ;
		}

		// mettre a jour les noeuds
		localTemp->noeudU->suiv = localTemp->y ;
		localTemp->y->pred = localTemp->noeudU ;
		localTemp->noeudV->suiv = localTemp->x ;
		localTemp->x->pred = localTemp->noeudV ;

		// mettre � jour les extr�mit�s
		if (localTemp->x->estUnDepot)
		{
			depotUFin->pred = depotVFin->pred ;
			depotUFin->pred->suiv = depotUFin ;
			localTemp->noeudV->suiv = depotVFin ;
			depotVFin->pred = localTemp->noeudV ;
		}
		else
		{
			depotUFin->pred = depotVFin->pred ;
			depotUFin->pred->suiv = depotUFin ;
			depotVFin->pred = depotUpred ;
			depotVFin->pred->suiv = depotVFin ;
		}

		localTemp->routeU->updateRouteData();
		localTemp->routeV->updateRouteData();
	}

	local->rechercheTerminee = false ;
	return 1 ;
}
