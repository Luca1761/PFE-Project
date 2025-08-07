#include "Individual.h"

// constructeur d'un Individu comme simple conteneur

Individual::Individual(Params* _params) : params(_params)
{
	nbScenario = params->nbScenarios;
	coutSol.evaluation = 0;
	coutSol.fitness = 0;
	coutSol.capacityViol = 0;

	//not same lengths because quantities at day 1 are second stage variables while roads are first stage
	chromT = vector<vector<unsigned int>>(nbScenario * (params->nbDays - 1) + 1 + 1);
	chromL = vector<vector<double>>(nbScenario * (params->nbDays) + 1, vector<double>(params->nbClients + params->nbDepots, 0.));

	// OPTION 2 -- JUST IN TIME POLICY //
	double dailyDelivery;
	vector<vector<double>> startInventory;
	// DAY 1
	bool isFirstOption = (params->rng->genrand64_real1() < 0.5) || (params->nbDays == 1);
	for (unsigned int i = params->nbDepots; i < params->nbClients + params->nbDepots; i++) {
		double initialInventory = params->cli[i].startingInventory;
		vector<double> scenariosInventory(nbScenario);
		if (isFirstOption) {
			unsigned nb = 0;
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				nb += (initialInventory >= params->cli[i].dailyDemand[scenario][1]);
			}
			if (nb < nbScenario && params->rng->genrand64_real1() < 0.5) {
				dailyDelivery = params->cli[i].maxInventory - initialInventory;
				if (dailyDelivery > 0) chromT[1].push_back(i);
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					chromL[1 + scenario * params->nbDays][i] = dailyDelivery;
					scenariosInventory[scenario] = params->cli[i].maxInventory - params->cli[i].dailyDemand[scenario][1];
					if (scenariosInventory[scenario] < 0 || scenariosInventory[scenario] > params->cli[i].maxInventory) {
						std::cout << "WRONG INVENTORY" << std::endl;
						throw std::string("WRONG INVENTORY");
					}
				}
			} else {
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					scenariosInventory[scenario] = std::max(initialInventory - params->cli[i].dailyDemand[scenario][1], 0.0);
					if (scenariosInventory[scenario] < 0 || scenariosInventory[scenario] > params->cli[i].maxInventory) {
						std::cout << "WRONG INVENTORY" << std::endl;
						throw std::string("WRONG INVENTORY");
					}
				}
			}
		} else {
			unsigned int nb = 0;
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				nb += (initialInventory >= params->cli[i].dailyDemand[scenario][1]);
			}
			if (nb == nbScenario) {
				unsigned int nb1 = 0;
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					double nextDayClientDemand = params->cli[i].dailyDemand[scenario][2];
					nb1 += (initialInventory - params->cli[i].dailyDemand[scenario][1]) >= nextDayClientDemand;
				}
				if (nb1 != nbScenario) {
					bool shouldDeliveryForNextDay = params->rng->genrand64_real1() <= 0.7;
					if (shouldDeliveryForNextDay) {
						chromT[1].push_back(i);
						for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
							double nextDayClientDemand = params->cli[i].dailyDemand[scenario][2];
							chromL[1 + scenario * params->nbDays][i] = std::max(nextDayClientDemand - (initialInventory - params->cli[i].dailyDemand[scenario][1]), 1.0);
						}
					}
				}
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					scenariosInventory[scenario] = initialInventory - params->cli[i].dailyDemand[scenario][1] + chromL[1 + scenario * params->nbDays][i];
					if (scenariosInventory[scenario] < 0 || scenariosInventory[scenario] > params->cli[i].maxInventory) {
						std::cout << "WRONG INVENTORY" << std::endl;
						throw std::string("WRONG INVENTORY");
					}
				}	
			} else {
				chromT[1].push_back(i);
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					chromL[1 + scenario * params->nbDays][i] = std::max(params->cli[i].dailyDemand[scenario][1] - initialInventory, 1.0);
					scenariosInventory[scenario] = initialInventory - params->cli[i].dailyDemand[scenario][1] + chromL[1 + scenario * params->nbDays][i];
					if (scenariosInventory[scenario] < 0 || scenariosInventory[scenario] > params->cli[i].maxInventory) {
						std::cout << "WRONG INVENTORY" << std::endl;
						throw std::string("WRONG INVENTORY");
					}
				}
			}

		}
		startInventory.push_back(scenariosInventory);
	}
	//TO CHECK
	//OTHER DAYS
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		for (unsigned int i = params->nbDepots; i < params->nbClients + params->nbDepots; i++) {
			double tempInventory = startInventory[i - params->nbDepots][scenario];
			if (isFirstOption) {
				for (unsigned int day = 2; day <= params->nbDays; day++) {
					bool doNotDeliver = (tempInventory >= params->cli[i].dailyDemand[scenario][day] || params->rng->genrand64_real1() < 0.5);
					if (doNotDeliver) {
						tempInventory = std::max<double>(0., tempInventory - params->cli[i].dailyDemand[scenario][day]);
					} else {
						dailyDelivery = params->cli[i].maxInventory - tempInventory; 
						tempInventory = params->cli[i].maxInventory - params->cli[i].dailyDemand[scenario][day];
						chromL[day + scenario * params->nbDays][i] = dailyDelivery;
						if (dailyDelivery > 0) chromT[day + scenario * (params->nbDays - 1)].push_back(i);
					}
				}
			} else {
				for (unsigned int k = 2; k <= params->nbDays; k++) {
					double currentDayClientDemand = params->cli[i].dailyDemand[scenario][k];
					double nextDayClientDemand = params->cli[i].dailyDemand[scenario][(k + 1) % params->nbDays];

					if (tempInventory >= currentDayClientDemand) {
						// enough initial inventory, no need to service
						chromL[k + scenario * params->nbDays][i] = 0.0;

						bool isInventoryEnoughForNextDay = tempInventory - currentDayClientDemand >= nextDayClientDemand;
						if (k < params->nbDays && !isInventoryEnoughForNextDay) {
							bool shouldDeliveryForNextDay = params->rng->genrand64_real1() <= 0.3;
							if (shouldDeliveryForNextDay) {
								chromL[k + scenario * params->nbDays][i] = nextDayClientDemand - (tempInventory - currentDayClientDemand);
								chromT[k + scenario * (params->nbDays - 1)].push_back(i);
							}
						}
						tempInventory = tempInventory + chromL[k + scenario * params->nbDays][i] - currentDayClientDemand;
						if (tempInventory < 0 || tempInventory > params->cli[i].maxInventory) {
							std::cout << "WRONG INVENTORY" << std::endl;
							throw std::string("WRONG INVENTORY");
						}
					} else {
						// not enough initial inventory, just in time policy for the initial solution
						dailyDelivery = currentDayClientDemand - tempInventory;
						tempInventory = 0;
						chromL[k + scenario * params->nbDays][i] = dailyDelivery;
						chromT[k + scenario * (params->nbDays - 1)].push_back(i);
					}
				}
			}
		}
	}
	// And shuffle the whole solution
	if (params->nbDays > 1) {
		for (unsigned int day = 1; day <= params->nbDays; day++) {
			for (unsigned int i = 0; i < chromT[day].size(); i++) {
				unsigned int j = i + (unsigned int) (params->rng->genrand64_int64() % (chromT[day].size() - i)); 
				// swap i and j elements
				std::swap(chromT[day][i], chromT[day][j]);
			}
		}
	}

	// initialisation of the other structures
	pred = vector<vector<vector<unsigned int>>>(nbScenario * (params->nbDays - 1) + 1 + 1, vector<vector<unsigned int>>(params->vehicleNumber[1] + 1, vector<unsigned int>(params->nbClients + params->nbDepots + 1, 0)));
	potentiels = vector<vector<double>>(params->vehicleNumber[1] + 1, vector<double>(params->nbClients + params->nbDepots + 1, 1.e30));
	potentiels[0][0] = 0;

	localSearchList = vector<LocalSearch*>(nbScenario, nullptr);
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario] = new LocalSearch(this, params, scenario);
	}
}

// destructeur
Individual::~Individual() {
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario] != NULL) {
			delete localSearchList[scenario];
		}
	}
}

// The Split is done as follows, we test if it's possible to split without setting a fixed limit on the number of vehicles
// If the resulting solution satisfies the number of vehicle, it's perfectly fine, we return it
// Otherwise we call the Split version with limited fleet (which is a bit slower).
void Individual::generalSplit_scenario() {
	// lancement de la procedure split pour chaque jour
	// on essaye deja le split simple, si c'est sans succes , le split LF
	if (chromT[1].size() > 0) {
		while (!splitSimpleDay1() && !splitLF_scenario_day1()) {
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				params->splitBounds[scenario] *= 1.1;
			}
		}
	}
	for (unsigned int k = 2; k < chromT.size(); k++) {
		if (params->nbDays == 1) throw std::string("ERROR");
		unsigned int consideredScenario = (k-2) / (params->nbDays - 1);
		if (chromT[k].size() > 0) {
			while(!splitSimple_scenario(k, consideredScenario) && !splitLF_scenario(k, consideredScenario)) {
				params->splitBounds[consideredScenario] *= 1.1;
			}
		}
	}
	measureSol_scenario();

	params->splitBounds = std::vector<double> (params->nbScenarios, 1.0);
}

// function split which does not consider a limit on the number of vehicles
// just uses the line "1" of the "potentiels" table.
int Individual::splitSimple_scenario(unsigned int k, unsigned int scenario) {
	// on va utiliser la ligne 1 des potentiels et structures pred
	double load, distance, cost;
	unsigned int s0, s1, sb;
	potentiels[1][0] = 0;
	unsigned int day = k - scenario * (params->nbDays - 1);
	unsigned int indexL = scenario * (params->nbDays) + day;
	s0 = params->vehicleOrder[day][0].depotNumber;
	for (unsigned int i = 0; i < chromT[k].size(); i++) {
		load = 0;
		distance = 0;
		for (unsigned int j = i; j < chromT[k].size() && load <= params->vehicleOrder[day][0].capacity * params->splitBounds[scenario]; j++) {
			s1 = chromT[k][j];
			load += chromL[indexL][s1];
			sb = (i == j) ? s0 : chromT[k][j - 1];
			distance += params->timeCost[sb][s1];

			// computing the penalized cost
			cost = distance + params->timeCost[s1][s0];
			if (load > params->vehicleOrder[day][0].capacity)
				cost += (load - params->vehicleOrder[day][0].capacity) * params->penalityCapa[scenario];

			if (potentiels[1][i] + cost < potentiels[1][j + 1]) // basic Bellman algorithm
			{
				potentiels[1][j + 1] = potentiels[1][i] + cost;
				pred[k][1][j + 1] = i;
			}
		}
	}

	// testing if le number of vehicles is correct
	// in addition, the table pred is updated to keep track of everything
	unsigned int l = (unsigned int) chromT[k].size();
	for (unsigned int jj = 0; jj < params->vehicleNumber[day]; jj++) {
		pred[k][params->vehicleNumber[day] - jj][l] = pred[k][1][l];
		l = pred[k][1][l];
	}

	coutSol.evaluation = -1.e30; // just for security, making sure this value is not used (as it does not contain the constants)
	initPot_scenario(k, scenario);
	return (l == 0);
}

bool Individual::splitSimpleDay1() {
	// on va utiliser la ligne 1 des potentiels et structures pred
	double distance, averageCost;
	vector<double> loadScenario;
	unsigned int s0, s1, sb;
	potentiels[1][0] = 0;
	s0 = params->vehicleOrder[1][0].depotNumber;
	for (unsigned int i = 0; i < chromT[1].size(); i++) {
		loadScenario = std::vector<double>(nbScenario, 0.0);
		distance = 0;
		bool isCapacityOk = true;
		for (unsigned int j = i; j < chromT[1].size() && isCapacityOk; j++) {
			s1 = chromT[1][j];
			sb = (i == j) ? s0 : chromT[1][j - 1];
			distance += params->timeCost[sb][s1];
			averageCost = 0.0;
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				loadScenario[scenario] += chromL[1 + scenario * (params->nbDays)][s1];
				if (loadScenario[scenario] > params->vehicleOrder[1][0].capacity * params->splitBounds[scenario]) isCapacityOk = false;
				if (loadScenario[scenario] > params->vehicleOrder[1][0].capacity) {
					averageCost += (loadScenario[scenario] - params->vehicleOrder[1][0].capacity) * params->penalityCapa[scenario];
				}
			}
			averageCost /= (double) nbScenario;
			averageCost += distance + params->timeCost[s1][s0];

			if (potentiels[1][i] + averageCost < potentiels[1][j + 1]) // basic Bellman algorithm
			{
				potentiels[1][j + 1] = potentiels[1][i] + averageCost;
				pred[1][1][j + 1] = i;
			}
		}
	}
	unsigned int l = (unsigned int) chromT[1].size();
	for (unsigned int jj = 0; jj < params->vehicleNumber[1]; jj++) {
		pred[1][params->vehicleNumber[1] - jj][l] = pred[1][1][l];
		l = pred[1][1][l];
	}

	coutSol.evaluation = -1.e30; // just for security, making sure this value is not used (as it does not contain the constants)
	initPot_scenario(1, 0);
	return (l == 0);
}

bool Individual::splitLF_scenario_day1() {
	double distance, averageCost;
	unsigned int sb, s0, s1;

	for (unsigned int cam = 0; cam < params->vehicleNumber[1]; cam++) {
		s0 = params->vehicleOrder[1][cam].depotNumber;
		for (unsigned int i = 0; i < chromT[1].size() && potentiels[cam][i] < 1.e29; i++) {
			if (potentiels[cam][i] < potentiels[cam + 1][i]) {
				potentiels[cam + 1][i] = potentiels[cam][i];
				pred[1][cam + 1][i] = i;
			}
			std::vector<double> loadScenario = std::vector<double>(nbScenario, 0.0);
			distance = 0;
			bool isCapacityOk = true;

			for (unsigned int j = i ; j < chromT[1].size() && isCapacityOk; j++) {
				s1 = chromT[1][j];
				sb = (i == j) ? s0 : chromT[1][j - 1];
				distance +=  params->timeCost[sb][s1];

				// computing the penalized cost
				averageCost = 0.0;
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					loadScenario[scenario] += chromL[1 + scenario * (params->nbDays)][s1];
					if (loadScenario[scenario] > params->vehicleOrder[1][cam].capacity * params->splitBounds[scenario]) isCapacityOk = false;
					if (loadScenario[scenario] > params->vehicleOrder[1][cam].capacity)
						averageCost += (loadScenario[scenario] - params->vehicleOrder[1][cam].capacity) * params->penalityCapa[scenario];
				}
				averageCost /= (double)nbScenario;
				averageCost += distance + params->timeCost[s1][s0];

				if (potentiels[cam][i] + averageCost < potentiels[cam + 1][j + 1]) // Basic Bellman iteration
				{
					potentiels[cam + 1][j + 1] = potentiels[cam][i] + averageCost;
					pred[1][cam + 1][j + 1] = i;
				}
			}
		}
	}
	unsigned int l = (unsigned int) chromT[1].size();
	bool isEveryonePlaced = false;
	for (unsigned int cam = 0; cam < params->vehicleNumber[1]; cam++) {
		isEveryonePlaced |= (potentiels[cam + 1][l] < 1.e29);
	}

	// on nettoie ce que l'on a déplacé
	initPot_scenario(1, 0);
	return isEveryonePlaced;

}

// fonction split pour probl�mes � flotte limit�e
bool Individual::splitLF_scenario(unsigned int k, unsigned int scenario)
{
	double load, distance, cost;
	unsigned int sb, s0, s1;

	unsigned int day = k - scenario * (params->nbDays - 1);
	// pour chaque camion
	for (unsigned int cam = 0; cam < params->vehicleNumber[day]; cam++) {
		s0 = params->vehicleOrder[day][cam].depotNumber;
		for (unsigned int i = 0; i < chromT[k].size() && potentiels[cam][i] < 1.e29; i++) {	
			if (potentiels[cam][i] < potentiels[cam + 1][i]) {
				potentiels[cam + 1][i] = potentiels[cam][i];
				pred[k][cam + 1][i] = i;
			}
			load = 0;
			distance = 0;

			for (unsigned int j = i; j < chromT[k].size() && load <= params->vehicleOrder[day][cam].capacity * params->splitBounds[scenario]; j++) {
				s1 = chromT[k][j];
				load += chromL[k][s1];
				sb = (i == j) ? s0 : chromT[k][j - 1];
				distance += params->timeCost[sb][s1];

				// computing the penalized cost
				cost = distance + params->timeCost[s1][s0];
				if (load > params->vehicleOrder[day][cam].capacity)
					cost += (load - params->vehicleOrder[day][cam].capacity) * params->penalityCapa[scenario];

				if (potentiels[cam][i] + cost < potentiels[cam + 1][j + 1]) // Basic Bellman iteration
				{
					potentiels[cam + 1][j + 1] = potentiels[cam][i] + cost;
					pred[k][cam + 1][j + 1] = i;
				}
			}
		}
	}

	unsigned int l = (unsigned int) chromT[k].size();
	bool isEveryonePlaced = false;
	for (unsigned int cam = 0; cam < params->vehicleNumber[day]; cam++) {
		isEveryonePlaced |= (potentiels[cam + 1][l] < 1.e29);
	}
	// on nettoie ce que l'on a déplacé
	initPot_scenario(k, scenario);
	return isEveryonePlaced;
}

//TO CHECK
double Individual::measureSol(std::vector<double> &delivers, unsigned int idxDay) {
	delivers.clear();
	unsigned int depot;
	unsigned int i, j;
	double distance, load;
	double inventoryCost = 0.0;
	double stockoutCost = 0.0;
	double routeCost = 0.0;
	double capaViol = 0.0;
	double fitness = 0.0;
	if (params->traces) std::cout << std::endl;
	std::vector<bool> isDelivered(params->nbDepots + params->nbClients, false);
	if (params->traces) std::cout << "Choosen tour for this day: " << std::endl;
	for (unsigned int a : chromT[1]) {
		if (params->traces) std::cout << a << " ";
		isDelivered[a] = true;
	}
	if (params->traces) std::cout << std::endl;
	if (chromT[1].empty()) {
		if (params->traces) std::cout << "NO DELIVERY" << std::endl;
	}
	vector<double> I_end(params->nbDepots + params->nbClients);
	for (unsigned int l = params->nbDepots; l < params->nbDepots + params->nbClients; l++){
		I_end[l] = params->cli[l].startingInventory;
	}
	double supplyCost = params->availableSupply[1] * params->inventoryCostSupplier;
	for (unsigned int cus = params->nbDepots; cus < params->nbDepots + params->nbClients; cus++) {
		double toDeliver = 0.0;
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			if (chromL[1 + scenario * (params->nbDays)][cus] != 0 && !isDelivered[cus]) {
				std::cout << "SHOULD NOT DELIVER, ERROR IN THE CODE" << std::endl;
				throw std::string("SHOULD NOT DELIVER, ERROR IN THE CODE");
			}
			toDeliver += chromL[1 + scenario * (params->nbDays)][cus];
		}
		toDeliver /= (double)nbScenario;
		toDeliver = round(toDeliver);
		inventoryCost += params->cli[cus].inventoryCost * std::max<double>(0, I_end[cus] + toDeliver - params->cli[cus].trueDemand[idxDay]);
		stockoutCost += params->cli[cus].stockoutCost * std::max<double>(0, params->cli[cus].trueDemand[idxDay] - I_end[cus] - toDeliver);
		if (params->endDayInventories && I_end[cus] + toDeliver - params->cli[cus].trueDemand[idxDay] > params->cli[cus].maxInventory) {
			std::cout << "INVALID INVENTORY" << std::endl;
			throw std::string("INVALID INVENTORY");
		}
		if (!params->endDayInventories && I_end[cus] + toDeliver > params->cli[cus].maxInventory) {
			std::cout << "INVALID INVENTORY" << std::endl;
			throw std::string("INVALID INVENTORY");
		}
		supplyCost -= toDeliver * params->inventoryCostSupplier;
		delivers.push_back(toDeliver);
	}

	j = (unsigned int) chromT[1].size();
	for (unsigned int jj = 0; jj < params->vehicleNumber[1]; jj++) {
		depot = params->vehicleOrder[1][params->vehicleNumber[1] - jj - 1].depotNumber;
		distance = 0;
		load = 0;
		i = pred[1][params->vehicleNumber[1] - jj][j];
		
		if (j == i) {
			distance = 0;
			load = 0;
		} else if (j == i + 1) {
			distance = params->timeCost[depot][chromT[1][i]] + params->timeCost[chromT[1][i]][depot];
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				load += chromL[1 + scenario * (params->nbDays)][chromT[1][i]];
			}
			load /= (double) nbScenario;
		} else {
			distance = params->timeCost[depot][chromT[1][i]];
			load = 0;

			// infos sommets milieu
			for (unsigned int k = i; k <= j - 2; k++) {
				distance += params->timeCost[chromT[1][k]][chromT[1][k + 1]];
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					load += chromL[1 + scenario * (params->nbDays)][chromT[1][k]];
				}
			}

			// infos sommet fin
			distance += params->timeCost[chromT[1][j - 1]][depot];
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++)
				load += chromL[1 + scenario * (params->nbDays)][chromT[1][j - 1]];
			load /= (double) nbScenario;
		}

		routeCost += distance;

		if (load > params->vehicleOrder[1][params->vehicleNumber[1] - jj - 1].capacity) {
			capaViol += load - params->vehicleOrder[1][params->vehicleNumber[1] - jj - 1].capacity;
		}
		j = i;	
	}

	fitness = routeCost + inventoryCost + stockoutCost + supplyCost;
	if (params->traces) std::cout << "Choosen deliveries for this day: " << std::endl;
	if (params->traces) for (double a : delivers) std::cout << a << " ";
	if (params->traces) std::cout << std::endl;
	if (params->traces) std::cout << std::endl;
	if (params->traces) std::cout << "Supply inventory cost: " << supplyCost << std::endl;
	if (params->traces) std::cout << "Routing cost: " << routeCost << std::endl;
	if (params->traces) std::cout << "Client inventory cost: " << inventoryCost << std::endl;
	if (params->traces) std::cout << "Client stockout cost: " << stockoutCost << std::endl;
	if (params->traces) std::cout << std::endl;
	return fitness;
}

void Individual::measureSol_scenario() {
	unsigned int depot;
	unsigned int i, j;
	double distance, load;
	vector<double> inventoryCost(nbScenario, 0);
	vector<double> routeCost(nbScenario, 0);
	vector<double> capaViol(nbScenario, 0);
	vector<double> fitness(nbScenario, 0);
	coutSol_scenario.fitness = fitness;
	coutSol_scenario.evaluation = fitness;
	coutSol_scenario.capacityViol = capaViol;

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		vector<vector<double>> I_end(params->nbDays + 1, vector<double>(params->nbDepots + params->nbClients));
		for (unsigned int l = params->nbDepots; l < params->nbDepots + params->nbClients; l++){
			I_end[0][l] = params->cli[l].startingInventory;
		}

		unsigned int startIndexL = scenario * (params->nbDays);

		vector<unsigned int> dayIndexL(params->nbDays + 1, 0);
		for (unsigned int k = 1; k <= params->nbDays; k++){
			dayIndexL[k] = startIndexL + k;
		}

		unsigned int startIndexT = scenario * (params->nbDays - 1);

		vector<unsigned int> dayIndexT(params->nbDays + 1, 0);
		dayIndexT[1] = 1;
		for (unsigned int k = 2; k <= params->nbDays; k++){
			dayIndexT[k] = startIndexT + k;
		}

		for (unsigned int k = 1; k <= params->nbDays; k++) {
			unsigned int day = dayIndexL[k];
			for (unsigned int cus = params->nbDepots; cus < params->nbDepots + params->nbClients; cus++) {
				inventoryCost[scenario] += params->cli[cus].inventoryCost * std::max<double>(0, I_end[k-1][cus]+chromL[day][cus]-params->cli[cus].dailyDemand[scenario][k]);
				inventoryCost[scenario] += params->cli[cus].stockoutCost * std::max<double>(0, params->cli[cus].dailyDemand[scenario][k]-I_end[k-1][cus]-chromL[day][cus]);
				inventoryCost[scenario] -= chromL[day][cus] * ((double) params->nbDays + 1 - k) * params->inventoryCostSupplier;

				I_end[k][cus] = std::max<double>(0,I_end[k-1][cus] + chromL[day][cus] - params->cli[cus].dailyDemand[scenario][k]);
			}
		}
		for (unsigned int kk = 1; kk <= params->nbDays; kk++) {
			unsigned int dayT = dayIndexT[kk];
			unsigned int dayL = dayIndexL[kk];

			j = (unsigned int) chromT[dayT].size();
			for (unsigned int jj = 0; jj < params->vehicleNumber[kk]; jj++) {
				depot = params->vehicleOrder[kk][params->vehicleNumber[kk] - jj - 1].depotNumber;
				distance = 0;
				load = 0;
				i = pred[kk][params->vehicleNumber[kk] - jj][j];

				if (j == i) {
					distance = 0;
					load = 0;
				} else if (j == i + 1) {
					distance = params->timeCost[depot][chromT[dayT][i]] + params->timeCost[chromT[dayT][i]][depot];
					load = chromL[dayL][chromT[dayT][i]];
				} else {
					distance = params->timeCost[depot][chromT[dayT][i]];
					load = 0;

					// infos sommets milieu
					for (unsigned int k = i; k <= j - 2; k++) {
						distance += params->timeCost[chromT[dayT][k]][chromT[dayT][k + 1]];
						load += chromL[dayL][chromT[dayT][k]];
					}

					// infos sommet fin
					distance += params->timeCost[chromT[dayT][j - 1]][depot];
					load += chromL[dayL][chromT[dayT][j - 1]];
				}

				routeCost[scenario] += distance;
				
				if (load > params->vehicleOrder[kk][params->vehicleNumber[kk] - jj - 1].capacity) {
					capaViol[scenario] += load - params->vehicleOrder[kk][params->vehicleNumber[kk] - jj - 1].capacity;
				}
				j = i;	
			}
		}

		fitness[scenario] = routeCost[scenario] + inventoryCost[scenario] + params->objectiveConstant;
	}

	estValide = true;
	coutSol.fitness = 0;
	coutSol.evaluation = 0;
	coutSol.capacityViol = 0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		coutSol_scenario.fitness[scenario] = fitness[scenario];
		coutSol.fitness += fitness[scenario];
		coutSol_scenario.evaluation[scenario] = params->penalityCapa[scenario] * capaViol[scenario] + fitness[scenario];
		coutSol.evaluation += coutSol_scenario.evaluation[scenario];
		coutSol_scenario.capacityViol[scenario] = capaViol[scenario];
		if (capaViol[scenario] > 0.0001) {
			coutSol.capacityViol += capaViol[scenario];
			estValide = false;
		}
	}
	coutSol.evaluation /= (double)nbScenario;
	coutSol.fitness /= (double)nbScenario;
	coutSol.capacityViol /= (double)nbScenario;
}

// initialisation du vecteur potentiels
void Individual::initPot_scenario(unsigned int k, unsigned int scenario)
{
	unsigned int day = k - scenario * (params->nbDays - 1);
	for (unsigned int i = 0; i < params->vehicleNumber[day] + 1; i++) {
		for (size_t j = 0; j <= chromT[k].size() + 1; j++) {
			potentiels[i][j] = 1.e30;
		}
	}
	potentiels[0][0] = 0;
	potentiels[1][0] = 0;
}

void Individual::updateLS_scenario() {
	unsigned int i, j;
	vector<Node*> myDepot(nbScenario);
	vector<Node*> myDepotFin(nbScenario);
	vector<Node*> myClient(nbScenario);
	vector<Route*> myRoute(nbScenario);

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		Node * tempDepot = myDepot[scenario];
		Node * tempDepotFin = myDepotFin[scenario];
		Node * tempClient = myClient[scenario];
		Route * tempRoute = myRoute[scenario];

		unsigned int startIndex = scenario * (params->nbDays);
		unsigned int startIndexT = scenario * (params->nbDays - 1);
		vector<vector<double>> deliveryPerDay(params->nbDays + 1, vector<double>(params->nbClients + params->nbDepots, 0.0));

		for (unsigned int day = 1; day <= params->nbDays; day++) {
			for (unsigned int client = 0; client < params->nbClients + params->nbDepots; client++) {
				deliveryPerDay[day][client] = chromL[startIndex + day][client];
			}
		}
		localSearchList[scenario]->deliveryPerDay = deliveryPerDay;

		for (unsigned int kk = 1; kk <= params->nbDays; kk++) {
			localSearchList[scenario]->ordreParcours[kk].clear();
			for (unsigned int l = params->nbDepots; l < localSearchList[scenario]->clients[kk].size(); l++) {
				localSearchList[scenario]->clients[kk][l]->isPresent = false;
			}

			unsigned int chromIndex = (kk == 1) ? 1 : startIndexT + kk;
			j = (unsigned int) chromT[chromIndex].size();

			for (unsigned int jj = 0; jj < params->vehicleNumber[kk]; jj++) {
				i = pred[kk][params->vehicleNumber[kk] - jj][j];

				tempDepot = localSearchList[scenario]->depots[kk][params->vehicleNumber[kk] - jj - 1];
				tempDepotFin = localSearchList[scenario]->depotsFin[kk][params->vehicleNumber[kk] - jj - 1];
				tempRoute = localSearchList[scenario]->routes[kk][params->vehicleNumber[kk] - jj - 1];

				tempDepot->next = tempDepotFin;
				tempDepot->prev = tempDepotFin;
				tempDepotFin->next = tempDepot;
				tempDepotFin->prev = tempDepot;

				// cas ou on a un seul sommet dans le cycle
				if (j == i + 1) {
					tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][i]];
					tempClient->prev = tempDepot;
					tempClient->next = tempDepotFin;
					tempClient->route = tempRoute;
					tempClient->isPresent = true;
					tempDepot->next = tempClient;
					tempDepotFin->prev = tempClient;
					localSearchList[scenario]->ordreParcours[kk].push_back(tempClient->idx);
				}
				else if (j > i + 1) {
					// on a au moins 2 sommets
					// infos sommet debut
					tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][i]];
					tempClient->prev = tempDepot;
					tempClient->next = localSearchList[scenario]->clients[kk][chromT[chromIndex][i + 1]];
					tempClient->route = tempRoute;
					tempClient->isPresent = true;
					tempDepot->next = tempClient;
					localSearchList[scenario]->ordreParcours[kk].push_back(tempClient->idx);

					// infos sommet fin
					tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][j - 1]];
					tempClient->prev = localSearchList[scenario]->clients[kk][chromT[chromIndex][j - 2]];
					tempClient->next = tempDepotFin;
					tempClient->route = tempRoute;
					tempClient->isPresent = true;
					tempDepotFin->prev = tempClient;
					localSearchList[scenario]->ordreParcours[kk].push_back(tempClient->idx);

					// infos sommets milieu
					for (unsigned int k = i + 1; k <= j - 2; k++) {
						tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][k]];
						tempClient->prev = localSearchList[scenario]->clients[kk][chromT[chromIndex][k - 1]];
						tempClient->next = localSearchList[scenario]->clients[kk][chromT[chromIndex][k + 1]];
						tempClient->route = tempRoute;
						tempClient->isPresent = true;
						localSearchList[scenario]->ordreParcours[kk].push_back(tempClient->idx);
					}
				}
				j = i;
			}
			// pour chaque route on met les charges partielles à jour
			for (unsigned r = 0; r < localSearchList[scenario]->routes[kk].size(); r++)
				localSearchList[scenario]->routes[kk][r]->updateRouteData();
		}
	}
}

void Individual::localSearchRunSearch_scenario() {
	const unsigned int GROUP_SIZE = 1 + nbScenario / params->nbCores;
	
	// Local search moves (mutation1-mutation9)
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->runSearchSameDay();
	}

	// mutation1-mutation9 for day 1 (average cost reduction)
	runSearchDay1();

	// Our brand new operator of local search using dynamic programming
	muterDifferentScenarioDP();
	
	// We repeat the process after dynamic programming
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->runSearchSameDay();
	}
	runSearchDay1();
}

void Individual::runSearchDay1() {
	int nbMoves = 1;
	int nbPhases = 0;
	while (nbMoves > 0 && nbPhases < 1000) {
		localSearchList[0]->updateMoves();
		nbMoves = 0;
		nbMoves += mutationSameDay1();
		nbPhases++;
	}
}

int Individual::mutationSameDay1() {
	localSearchList[0]->dayCour = 1;
	unsigned int size = (unsigned int) localSearchList[0]->ordreParcours[1].size();
	unsigned int size2;
	localSearchList[0]->rechercheTerminee = false;
	bool moveEffectue = false;
	int nbMoves = 0;
	localSearchList[0]->firstLoop = true;
	
	while (!localSearchList[0]->rechercheTerminee) {
		localSearchList[0]->rechercheTerminee = true;
		moveEffectue = false;
		for (unsigned int posU = 0; posU < size; posU++) {
			posU -= moveEffectue; // on retourne sur le dernier noeud si on a modifie
			nbMoves += moveEffectue;
			moveEffectue = false;
			
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				localSearchList[scenario]->noeudU = localSearchList[scenario]->clients[1][localSearchList[0]->ordreParcours[1][posU]];
				localSearchList[scenario]->noeudUPred = localSearchList[scenario]->noeudU->prev;
				localSearchList[scenario]->x = localSearchList[scenario]->noeudU->next;
				localSearchList[scenario]->noeudXSuiv = localSearchList[scenario]->x->next;
				localSearchList[scenario]->xSuivCour = localSearchList[scenario]->x->next->idx;
				localSearchList[scenario]->routeU = localSearchList[scenario]->noeudU->route;
				localSearchList[scenario]->noeudUCour = localSearchList[scenario]->noeudU->idx;
				localSearchList[scenario]->noeudUPredCour = localSearchList[scenario]->noeudUPred->idx;
				localSearchList[scenario]->xCour = localSearchList[scenario]->x->idx;
			}

			size2 = (unsigned int) localSearchList[0]->noeudU->moves.size();
			for (unsigned int posV = 0; posV < size2 && moveEffectue == 0; posV++) {
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) 
				localSearchList[scenario]->noeudV = localSearchList[scenario]->clients[1][localSearchList[0]->noeudU->moves[posV]];
				if (!localSearchList[0]->noeudV->route->nodeAndRouteTested[localSearchList[0]->noeudU->idx] ||
					!localSearchList[0]->noeudU->route->nodeAndRouteTested[localSearchList[0]->noeudU->idx] || localSearchList[0]->firstLoop)
					{
						for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
							localSearchList[scenario]->noeudVPred = localSearchList[scenario]->noeudV->prev;
							localSearchList[scenario]->y = localSearchList[scenario]->noeudV->next;
							localSearchList[scenario]->noeudYSuiv = localSearchList[scenario]->y->next;
							localSearchList[scenario]->ySuivCour = localSearchList[scenario]->y->next->idx;
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
							for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
								localSearchList[scenario]->routeU->reinitSingleDayMoves();
								localSearchList[scenario]->routeV->reinitSingleDayMoves();
							}
						}
					}
			}
			
			// c'est un depot on tente l'insertion derriere le depot de ce jour
			// si il ya correlation
			if (localSearchList[0]->params->isCorrelated[localSearchList[0]->noeudU->idx][localSearchList[0]->depots[1][0]->idx] &&
				!moveEffectue)
				for (unsigned int route = 0; route < localSearchList[0]->depots[1].size(); route++)
				{
					for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
						localSearchList[scenario]->noeudV = localSearchList[scenario]->depots[1][route];
					}
					if (!localSearchList[0]->noeudV->route->nodeAndRouteTested[localSearchList[0]->noeudU->idx] ||
						!localSearchList[0]->noeudU->route->nodeAndRouteTested[localSearchList[0]->noeudU->idx] || localSearchList[0]->firstLoop)
					{
						for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
							localSearchList[scenario]->noeudVPred = localSearchList[scenario]->noeudV->prev;
							localSearchList[scenario]->y = localSearchList[scenario]->noeudV->next;
							localSearchList[scenario]->noeudYSuiv = localSearchList[scenario]->y->next;
							localSearchList[scenario]->ySuivCour = localSearchList[scenario]->y->next->idx;
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

						if (!localSearchList[0]->noeudV->next->isADepot)
						{
						if (!moveEffectue)
							moveEffectue = mutation8_indiv();
						if (!moveEffectue)
							moveEffectue = mutation9_indiv();
						}

						if (moveEffectue) {
							for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
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

void Individual::muterDifferentScenarioDP() {
	vector<unsigned int> randomClients;

	for (unsigned int client = params->nbDepots; client < params->nbClients + params->nbDepots; client++) {
		randomClients.push_back(client);
	}

	std::mt19937 g((unsigned int) params->seed);
	shuffle(randomClients.begin(), randomClients.end(), g);

	bool rechercheTerminee = false;
   	int nbMoves = 0;
	int nbPhases = 0;
	while (!rechercheTerminee) {
		rechercheTerminee = true;		
		for (unsigned int client : randomClients) {
			nbMoves += mutationDP(client, rechercheTerminee);
		}
		nbPhases++;
	}
}

int Individual::mutationDP(unsigned int client, bool &rechercheTerminee) {
	Node *noeudTravail;
	double currentCost = 0.0;
	// First, make sure all insertion costs are computed
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		for (unsigned int k = 1; k <= params->nbDays; k++){
			noeudTravail = localSearchList[scenario]->clients[k][client]; //node* day k client
			localSearchList[scenario]->computeCoutInsertion(noeudTravail); // detour,place (dominated) for each route
		}
	}
	// Compute the current lot sizing solution cost (from the model point of view)
	//before optimization
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		currentCost += localSearchList[scenario]->evaluateCurrentCost_stockout(client);
	}
	currentCost /= (double)nbScenario;
	/* Generate the structures of the subproblem */
	vector<vector<vector<Insertion>>> insertions = vector<vector<vector<Insertion>>>(nbScenario, vector<vector<Insertion>>(params->nbDays));
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		for (unsigned int k = 1; k <= params->nbDays; k++) {
			insertions[scenario][k - 1] = localSearchList[scenario]->clients[k][client]->allInsertions;
		}
	}
	
	unique_ptr<LotSizingSolver> lotsizingSolver(
			make_unique<LotSizingSolver>(params, insertions, client));

		
	lotsizingSolver->solveStockoutBackward();
		
	vector<vector<double>> quantities = vector<vector<double>>(nbScenario, vector<double>(params->nbDays));
	vector<vector<int>> breakpoints = vector<vector<int>>(nbScenario, vector<int>(params->nbDays));
	vector<double> objectiveScenarios = lotsizingSolver->objective;
	double objective = 0.0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++){
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
		for (unsigned int k = 1; k <= params->nbDays; k++) {
			noeudTravail = localSearchList[scenario]->clients[k][client];
			if (noeudTravail->isPresent){
				localSearchList[scenario]->removeNoeud(noeudTravail);
			}
			localSearchList[scenario]->deliveryPerDay[k][client] = 0.;

		}
		// Then looking at the solution of the model and inserting in the good place
		for (unsigned int k = 1; k <= params->nbDays; k++) {
			if (quantities[scenario][k - 1] > 0.0001 || (lotsizingSolver->breakpoints[scenario][k - 1]&&gt(0, lotsizingSolver->breakpoints[scenario][k - 1]->detour) )) // don't forget that in the model the index      // goes from 0 to t-1
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
			std::cout << "scenario " << scenario << std::endl;
			for (unsigned int k = 1; k <= params->nbDays; k++) {
				std::cout << quantities[scenario][k-1] << " ";
			}
			std::cout << std::endl;
			realCost = localSearchList[scenario]->evaluateCurrentCost_stockout(client, true);
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
void Individual::updateIndiv_scenario() {
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		unsigned int startIndex = scenario * (params->nbDays) + 1;
		unsigned int endIndex = startIndex + (params->nbDays);
		for (unsigned int i = startIndex; i < endIndex; i++) {  
			if (localSearchList[scenario] != nullptr && i - startIndex + 1 < localSearchList[scenario]->deliveryPerDay.size()) {
				chromL[i] = localSearchList[scenario]->deliveryPerDay[i - startIndex + 1];
			}
		}
		Node *node;
		unsigned int startIndexT = scenario * (params->nbDays - 1);
		for (unsigned int day = 1; day <= params->nbDays; day++) {
			unsigned int chromIndex = (day == 1) ? 1 : startIndexT + day;
			chromT[chromIndex].clear();
			for (Route* temp : localSearchList[scenario]->routes[day]) {
				node = temp->depot->next;
				while (!node->isADepot) {
					chromT[chromIndex].push_back(node->idx);
					node = node->next;
				}
			}
		}
	}
	generalSplit_scenario();
}

// distance generale
double Individual::distance(Individual *indiv2) {
	// TO CHECK
	double note = 0;
	bool isIdentical;
	vector<unsigned int> dayIndexL(params->nbDays + 1, 0);

	// Inventory Routing
	// distance based on number of customers which have different service days
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		int noteScenario = 0;
		unsigned int startIndex = scenario * (params->nbDays);
		for (unsigned int k = 1; k <= params->nbDays; k++)
			dayIndexL[k] = startIndex + k;
		for (unsigned int client = params->nbDepots; client < params->nbClients + params->nbDepots; client++) {
			isIdentical = true;
			for (unsigned int k : dayIndexL)
				if ((chromL[k][client] < 0.0001 && indiv2->chromL[k][client] > 0.0001) || (indiv2->chromL[k][client] < 0.0001 && chromL[k][client] > 0.0001))
					isIdentical = false;
			if (isIdentical == false)
				noteScenario++;
		}
		note += ((double)noteScenario / (double)(params->nbClients));
	}

	return ((double)note / (double)(nbScenario));
}

// ajoute un element proche dans les structures de proximite
void Individual::addProche(Individual *indiv) {
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
void Individual::removeProche(Individual *indiv) {
	list<proxData>::iterator last = plusProches.end();
	for (list<proxData>::iterator it = plusProches.begin(); it != last;) {
		if (it->indiv == indiv) {
			it = plusProches.erase(it);
			return;
		}
		else ++it;
	}
}

// distance moyenne avec les n individus les plus proches
double Individual::distPlusProche(int n) {
	double result = 0;
	double compte = 0;
	list<proxData>::iterator it = plusProches.begin();

	for (int i = 0; i < n && it != plusProches.end(); i++)
	{
		result += it->dist;
		compte += 1.0;
		++it;
	}
	if (compte == 0) return 10000000;
	return result / compte;
}

int Individual::mutation1_indiv() {
	double costSuppU = params->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->xCour] 
	- params->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->noeudUCour]  
	- params->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->xCour];

	double costSuppV = params->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->noeudUCour] 
	+ params->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->yCour] 
	- params->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->yCour];

	// dans le cas ou l'on est dans la meme route , le cout n'est pas calcule correctement en realite
	// tout ce qu'on sait c'est que si il est negatif c'est qu'il est bien reellement negatif
	// pas d'incidence pour l'instant mais attention
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV) {
			costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour])*params->penalityCapa[scenario]
			- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double)nbScenario ;

			costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour])*params->penalityCapa[scenario]
			- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double)nbScenario ;
		}
	}

	if (costSuppU + costSuppV > -0.0001) return 0 ;
	if (localSearchList[0]->noeudUCour == localSearchList[0]->yCour) return 0;

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->noeudU,localSearchList[scenario]->noeudV);
	}

	localSearchList[0]->rechercheTerminee = false;
	return 1 ;
}

int Individual::mutation2_indiv() {
	double costSuppU = params->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->xSuivCour] 
	- params->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->noeudUCour] 
	- params->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->xCour] 
	- params->timeCost[localSearchList[0]->xCour][localSearchList[0]->xSuivCour];

	double costSuppV = params->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->noeudUCour] 
	+ params->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->xCour] 
	+ params->timeCost[localSearchList[0]->xCour][localSearchList[0]->yCour] 
	- params->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->yCour];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV) {
			costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*params->penalityCapa[scenario]
			- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double) nbScenario;
			
			costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*params->penalityCapa[scenario]
			- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( localSearchList[0]->noeudU == localSearchList[0]->y || localSearchList[0]->noeudV == localSearchList[0]->x || localSearchList[0]->x->isADepot ) { return 0 ;}

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->noeudU,localSearchList[scenario]->noeudV);
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->x,localSearchList[scenario]->noeudU);
	}

	localSearchList[0]->rechercheTerminee = false ; 
	return 1 ;
}

int Individual::mutation3_indiv() {
	double costSuppU = params->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->xSuivCour] 
	- params->timeCost[localSearchList[0]->noeudUPredCour][localSearchList[0]->noeudUCour] 
	- params->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->xCour] 
	- params->timeCost[localSearchList[0]->xCour][localSearchList[0]->xSuivCour];

	double costSuppV = params->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->xCour] 
	+ params->timeCost[localSearchList[0]->xCour][localSearchList[0]->noeudUCour] 
	+ params->timeCost[localSearchList[0]->noeudUCour][localSearchList[0]->yCour] 
	- params->timeCost[localSearchList[0]->noeudVCour][localSearchList[0]->yCour];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		
		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) return 0;
	if (localSearchList[0]->noeudU == localSearchList[0]->y ||  localSearchList[0]->x == localSearchList[0]->noeudV || localSearchList[0]->x->isADepot ) return 0;

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->x,localSearchList[scenario]->noeudV);
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->noeudU,localSearchList[scenario]->x);
	}

	localSearchList[0]->rechercheTerminee = false ; 
	return 1 ;
}

int Individual::mutation4_indiv() {
	LocalSearch* local = localSearchList[0];
	double costSuppU = params->timeCost[local->noeudUPredCour][local->noeudVCour] 
	+ params->timeCost[local->noeudVCour][local->xCour]
	- params->timeCost[local->noeudUPredCour][local->noeudUCour] 
	- params->timeCost[local->noeudUCour][local->xCour];

	double costSuppV = params->timeCost[local->noeudVPredCour][local->noeudUCour] 
	+ params->timeCost[local->noeudUCour][local->yCour]
	- params->timeCost[local->noeudVPredCour][local->noeudVCour] 
	- params->timeCost[local->noeudVCour][local->yCour];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		
		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
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

int Individual::mutation5_indiv() {
	LocalSearch* local = localSearchList[0];
	double costSuppU = params->timeCost[local->noeudUPredCour][local->noeudVCour] 
	+ params->timeCost[local->noeudVCour][local->xSuivCour]
	- params->timeCost[local->noeudUPredCour][local->noeudUCour] 
	- params->timeCost[local->noeudUCour][local->xCour] 
	- params->timeCost[local->xCour][local->xSuivCour];

	double costSuppV = params->timeCost[local->noeudVPredCour][local->noeudUCour] 
	+ params->timeCost[local->xCour][local->yCour]
	+ params->timeCost[local->noeudUCour][local->xCour]
	- params->timeCost[local->noeudVPred->idx][local->noeudVCour] 
	- params->timeCost[local->noeudVCour][local->yCour];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double) nbScenario;

		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( local->noeudU == local->noeudVPred || local->x == local->noeudVPred || local->noeudU == local->y || local->x->isADepot ) { return 0 ;}

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->swapNoeud(localSearchList[scenario]->noeudU, localSearchList[scenario]->noeudV) ;
		localSearchList[scenario]->insertNoeud(localSearchList[scenario]->x, localSearchList[scenario]->noeudU);
	}
	local->rechercheTerminee = false ; 
	return 1 ;
}

int Individual::mutation6_indiv() {
	LocalSearch* local = localSearchList[0];
	double costSuppU = params->timeCost[local->noeudUPredCour][local->noeudVCour]  
	+ params->timeCost[local->noeudVCour][local->yCour]
	+ params->timeCost[local->yCour][local->xSuivCour]
	- params->timeCost[local->noeudUPredCour][local->noeudUCour] 
	- params->timeCost[local->noeudUCour][local->xCour] 
	- params->timeCost[local->xCour][local->xSuivCour];

	double costSuppV = params->timeCost[local->noeudVPredCour][local->noeudUCour] 
	+ params->timeCost[local->noeudUCour][local->xCour]
	+ params->timeCost[local->xCour][local->ySuivCour]
	- params->timeCost[local->noeudVPredCour][local->noeudVCour] 
	- params->timeCost[local->noeudVCour][local->yCour]
	- params->timeCost[local->yCour][local->ySuivCour];
	
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->yCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		
		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudUCour] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->xCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->noeudVCour] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->yCour])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return 0 ;}
	if ( local->x->isADepot || local->y->isADepot || local->y == local->noeudUPred || local->noeudU == local->y || local->x == local->noeudV || local->noeudV == local->noeudXSuiv ) { return 0 ;}

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->swapNoeud(localSearchList[scenario]->noeudU, localSearchList[scenario]->noeudV) ;
		localSearchList[scenario]->swapNoeud(localSearchList[scenario]->x,localSearchList[scenario]->y) ;
	}

	local->rechercheTerminee = false ; 
	return 1 ;
	
}

int Individual::mutation7_indiv() {
	LocalSearch* local = localSearchList[0];
	if  ((local->routeU->idx != local->routeV->idx) || local->noeudU->next == local->noeudV || local->noeudU->place > local->noeudV->place ) {  return 0 ; }
	
	double cost = params->timeCost[local->noeudUCour][local->noeudVCour] + params->timeCost[local->xCour][local->yCour]
	- params->timeCost[local->noeudUCour][local->xCour] - params->timeCost[local->noeudVCour][local->yCour] ;
	
	if ( cost > -0.0001 ) { return 0 ;}
	
	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		Node * nodeNum = localSearchList[scenario]->noeudXSuiv ;
		Node * temp ;
		localSearchList[scenario]->x->prev = nodeNum ;
		localSearchList[scenario]->x->next = localSearchList[scenario]->y ;

		while ( nodeNum != localSearchList[scenario]->noeudV )
		{
			temp = nodeNum->next ;
			nodeNum->next = nodeNum->prev ;
			nodeNum->prev = temp ;
			nodeNum = temp ;
		}

		localSearchList[scenario]->noeudV->next = localSearchList[scenario]->noeudV->prev ;
		localSearchList[scenario]->noeudV->prev = localSearchList[scenario]->noeudU ;
		localSearchList[scenario]->noeudU->next = localSearchList[scenario]->noeudV ;
		localSearchList[scenario]->y->prev = localSearchList[scenario]->x ;

		// et mettre a jour les routes
		localSearchList[scenario]->routeU->updateRouteData();
	}

	local->rechercheTerminee = false ; 
	return 1 ;
	
}

int Individual::mutation8_indiv() {
	LocalSearch* local = localSearchList[0];
	if  (local->routeU->idx == local->routeV->idx || local->routeU->depot->idx != local->routeV->depot->idx) { return 0 ; }
	double cost = 0.0;

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		LocalSearch* localTemp = localSearchList[scenario];
		double chargeResteU = localTemp->routeU->load - localTemp->noeudU->previousLoad ;
		double chargeResteV = localTemp->routeV->load - localTemp->noeudV->previousLoad ;

		cost += (params->timeCost[localTemp->noeudUCour][localTemp->noeudVCour] 
		+ params->timeCost[localTemp->xCour][localTemp->yCour]
		- params->timeCost[localTemp->noeudUCour][localTemp->xCour] 
		- params->timeCost[localTemp->noeudVCour][localTemp->yCour]
		+ localTemp->routeU->excedentCharge(localTemp->noeudU->previousLoad + localTemp->noeudV->previousLoad)*params->penalityCapa[scenario]
		+ localTemp->routeV->excedentCharge(chargeResteV + chargeResteU)*params->penalityCapa[scenario]
		- localTemp->routeU->excedentCharge(localTemp->routeU->load)*params->penalityCapa[scenario]
		- localTemp->routeV->excedentCharge(localTemp->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
	}

	if ( cost > -0.0001 ) { return 0 ; } 

	/////////////////////////// ON EFFECTUE LA MUTATION ///////////////////////////////
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		LocalSearch* localTemp = localSearchList[scenario];
		Node * depotU = localTemp->routeU->depot ;
		Node * depotV = localTemp->routeV->depot ;
		Node * depotUFin = localTemp->routeU->depot->prev ;
		Node * depotVFin = localTemp->routeV->depot->prev ;
		Node * depotVSuiv = depotV->next ;

		// on inverse le sens et on change le nom des routes
		Node * temp ;
		Node * xx = localTemp->x ;
		Node * vv = localTemp->noeudV ;

		while ( !xx->isADepot )
		{
			temp = xx->next ;
			xx->next = xx->prev ;
			xx->prev = temp ;
			xx->route = localTemp->routeV ;
			xx = temp ;
		}

		while ( !vv->isADepot )
		{
			temp = vv->prev ;
			vv->prev = vv->next ;
			vv->next = temp ;
			vv->route = localTemp->routeU ;
			vv = temp ;
		}

		// mettre a jour les noeuds
		localTemp->noeudU->next = localTemp->noeudV ;
		localTemp->noeudV->prev = localTemp->noeudU ;
		localTemp->x->next = localTemp->y ;
		localTemp->y->prev = localTemp->x ;

		// mettre � jour les extr�mit�s
		if (localTemp->x->isADepot)
		{
			depotUFin->next = depotU ;
			depotUFin->prev = depotVSuiv ;
			depotUFin->prev->next = depotUFin ;
			depotV->next = localTemp->y ;
			localTemp->y->prev = depotV ;
		}
		else if ( localTemp->noeudV->isADepot )
		{
			depotV->next = depotUFin->prev ;
			depotV->next->prev = depotV ;
			depotV->prev = depotVFin ;
			depotUFin->prev = localTemp->noeudU ;
			localTemp->noeudU->next = depotUFin ;
		}
		else
		{
			depotV->next = depotUFin->prev ;
			depotV->next->prev = depotV ;
			depotUFin->prev = depotVSuiv ;
			depotUFin->prev->next = depotUFin ;
		}

		// et mettre a jour les routes
		localTemp->routeU->updateRouteData();
		localTemp->routeV->updateRouteData();
	}

	local->rechercheTerminee = false ; 
	return 1 ;
	
}

int Individual::mutation9_indiv() {
	LocalSearch* local = localSearchList[0];
	if  (local->routeU->idx == local->routeV->idx || local->routeU->depot->idx != local->routeV->depot->idx) { return 0 ; }

	double cost = 0.0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		LocalSearch* localTemp = localSearchList[scenario];
		double chargeResteU = localTemp->routeU->load - localTemp->noeudU->previousLoad ;
		double chargeResteV = localTemp->routeV->load - localTemp->noeudV->previousLoad ;

		cost += (params->timeCost[localTemp->noeudUCour][localTemp->yCour] 
		+ params->timeCost[localTemp->noeudVCour][localTemp->xCour]
		- params->timeCost[localTemp->noeudUCour][localTemp->xCour] 
		- params->timeCost[localTemp->noeudVCour][localTemp->yCour]
		+ (localTemp->routeU->excedentCharge(localTemp->noeudU->previousLoad + chargeResteV)
		+ localTemp->routeV->excedentCharge(localTemp->noeudV->previousLoad + chargeResteU)
		- localTemp->routeU->excedentCharge(localTemp->routeU->load)
		- localTemp->routeV->excedentCharge(localTemp->routeV->load))*params->penalityCapa[scenario]) / (double) nbScenario;
	}

	if (cost > -0.0001) {return 0;} 

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		LocalSearch* localTemp = localSearchList[scenario];
		Node * count ;
		/////////////////////////// ON EFFECTUE LA MUTATION ///////////////////////////////
		// on parcourt les noeuds pour les associer aux bonnes routes

		Node * depotU = localTemp->routeU->depot ;
		Node * depotV = localTemp->routeV->depot ;
		Node * depotUFin = depotU->prev ;
		Node * depotVFin = depotV->prev ;
		Node * depotUpred = depotUFin->prev ;

		count = localTemp->y ;
		while (!count->isADepot)
		{
			count->route = localTemp->routeU ;
			count = count->next ;
		}

		count = localTemp->x ;
		while ( !count->isADepot )
		{
			count->route = localTemp->routeV ;
			count = count->next ;
		}

		// mettre a jour les noeuds
		localTemp->noeudU->next = localTemp->y ;
		localTemp->y->prev = localTemp->noeudU ;
		localTemp->noeudV->next = localTemp->x ;
		localTemp->x->prev = localTemp->noeudV ;

		// mettre � jour les extr�mit�s
		if (localTemp->x->isADepot)
		{
			depotUFin->prev = depotVFin->prev ;
			depotUFin->prev->next = depotUFin ;
			localTemp->noeudV->next = depotVFin ;
			depotVFin->prev = localTemp->noeudV ;
		}
		else
		{
			depotUFin->prev = depotVFin->prev ;
			depotUFin->prev->next = depotUFin ;
			depotVFin->prev = depotUpred ;
			depotVFin->prev->next = depotVFin ;
		}

		localTemp->routeU->updateRouteData();
		localTemp->routeV->updateRouteData();
	}

	local->rechercheTerminee = false ;
	return 1 ;
}
