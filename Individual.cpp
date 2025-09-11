#include "Individual.h"

Individual::Individual(Params* _params) : params(_params), nbScenario(params->nbScenarios)
{
	coutSol.evaluation = 0;
	coutSol.fitness = 0;
	coutSol.capacityViol = 0;

	//not same lengths because quantities at day 1 are second stage variables while roads are first stage
	chromT = vector<vector<unsigned int>>(nbScenario * (params->nbDays - 1) + 1 + 1);
	chromL = vector<vector<double>>(nbScenario * (params->nbDays) + 1, vector<double>(params->nbClients + params->nbDepots, 0.));

	double dailyDelivery;
	vector<vector<double>> startInventory(params->nbClients + params->nbDepots, vector<double>(nbScenario));
	// DAY 1

	// First policy is: if we have enough to satisfy demand, we don't deliver. Otherwise, we deliver what's missing to fill the inventory to max stockage
	// Second policy is: if we have enough to satisfy demand, we check if we also have enough to satisfy next day demand
	// if Yes, we don't deliver, otherwise we deliver just what's missing. If we don't have enough for current day, we also deliver what's missing.
	bool firstPolicy = (params->rng->genrand64_real1() < 0.5) || (params->nbDays == 1);

	for (unsigned int i = params->nbDepots; i < params->nbClients + params->nbDepots; i++) {
		double initialInventory = params->cli[i].startingInventory;
		bool areEveryInventoryEnough = true;
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			if (initialInventory < params->cli[i].dailyDemand[scenario][1]) {
				areEveryInventoryEnough = false;
				break;
			}
		}
		if (firstPolicy) {
			bool doWeDeliver = (!areEveryInventoryEnough && params->rng->genrand64_real1() < 0.5); 
			if (doWeDeliver) { // if it's not always enough or with a proba of 0.5, we deliver
				dailyDelivery = params->cli[i].maxInventory - initialInventory;
				if (dailyDelivery > 0) chromT[1].push_back(i);
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					chromL[1 + scenario * params->nbDays][i] = dailyDelivery;
					startInventory[i][scenario] = std::max(0.0, params->cli[i].maxInventory - params->cli[i].dailyDemand[scenario][1]);
				}
			} else {
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					startInventory[i][scenario] = std::max(initialInventory - params->cli[i].dailyDemand[scenario][1], 0.0);
					if (startInventory[i][scenario] > params->cli[i].maxInventory) {
						std::cout << "WRONG INVENTORY" << std::endl;
						throw std::string("WRONG INVENTORY");
					}
				}
			}
		} else {
			if (areEveryInventoryEnough) { // if enough inventory to satisfy each demand 
				unsigned int count1 = 0; // count the number of scenarios in which current inventory is enough to satisfy the first two days
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					double nextDayClientDemand = params->cli[i].dailyDemand[scenario][2];
					count1 += (initialInventory - params->cli[i].dailyDemand[scenario][1]) >= nextDayClientDemand;
				}
				if (count1 < nbScenario) { // if some scenarios don't have enough inventory
					bool shouldDeliverForNextDay = params->rng->genrand64_real1() <= 0.7;
					if (shouldDeliverForNextDay) {
						chromT[1].push_back(i);
						for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
							double nextDayClientDemand = params->cli[i].dailyDemand[scenario][2];
							chromL[1 + scenario * params->nbDays][i] = std::max(std::min(nextDayClientDemand + params->cli[i].dailyDemand[scenario][1], params->cli[i].maxInventory) - initialInventory, 1.0);
						}
					}
				}
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					startInventory[i][scenario] = initialInventory - params->cli[i].dailyDemand[scenario][1] + chromL[1 + scenario * params->nbDays][i];
					if (initialInventory +  chromL[1 + scenario * params->nbDays][i] > params->cli[i].maxInventory && !params->endDayInventories) {
						std::cout << "WRONG INVENTORY" << std::endl;
						throw std::string("WRONG INVENTORY");
					}
				}	
			} else { // if not enough
				chromT[1].push_back(i);
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					chromL[1 + scenario * params->nbDays][i] = std::max(std::min(params->cli[i].dailyDemand[scenario][1], params->cli[i].maxInventory) - initialInventory, 1.0);
					startInventory[i][scenario] = std::max(0.0, initialInventory - params->cli[i].dailyDemand[scenario][1] + chromL[1 + scenario * params->nbDays][i]);
					if (initialInventory +  chromL[1 + scenario * params->nbDays][i] > params->cli[i].maxInventory && !params->endDayInventories) {
						std::cout << "WRONG INVENTORY" << std::endl;
						throw std::string("WRONG INVENTORY");
					}
				}
			}
		}
	}

	//OTHER DAYS
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		for (unsigned int i = params->nbDepots; i < params->nbClients + params->nbDepots; i++) {
			double tempInventory = startInventory[i][scenario];
			if (firstPolicy) {
				for (unsigned int day = 2; day <= params->nbDays; day++) {
					bool doNotDeliver = (tempInventory >= params->cli[i].dailyDemand[scenario][day] || params->rng->genrand64_real1() < 0.5);
					if (doNotDeliver) {
						tempInventory = std::max<double>(0.0, tempInventory - params->cli[i].dailyDemand[scenario][day]);
					} else {
						dailyDelivery = params->cli[i].maxInventory - tempInventory; 
						tempInventory = std::max(0.0, params->cli[i].maxInventory - params->cli[i].dailyDemand[scenario][day]);
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

						bool isInventoryEnoughForNextDay = tempInventory - currentDayClientDemand >= nextDayClientDemand;
						if (k < params->nbDays && !isInventoryEnoughForNextDay) {
							bool shouldDeliveryForNextDay = params->rng->genrand64_real1() <= 0.3;
							if (shouldDeliveryForNextDay) {
								chromL[k + scenario * params->nbDays][i] = std::min(nextDayClientDemand + currentDayClientDemand, params->cli[i].maxInventory) - tempInventory;
								if (chromL[k + scenario * params->nbDays][i] > 0)
									chromT[k + scenario * (params->nbDays - 1)].push_back(i);
							}
						}
						tempInventory = tempInventory + chromL[k + scenario * params->nbDays][i] - currentDayClientDemand;
						if (tempInventory + currentDayClientDemand > params->cli[i].maxInventory) {
							std::cout << "WRONG INVENTORY" << std::endl;
							throw std::string("WRONG INVENTORY");
						}
					} else {
						// not enough initial inventory, just in time policy for the initial solution
						dailyDelivery = std::min(currentDayClientDemand, params->cli[i].maxInventory) - tempInventory;
						tempInventory = 0;
						chromL[k + scenario * params->nbDays][i] = dailyDelivery;
						if (dailyDelivery > 0)
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

	// initialization of the other structures
	pred = vector<vector<vector<unsigned int>>>(nbScenario * (params->nbDays - 1) + 1 + 1, vector<vector<unsigned int>>(params->vehicleNumber[1] + 1, vector<unsigned int>(params->nbClients + params->nbDepots + 1, 0)));
	potentials = vector<vector<double>>(params->vehicleNumber[1] + 1, vector<double>(params->nbClients + params->nbDepots + 1, 1.e30));
	potentials[0][0] = 0;

	localSearchList = vector<LocalSearch*>(nbScenario, nullptr);
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario] = new LocalSearch(this, params, scenario);
	}
}

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
void Individual::split() {
	// launch the split procedure for every day (different procedure for first day, as the route is shared by every scenario)
	// we try the best split, otherwise we use the limited fleet split
	if (chromT[1].size() > 0) {
		while (!bestSplitFirstDay() && !splitLimitedFleetFirstDay()) {
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				params->splitBounds[scenario] *= 1.1; // we must split, even if it's not feasible and if we need to overcharge capacity (but we try without it)
			}
		}
	}

	// day >= 2
	for (unsigned int k = 2; k < chromT.size(); k++) {
		if (params->nbDays == 1) throw std::string("ERROR");
		unsigned int scenario = (k-2) / (params->nbDays - 1);
		if (chromT[k].size() > 0) {
			while(!bestSplit(k, scenario) && !splitLimitedFleet(k, scenario)) {
				params->splitBounds[scenario] *= 1.1;
			}
		}
	}

	measureSol(); // we compute the solution cost

	params->splitBounds = std::vector<double> (params->nbScenarios, 1.0); // reinitialize splitbounds (used only here)
}

// function split which does not consider a limit on the number of vehicles and find the best possible split
// just uses the line "1" of the "potentials" table.
int Individual::bestSplit(unsigned int k, unsigned int scenario) {
	double load, distance, cost;
	unsigned int s0, s1, sb;
	potentials[1][0] = 0;
	unsigned int day = k - scenario * (params->nbDays - 1);
	unsigned int indexL = scenario * (params->nbDays) + day;
	s0 = params->vehicleOrder[day][0].depotNumber;
	for (unsigned int i = 0; i < chromT[k].size(); i++) {
		load = 0;
		distance = 0;
		for (unsigned int j = i; j < chromT[k].size() && load <= params->vehicleOrder[day][0].capacity * params->splitBounds[scenario]; j++) {
			s1 = chromT[k][j];
			load += chromL[indexL][s1];
			if (load > params->vehicleOrder[day][0].capacity * params->splitBounds[scenario]) break;

			sb = (i == j) ? s0 : chromT[k][j - 1];
			distance += params->timeCost[sb][s1];

			// computing the penalized cost
			cost = distance + params->timeCost[s1][s0];
			if (load > params->vehicleOrder[day][0].capacity)
				cost += (load - params->vehicleOrder[day][0].capacity) * params->penalityCapa[scenario];

			if (potentials[1][i] + cost < potentials[1][j + 1]) // basic Bellman algorithm
			{
				potentials[1][j + 1] = potentials[1][i] + cost;
				pred[k][1][j + 1] = i;
			}
		}
	}

	// testing if the number of vehicles is correct
	// in addition, the table pred is updated to keep track of everything
	unsigned int l = (unsigned int) chromT[k].size();
	for (unsigned int jj = 0; jj < params->vehicleNumber[day]; jj++) {
		pred[k][params->vehicleNumber[day] - jj][l] = pred[k][1][l];
		l = pred[k][1][l];
	}

	coutSol.evaluation = -1.e30; // just for security, making sure this value is not used (as it does not contain the constants)
	initPotentials(k, scenario);
	return (l == 0); // if l == 0, means that every client has been placed in one of the routes
}

bool Individual::bestSplitFirstDay() {
	double distance, averageCost; // this time, it's an average cost because the first route is shared by every scenario
	vector<double> loadScenario;
	unsigned int s0, s1, sb;
	potentials[1][0] = 0;
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
			if (!isCapacityOk) break;
			averageCost /= (double) nbScenario;
			averageCost += distance + params->timeCost[s1][s0];

			if (potentials[1][i] + averageCost < potentials[1][j + 1]) // basic Bellman algorithm
			{
				potentials[1][j + 1] = potentials[1][i] + averageCost;
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
	initPotentials(1, 0);
	return (l == 0); // if l == 0, means that every client has been placed in one of the routes
}

bool Individual::splitLimitedFleetFirstDay() {
	double distance, averageCost;
	unsigned int sb, s0, s1;

	for (unsigned int cam = 0; cam < params->vehicleNumber[1]; cam++) {
		s0 = params->vehicleOrder[1][cam].depotNumber;
		for (unsigned int i = 0; i < chromT[1].size() && potentials[cam][i] < 1.e29; i++) {
			if (potentials[cam][i] < potentials[cam + 1][i]) {
				potentials[cam + 1][i] = potentials[cam][i];
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
				if (!isCapacityOk) break;
				averageCost /= (double)nbScenario;
				averageCost += distance + params->timeCost[s1][s0];

				if (potentials[cam][i] + averageCost < potentials[cam + 1][j + 1]) // Basic Bellman iteration
				{
					potentials[cam + 1][j + 1] = potentials[cam][i] + averageCost;
					pred[1][cam + 1][j + 1] = i;
				}
			}
		}
	}
	unsigned int l = (unsigned int) chromT[1].size();
	bool isEveryonePlaced = false;
	for (unsigned int cam = 0; cam < params->vehicleNumber[1]; cam++) {
		isEveryonePlaced |= (potentials[cam + 1][l] < 1.e29); // we check if the last client of the tour is visited by a truck
															  // if it is, it means everyone has been placed, since one of the loop
															  // avoid to do next customer while the current one is not visited
	}

	initPotentials(1, 0);
	return isEveryonePlaced;

}

bool Individual::splitLimitedFleet(unsigned int k, unsigned int scenario)
{
	double load, distance, cost;
	unsigned int sb, s0, s1;

	unsigned int day = k - scenario * (params->nbDays - 1);
	for (unsigned int cam = 0; cam < params->vehicleNumber[day]; cam++) {
		s0 = params->vehicleOrder[day][cam].depotNumber;
		for (unsigned int i = 0; i < chromT[k].size() && potentials[cam][i] < 1.e29; i++) {	
			if (potentials[cam][i] < potentials[cam + 1][i]) {
				potentials[cam + 1][i] = potentials[cam][i];
				pred[k][cam + 1][i] = i;
			}
			load = 0;
			distance = 0;

			for (unsigned int j = i; j < chromT[k].size() && load <= params->vehicleOrder[day][cam].capacity * params->splitBounds[scenario]; j++) {
				s1 = chromT[k][j];
				load += chromL[k][s1];
				if (load > params->vehicleOrder[day][cam].capacity * params->splitBounds[scenario]) break;
				sb = (i == j) ? s0 : chromT[k][j - 1];
				distance += params->timeCost[sb][s1];

				// computing the penalized cost
				cost = distance + params->timeCost[s1][s0];
				if (load > params->vehicleOrder[day][cam].capacity)
					cost += (load - params->vehicleOrder[day][cam].capacity) * params->penalityCapa[scenario];

				if (potentials[cam][i] + cost < potentials[cam + 1][j + 1]) // Basic Bellman iteration
				{
					potentials[cam + 1][j + 1] = potentials[cam][i] + cost;
					pred[k][cam + 1][j + 1] = i;
				}
			}
		}
	}

	unsigned int l = (unsigned int) chromT[k].size();
	bool isEveryonePlaced = false;
	for (unsigned int cam = 0; cam < params->vehicleNumber[day]; cam++) {
		isEveryonePlaced |= (potentials[cam + 1][l] < 1.e29); // same that for previous function
	}

	initPotentials(k, scenario);
	return isEveryonePlaced;
}

double Individual::measureTrueCost(std::vector<double> &delivers) {
	delivers.clear();
	unsigned int depot;
	unsigned int i, j;
	double distance, load;
	double inventoryCost = 0.0;
	double stockoutCost = 0.0;
	double routeCost = 0.0;
	double capaViol = 0.0;

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

	double supplyCost = params->availableSupply[1] * params->inventoryCostSupplier; // cost from supplier inventory

	for (unsigned int cus = params->nbDepots; cus < params->nbDepots + params->nbClients; cus++) {
		double toDeliver = 0.0;
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			if (chromL[1 + scenario * (params->nbDays)][cus] > 0 && !isDelivered[cus]) { // check that a client is not delivered if he's not in the tour
				std::cout << "SHOULD NOT DELIVER, ERROR IN THE CODE" << std::endl;
				throw std::string("SHOULD NOT DELIVER, ERROR IN THE CODE");
			}
			toDeliver += chromL[1 + scenario * (params->nbDays)][cus];
		}
		toDeliver /= (double) nbScenario; // we take the average delivered quantity on every scenario
		toDeliver = round(toDeliver);

		inventoryCost += params->cli[cus].inventoryCost * std::max<double>(0, params->cli[cus].startingInventory + toDeliver - params->cli[cus].trueDemand[params->jVal]);
		stockoutCost += params->cli[cus].stockoutCost * std::max<double>(0, params->cli[cus].trueDemand[params->jVal] - params->cli[cus].startingInventory - toDeliver);

		if (params->endDayInventories && params->cli[cus].startingInventory + toDeliver - params->cli[cus].trueDemand[params->jVal] > params->cli[cus].maxInventory) {
			std::cout << "INVALID INVENTORY" << std::endl;
			throw std::string("INVALID INVENTORY");
		}
		if (!params->endDayInventories && params->cli[cus].startingInventory + toDeliver > params->cli[cus].maxInventory) {
			std::cout << "INVALID INVENTORY" << std::endl;
			throw std::string("INVALID INVENTORY");
		}

		supplyCost -= toDeliver * params->inventoryCostSupplier; // we remove from supplier cost what's used to deliver
		delivers.push_back(toDeliver);
	}

	// we compute time costs
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

			for (unsigned int k = i; k <= j - 2; k++) {
				distance += params->timeCost[chromT[1][k]][chromT[1][k + 1]];
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					load += chromL[1 + scenario * (params->nbDays)][chromT[1][k]];
				}
			}

			distance += params->timeCost[chromT[1][j - 1]][depot];
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++)
				load += chromL[1 + scenario * (params->nbDays)][chromT[1][j - 1]];
			load /= (double) nbScenario;
		}

		routeCost += distance;

		if (load > params->vehicleOrder[1][params->vehicleNumber[1] - jj - 1].capacity) {
			capaViol += (load - params->vehicleOrder[1][params->vehicleNumber[1] - jj - 1].capacity);
		}
		j = i;	
	}

	if (params->traces) {
		std::cout << "Choosen deliveries for this day: " << std::endl;
		for (double a : delivers) 
			std::cout << a << " ";
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "Supply inventory cost: " << supplyCost << std::endl;
		std::cout << "Routing cost: " << routeCost << std::endl;
		std::cout << "Client inventory cost: " << inventoryCost << std::endl;
		std::cout << "Client stockout cost: " << stockoutCost << std::endl;
		std::cout << std::endl;
	}

	return routeCost + inventoryCost + stockoutCost + supplyCost;
}

void Individual::measureSol() {
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
	isFeasible = true;

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		vector<vector<double>> I_end(params->nbDays + 1, vector<double>(params->nbDepots + params->nbClients));
		for (unsigned int l = params->nbDepots; l < params->nbDepots + params->nbClients; l++){
			I_end[0][l] = params->cli[l].startingInventory;
		}

		vector<unsigned int> dayIndexL(params->nbDays + 1);
		for (unsigned int k = 1; k <= params->nbDays; k++){
			dayIndexL[k] = scenario * (params->nbDays) + k;
		}

		vector<unsigned int> dayIndexT(params->nbDays + 1, 1);
		for (unsigned int k = 2; k <= params->nbDays; k++){
			dayIndexT[k] = scenario * (params->nbDays - 1) + k;
		}

		for (unsigned int k = 1; k <= params->nbDays; k++) {
			unsigned int day = dayIndexL[k];
			for (unsigned int cus = params->nbDepots; cus < params->nbDepots + params->nbClients; cus++) {
				inventoryCost[scenario] += params->cli[cus].inventoryCost * std::max<double>(0, I_end[k-1][cus]+chromL[day][cus]-params->cli[cus].dailyDemand[scenario][k]);
				inventoryCost[scenario] += params->cli[cus].stockoutCost * std::max<double>(0, params->cli[cus].dailyDemand[scenario][k]-I_end[k-1][cus]-chromL[day][cus]);
				inventoryCost[scenario] -= chromL[day][cus] * ((double) params->nbDays + 1 - k) * params->inventoryCostSupplier;
                
				double maxDeliverable = (params->endDayInventories) ? params->cli[cus].dailyDemand[scenario][k] + params->cli[cus].maxInventory : params->cli[cus].maxInventory;
				double stockOvercharge = 1000000 * std::max<double>(I_end[k-1][cus] + chromL[day][cus] - maxDeliverable, 0.0);
				if (stockOvercharge > 0) isFeasible = false;
				inventoryCost[scenario] += stockOvercharge;
				
				I_end[k][cus] = std::max<double>(0, I_end[k-1][cus] + chromL[day][cus] - params->cli[cus].dailyDemand[scenario][k]);
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

					for (unsigned int k = i; k <= j - 2; k++) {
						distance += params->timeCost[chromT[dayT][k]][chromT[dayT][k + 1]];
						load += chromL[dayL][chromT[dayT][k]];
					}

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
			isFeasible = false;
		}
	}
	coutSol.evaluation /= (double)nbScenario;
	coutSol.fitness /= (double)nbScenario;
	coutSol.capacityViol /= (double)nbScenario;
}

void Individual::initPotentials(unsigned int k, unsigned int scenario) {
	unsigned int day = k - scenario * (params->nbDays - 1);
	for (unsigned int i = 0; i < params->vehicleNumber[day] + 1; i++) {
		for (size_t j = 0; j <= chromT[k].size() + 1; j++) {
			potentials[i][j] = 1.e30;
		}
	}
	potentials[0][0] = 0;
	potentials[1][0] = 0;
}

void Individual::updateLocalSearch() {
	unsigned int i, j;
	vector<Node*> myDepot(nbScenario);
	vector<Node*> myDepotFin(nbScenario);
	vector<Node*> myClient(nbScenario);
	vector<Route*> myRoute(nbScenario);

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		Node* tempDepot = myDepot[scenario];
		Node* tempDepotFin = myDepotFin[scenario];
		Node* tempClient = myClient[scenario];
		Route* tempRoute = myRoute[scenario];

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
			localSearchList[scenario]->clientOrder[kk].clear();
			for (unsigned int l = params->nbDepots; l < localSearchList[scenario]->clients[kk].size(); l++) {
				localSearchList[scenario]->clients[kk][l]->isPresent = false;
			}

			unsigned int chromIndex = (kk == 1) ? 1 : startIndexT + kk;
			j = (unsigned int) chromT[chromIndex].size();

			for (unsigned int jj = 0; jj < params->vehicleNumber[kk]; jj++) {
				i = pred[kk][params->vehicleNumber[kk] - jj][j];

				tempDepot = localSearchList[scenario]->depots[kk][params->vehicleNumber[kk] - jj - 1];
				tempDepotFin = localSearchList[scenario]->endDepots[kk][params->vehicleNumber[kk] - jj - 1];
				tempRoute = localSearchList[scenario]->routes[kk][params->vehicleNumber[kk] - jj - 1];

				tempDepot->next = tempDepotFin;
				tempDepot->prev = tempDepotFin;
				tempDepotFin->next = tempDepot;
				tempDepotFin->prev = tempDepot;

				// first case: only one client in the tour
				if (j == i + 1) {
					tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][i]];
					tempClient->prev = tempDepot;
					tempClient->next = tempDepotFin;
					tempClient->route = tempRoute;
					tempClient->isPresent = true;
					tempDepot->next = tempClient;
					tempDepotFin->prev = tempClient;
					localSearchList[scenario]->clientOrder[kk].push_back(tempClient->idx);
				}
				else if (j > i + 1) {
					// second case: at least two clients
					// first client
					tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][i]];
					tempClient->prev = tempDepot;
					tempClient->next = localSearchList[scenario]->clients[kk][chromT[chromIndex][i + 1]];
					tempClient->route = tempRoute;
					tempClient->isPresent = true;
					tempDepot->next = tempClient;
					localSearchList[scenario]->clientOrder[kk].push_back(tempClient->idx);

					// last client
					tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][j - 1]];
					tempClient->prev = localSearchList[scenario]->clients[kk][chromT[chromIndex][j - 2]];
					tempClient->next = tempDepotFin;
					tempClient->route = tempRoute;
					tempClient->isPresent = true;
					tempDepotFin->prev = tempClient;
					localSearchList[scenario]->clientOrder[kk].push_back(tempClient->idx);

					// intermediate clients
					for (unsigned int k = i + 1; k <= j - 2; k++) {
						tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][k]];
						tempClient->prev = localSearchList[scenario]->clients[kk][chromT[chromIndex][k - 1]];
						tempClient->next = localSearchList[scenario]->clients[kk][chromT[chromIndex][k + 1]];
						tempClient->route = tempRoute;
						tempClient->isPresent = true;
						localSearchList[scenario]->clientOrder[kk].push_back(tempClient->idx);
					}
				}
				j = i;
			}
			// for each route, we update the partial charges
			for (unsigned r = 0; r < localSearchList[scenario]->routes[kk].size(); r++)
				localSearchList[scenario]->routes[kk][r]->updateRouteData();
		}
	}
}

void Individual::runLocalSearch() {
	// const unsigned int GROUP_SIZE = 1 + nbScenario / params->nbCores; // to use when multithreading
	
	// Local search moves (mutation1-mutation9) (no multithreading yet because the random generator is shared by scenarios...)
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->runSearchSameDay();
	}

	// mutation1-mutation9 for day 1 (average cost reduction)
	runSearchFirstDay();

	// Our brand new operator of local search using dynamic programming (multithreading)
	backwardDPOperator();
	
	// We repeat the process after dynamic programming
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->runSearchSameDay();
	}
	runSearchFirstDay();
}

void Individual::runSearchFirstDay() {
	bool stop = false;
	int nbPhases = 0;
	while (!stop && nbPhases < 1000) { // limit on number of phases because big instances can be difficult to solve
		localSearchList[0]->updateMoves();
		stop = (mutationFirstDay() == 0);
		nbPhases++;
	}
}

int Individual::mutationFirstDay() {
	localSearchList[0]->currDay = 1;
	unsigned int size = (unsigned int) localSearchList[0]->clientOrder[1].size();
	unsigned int size2;
	localSearchList[0]->stopResearch = false;
	bool moveEffectue = false;
	int nbMoves = 0;
	localSearchList[0]->firstLoop = true;
	
	while (!localSearchList[0]->stopResearch) {
		localSearchList[0]->stopResearch = true;
		moveEffectue = false;
		for (unsigned int posU = 0; posU < size; posU++) {
			posU -= moveEffectue; // come back on last node if modified
			nbMoves += moveEffectue;
			moveEffectue = false;
			
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				localSearchList[scenario]->nodeU = localSearchList[scenario]->clients[1][localSearchList[0]->clientOrder[1][posU]];
				localSearchList[scenario]->nodeUPrev = localSearchList[scenario]->nodeU->prev;
				localSearchList[scenario]->x = localSearchList[scenario]->nodeU->next;
				localSearchList[scenario]->nodeXNext = localSearchList[scenario]->x->next;
				localSearchList[scenario]->idxXNext = localSearchList[scenario]->x->next->idx;
				localSearchList[scenario]->routeU = localSearchList[scenario]->nodeU->route;
				localSearchList[scenario]->idxNodeU = localSearchList[scenario]->nodeU->idx;
				localSearchList[scenario]->idxNodeUPrev = localSearchList[scenario]->nodeUPrev->idx;
				localSearchList[scenario]->idxX = localSearchList[scenario]->x->idx;
			}

			size2 = (unsigned int) localSearchList[0]->nodeU->moves.size();
			for (unsigned int posV = 0; posV < size2 && !moveEffectue; posV++) {
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) 
				localSearchList[scenario]->nodeV = localSearchList[scenario]->clients[1][localSearchList[0]->nodeU->moves[posV]];
				if (!localSearchList[0]->nodeV->route->nodeAndRouteTested[localSearchList[0]->nodeU->idx] ||
					!localSearchList[0]->nodeU->route->nodeAndRouteTested[localSearchList[0]->nodeU->idx] || localSearchList[0]->firstLoop)
					{
						for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
							localSearchList[scenario]->nodeVPrev = localSearchList[scenario]->nodeV->prev;
							localSearchList[scenario]->y = localSearchList[scenario]->nodeV->next;
							localSearchList[scenario]->nodeYNext = localSearchList[scenario]->y->next;
							localSearchList[scenario]->idxYNext = localSearchList[scenario]->y->next->idx;
							localSearchList[scenario]->routeV = localSearchList[scenario]->nodeV->route;
							localSearchList[scenario]->idxNodeV = localSearchList[scenario]->nodeV->idx;
							localSearchList[scenario]->idxNodeVPrev = localSearchList[scenario]->nodeVPrev->idx;
							localSearchList[scenario]->idxY = localSearchList[scenario]->y->idx;
						}
						
						if (!moveEffectue)
						moveEffectue = mutation1_indiv();
						if (!moveEffectue)
						moveEffectue = mutation2_indiv();
						if (!moveEffectue)
						moveEffectue = mutation3_indiv();
						
						// mutations 4 and 6 (switch) are symetrical
						if (localSearchList[0]->nodeU->idx <= localSearchList[0]->nodeV->idx) {
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
			
			if (localSearchList[0]->params->isCorrelated[localSearchList[0]->nodeU->idx][localSearchList[0]->depots[1][0]->idx] &&
				!moveEffectue)
				for (unsigned int route = 0; route < localSearchList[0]->depots[1].size(); route++)
				{
					for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
						localSearchList[scenario]->nodeV = localSearchList[scenario]->depots[1][route];
					}
					if (!localSearchList[0]->nodeV->route->nodeAndRouteTested[localSearchList[0]->nodeU->idx] ||
						!localSearchList[0]->nodeU->route->nodeAndRouteTested[localSearchList[0]->nodeU->idx] || localSearchList[0]->firstLoop)
					{
						for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
							localSearchList[scenario]->nodeVPrev = localSearchList[scenario]->nodeV->prev;
							localSearchList[scenario]->y = localSearchList[scenario]->nodeV->next;
							localSearchList[scenario]->nodeYNext = localSearchList[scenario]->y->next;
							localSearchList[scenario]->idxYNext = localSearchList[scenario]->y->next->idx;
							localSearchList[scenario]->routeV = localSearchList[scenario]->nodeV->route;
							localSearchList[scenario]->idxNodeV = localSearchList[scenario]->nodeV->idx;
							localSearchList[scenario]->idxNodeVPrev = localSearchList[scenario]->nodeVPrev->idx;
							localSearchList[scenario]->idxY = localSearchList[scenario]->y->idx;
						}

						if (!moveEffectue)
						moveEffectue = mutation1_indiv();
						if (!moveEffectue)
						moveEffectue = mutation2_indiv();
						if (!moveEffectue)
						moveEffectue = mutation3_indiv();

						if (!localSearchList[0]->nodeV->next->isADepot)
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
		}
		localSearchList[0]->firstLoop = false;
	}
	return nbMoves;
}

void Individual::backwardDPOperator() {
	vector<unsigned int> randomClients;

	for (unsigned int client = params->nbDepots; client < params->nbClients + params->nbDepots; client++) {
		randomClients.push_back(client);
	}

	std::mt19937 g((unsigned int) params->seed);
	shuffle(randomClients.begin(), randomClients.end(), g);

	bool stopResearch = false;
	while (!stopResearch) { // while there's an improving move, we continue the algorithm
		stopResearch = true;		
		for (unsigned int client : randomClients) { // we remove client from every scenario and every day and try to reinsert it
			backwardDPSingleClient(client, stopResearch);
		}
	}
}

void Individual::backwardDPSingleClient(unsigned int client, bool &stopResearch) {
	Node *virtualNode;
	// First, make sure all insertion costs are computed
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		for (unsigned int k = 1; k <= params->nbDays; k++){
			virtualNode = localSearchList[scenario]->clients[k][client];
			localSearchList[scenario]->computeInsertionCost(virtualNode); // detour, place for each route
		}
	}

	// Compute the current lot sizing solution cost (from the model point of view)
	// before optimization
	double currentCost = 0.0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		currentCost += localSearchList[scenario]->evaluateCurrentClientCost(client);
	}
	currentCost /= (double)nbScenario; // we take the mean on all scenarios

	/* Generate the structures of the subproblem */
	vector<vector<vector<Insertion>>> insertions = vector<vector<vector<Insertion>>>(nbScenario, vector<vector<Insertion>>(params->nbDays));
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		for (unsigned int k = 1; k <= params->nbDays; k++) {
			insertions[scenario][k - 1] = localSearchList[scenario]->clients[k][client]->allInsertions;
		}
	}
	
	unique_ptr<LotSizingSolver> lotsizingSolver(make_unique<LotSizingSolver>(params, insertions, client));
	
	lotsizingSolver->solveStockoutBackward(); // solve by DP + backtrack
		
	// convert DP results
	double objective = 0.0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++){
		objective += lotsizingSolver->objective[scenario];
	}
	objective /= (double)nbScenario;
	
	if(lt(currentCost, objective - 0.01)) { // if the new cost if not better than the old one, we stop
		return;
	}
	
	/* APPLYING THE MOVEMENT */
	// Later on we will verify whether it's an improving move or not to trigger a
	// good termination.
	
	vector<vector<double>> quantities = lotsizingSolver->quantities;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		// First, removing all occurences of the node.
		for (unsigned int k = 1; k <= params->nbDays; k++) {
			virtualNode = localSearchList[scenario]->clients[k][client];
			if (virtualNode->isPresent){
				localSearchList[scenario]->removeNode(virtualNode);
			}
			localSearchList[scenario]->deliveryPerDay[k][client] = 0.;

		}
		// Then looking at the solution of the model and inserting in the good place
		for (unsigned int k = 1; k <= params->nbDays; k++) {
			if (quantities[scenario][k - 1] > 0.0001) { // if we deliver
				localSearchList[scenario]->deliveryPerDay[k][client] = round(quantities[scenario][k - 1]);
			
				localSearchList[scenario]->clients[k][client]->placeInsertion = lotsizingSolver->breakpoints[scenario][k - 1]->place;
		
				localSearchList[scenario]->addNode(localSearchList[scenario]->clients[k][client]);
			}
		}

		double realCost = localSearchList[scenario]->evaluateCurrentClientCost(client);
		double expectedCost = lotsizingSolver->objective[scenario];
		if (fabs(realCost - expectedCost) > 0.001) { // we check if backtracked solution has the good cost
			std::cout << "The solution doesn't give the expected cost for scenario " << scenario << std::endl;
			std::cout << "Cost: " << realCost << "; Expected cost: " << expectedCost << std::endl;
			std::cout << "Quantities: " << std::endl;
			for (unsigned int k = 1; k <= params->nbDays; k++) {
				std::cout << quantities[scenario][k-1] << " ";
			}
			std::cout << std::endl;
			realCost = localSearchList[scenario]->evaluateCurrentClientCost(client);
			throw string("Cost error");
		}
	}
	stopResearch &= ge(objective + 0.01, currentCost); // An improving move has been found, the search is not finished.
	return;
}

void Individual::updateIndividual() {
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		unsigned int startIndexL = scenario * (params->nbDays) + 1;
		for (unsigned int i = startIndexL; i < startIndexL + params->nbDays; i++) {
			chromL[i] = localSearchList[scenario]->deliveryPerDay[i - startIndexL + 1];
		}
		Node *tempNode;
		unsigned int startIndexT = scenario * (params->nbDays - 1);
		for (unsigned int day = 1; day <= params->nbDays; day++) {
			unsigned int chromIndex = (day == 1) ? 1 : startIndexT + day;
			chromT[chromIndex].clear();
			for (Route* temp : localSearchList[scenario]->routes[day]) {
				tempNode = temp->depot->next;
				while (!tempNode->isADepot) {
					chromT[chromIndex].push_back(tempNode->idx);
					tempNode = tempNode->next;
				}
			}
		}
	}
	split();
}

double Individual::distance(Individual *indiv2) {
	double note = 0.0;
	vector<unsigned int> dayIndexL(params->nbDays + 1, 0);

	// Inventory Routing
	// distance based on number of customers which have different service days
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		int noteScenario = 0;
		for (unsigned int k = 1; k <= params->nbDays; k++)
			dayIndexL[k] = scenario * (params->nbDays) + k;
		for (unsigned int client = params->nbDepots; client < params->nbClients + params->nbDepots; client++) {
			for (unsigned int k : dayIndexL) {
				if ((chromL[k][client] < 0.0001 && indiv2->chromL[k][client] > 0.0001) || (indiv2->chromL[k][client] < 0.0001 && chromL[k][client] > 0.0001)) {
					noteScenario++;
					break;
				}
			}
		}
		note += (double)noteScenario / (double)(params->nbClients);
	}

	return note / (double)(nbScenario);
}

void Individual::addNearest(Individual *indiv) {
	proxData data;
	data.indiv = indiv;
	data.dist = distance(indiv);

	if (nearestIndiv.empty()) {
		nearestIndiv.push_back(data);
	} else {
		for (list<proxData>::iterator it = nearestIndiv.begin(); it != nearestIndiv.end(); it++) {
			if (it->dist >= data.dist) {
				nearestIndiv.insert(it, data);
				return;
			}
		}
	}
}

void Individual::removeNearest(Individual *indiv) {
	for (list<proxData>::iterator it = nearestIndiv.begin(); it != nearestIndiv.end(); it++) {
		if (it->indiv == indiv) {
			it = nearestIndiv.erase(it);
			return;
		}
	}
}

double Individual::distNearest(int n) {
	if (nearestIndiv.empty()) return 10000000;

	list<proxData>::iterator it = nearestIndiv.begin();
	double result = 0.0;
	int compte = 0;
	for (int i = 0; i < n && it != nearestIndiv.end(); i++) {
		result += it->dist;
		compte ++;
		++it;
	}
	return result / (double)compte;
}

bool Individual::mutation1_indiv() {
	double costSuppU = params->timeCost[localSearchList[0]->idxNodeUPrev][localSearchList[0]->idxX] 
	- params->timeCost[localSearchList[0]->idxNodeUPrev][localSearchList[0]->idxNodeU]  
	- params->timeCost[localSearchList[0]->idxNodeU][localSearchList[0]->idxX];

	double costSuppV = params->timeCost[localSearchList[0]->idxNodeV][localSearchList[0]->idxNodeU] 
	+ params->timeCost[localSearchList[0]->idxNodeU][localSearchList[0]->idxY] 
	- params->timeCost[localSearchList[0]->idxNodeV][localSearchList[0]->idxY];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV) {
			costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU])*params->penalityCapa[scenario]
			- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double)nbScenario ;

			costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU])*params->penalityCapa[scenario]
			- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double)nbScenario ;
		}
	}

	if (costSuppU + costSuppV > -0.0001) return false;
	if (localSearchList[0]->idxNodeU == localSearchList[0]->idxY) return false;

	// update nodes
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->insertNode(localSearchList[scenario]->nodeU,localSearchList[scenario]->nodeV);
	}

	localSearchList[0]->stopResearch = false;
	return true;
}

bool Individual::mutation2_indiv() {
	double costSuppU = params->timeCost[localSearchList[0]->idxNodeUPrev][localSearchList[0]->idxXNext] 
	- params->timeCost[localSearchList[0]->idxNodeUPrev][localSearchList[0]->idxNodeU] 
	- params->timeCost[localSearchList[0]->idxNodeU][localSearchList[0]->idxX] 
	- params->timeCost[localSearchList[0]->idxX][localSearchList[0]->idxXNext];

	double costSuppV = params->timeCost[localSearchList[0]->idxNodeV][localSearchList[0]->idxNodeU] 
	+ params->timeCost[localSearchList[0]->idxNodeU][localSearchList[0]->idxX] 
	+ params->timeCost[localSearchList[0]->idxX][localSearchList[0]->idxY] 
	- params->timeCost[localSearchList[0]->idxNodeV][localSearchList[0]->idxY];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV) {
			costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxX])*params->penalityCapa[scenario]
			- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double) nbScenario;
			
			costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxX])*params->penalityCapa[scenario]
			- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return false;}
	if ( localSearchList[0]->nodeU == localSearchList[0]->y || localSearchList[0]->nodeV == localSearchList[0]->x || localSearchList[0]->x->isADepot ) { return false;}

	// update nodes
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->insertNode(localSearchList[scenario]->nodeU,localSearchList[scenario]->nodeV);
		localSearchList[scenario]->insertNode(localSearchList[scenario]->x,localSearchList[scenario]->nodeU);
	}

	localSearchList[0]->stopResearch = false ; 
	return true;
}

bool Individual::mutation3_indiv() {
	double costSuppU = params->timeCost[localSearchList[0]->idxNodeUPrev][localSearchList[0]->idxXNext] 
	- params->timeCost[localSearchList[0]->idxNodeUPrev][localSearchList[0]->idxNodeU] 
	- params->timeCost[localSearchList[0]->idxNodeU][localSearchList[0]->idxX] 
	- params->timeCost[localSearchList[0]->idxX][localSearchList[0]->idxXNext];

	double costSuppV = params->timeCost[localSearchList[0]->idxNodeV][localSearchList[0]->idxX] 
	+ params->timeCost[localSearchList[0]->idxX][localSearchList[0]->idxNodeU] 
	+ params->timeCost[localSearchList[0]->idxNodeU][localSearchList[0]->idxY] 
	- params->timeCost[localSearchList[0]->idxNodeV][localSearchList[0]->idxY];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxX])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		
		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxX])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) return false;
	if (localSearchList[0]->nodeU == localSearchList[0]->y ||  localSearchList[0]->x == localSearchList[0]->nodeV || localSearchList[0]->x->isADepot ) return false;

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->insertNode(localSearchList[scenario]->x,localSearchList[scenario]->nodeV);
		localSearchList[scenario]->insertNode(localSearchList[scenario]->nodeU,localSearchList[scenario]->x);
	}

	localSearchList[0]->stopResearch = false ; 
	return true;
}

bool Individual::mutation4_indiv() {
	LocalSearch* local = localSearchList[0];
	double costSuppU = params->timeCost[local->idxNodeUPrev][local->idxNodeV] 
	+ params->timeCost[local->idxNodeV][local->idxX]
	- params->timeCost[local->idxNodeUPrev][local->idxNodeU] 
	- params->timeCost[local->idxNodeU][local->idxX];

	double costSuppV = params->timeCost[local->idxNodeVPrev][local->idxNodeU] 
	+ params->timeCost[local->idxNodeU][local->idxY]
	- params->timeCost[local->idxNodeVPrev][local->idxNodeV] 
	- params->timeCost[local->idxNodeV][local->idxY];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeV] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		
		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeV])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return false;}
	if ( local->idxNodeU == local->idxNodeVPrev || local->idxNodeU == local->idxY) { return false;}

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->swapNode(localSearchList[scenario]->nodeU, localSearchList[scenario]->nodeV) ;
	}

	localSearchList[0]->stopResearch = false ; 
	return true;
}

bool Individual::mutation5_indiv() {
	LocalSearch* local = localSearchList[0];
	double costSuppU = params->timeCost[local->idxNodeUPrev][local->idxNodeV] 
	+ params->timeCost[local->idxNodeV][local->idxXNext]
	- params->timeCost[local->idxNodeUPrev][local->idxNodeU] 
	- params->timeCost[local->idxNodeU][local->idxX] 
	- params->timeCost[local->idxX][local->idxXNext];

	double costSuppV = params->timeCost[local->idxNodeVPrev][local->idxNodeU] 
	+ params->timeCost[local->idxX][local->idxY]
	+ params->timeCost[local->idxNodeU][local->idxX]
	- params->timeCost[local->nodeVPrev->idx][local->idxNodeV] 
	- params->timeCost[local->idxNodeV][local->idxY];

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeV] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxX])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double) nbScenario;

		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxX] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeV])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return false;}
	if ( local->nodeU == local->nodeVPrev || local->x == local->nodeVPrev || local->nodeU == local->y || local->x->isADepot ) { return false;}

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->swapNode(localSearchList[scenario]->nodeU, localSearchList[scenario]->nodeV) ;
		localSearchList[scenario]->insertNode(localSearchList[scenario]->x, localSearchList[scenario]->nodeU);
	}
	local->stopResearch = false ; 
	return true;
}

bool Individual::mutation6_indiv() {
	LocalSearch* local = localSearchList[0];
	double costSuppU = params->timeCost[local->idxNodeUPrev][local->idxNodeV]  
	+ params->timeCost[local->idxNodeV][local->idxY]
	+ params->timeCost[local->idxY][local->idxXNext]
	- params->timeCost[local->idxNodeUPrev][local->idxNodeU] 
	- params->timeCost[local->idxNodeU][local->idxX] 
	- params->timeCost[local->idxX][local->idxXNext];

	double costSuppV = params->timeCost[local->idxNodeVPrev][local->idxNodeU] 
	+ params->timeCost[local->idxNodeU][local->idxX]
	+ params->timeCost[local->idxX][local->idxYNext]
	- params->timeCost[local->idxNodeVPrev][local->idxNodeV] 
	- params->timeCost[local->idxNodeV][local->idxY]
	- params->timeCost[local->idxY][local->idxYNext];
	
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (localSearchList[scenario]->routeU != localSearchList[scenario]->routeV)
		{
		costSuppU += (localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeV] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxY] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxX])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeU->excedentCharge(localSearchList[scenario]->routeU->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		
		costSuppV += (localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeU] + localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxX] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxNodeV] - localSearchList[scenario]->deliveryPerDay[1][localSearchList[scenario]->idxY])*params->penalityCapa[scenario]
		- localSearchList[scenario]->routeV->excedentCharge(localSearchList[scenario]->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
		}
	}

	if ( costSuppU + costSuppV > -0.0001 ) { return false;}
	if ( local->x->isADepot || local->y->isADepot || local->y == local->nodeUPrev || local->nodeU == local->y || local->x == local->nodeV || local->nodeV == local->nodeXNext ) { return false;}

	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		localSearchList[scenario]->swapNode(localSearchList[scenario]->nodeU, localSearchList[scenario]->nodeV) ;
		localSearchList[scenario]->swapNode(localSearchList[scenario]->x,localSearchList[scenario]->y) ;
	}

	local->stopResearch = false ; 
	return true;
	
}

bool Individual::mutation7_indiv() {
	LocalSearch* local = localSearchList[0];
	if  ((local->routeU->idx != local->routeV->idx) || local->nodeU->next == local->nodeV || local->nodeU->place > local->nodeV->place ) {  return false; }
	
	double cost = params->timeCost[local->idxNodeU][local->idxNodeV] + params->timeCost[local->idxX][local->idxY]
	- params->timeCost[local->idxNodeU][local->idxX] - params->timeCost[local->idxNodeV][local->idxY] ;
	
	if ( cost > -0.0001 ) { return false;}
	
	// mettre a jour les noeuds
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		Node * nodeNum = localSearchList[scenario]->nodeXNext ;
		Node * temp ;
		localSearchList[scenario]->x->prev = nodeNum ;
		localSearchList[scenario]->x->next = localSearchList[scenario]->y ;

		while ( nodeNum != localSearchList[scenario]->nodeV )
		{
			temp = nodeNum->next ;
			nodeNum->next = nodeNum->prev ;
			nodeNum->prev = temp ;
			nodeNum = temp ;
		}

		localSearchList[scenario]->nodeV->next = localSearchList[scenario]->nodeV->prev ;
		localSearchList[scenario]->nodeV->prev = localSearchList[scenario]->nodeU ;
		localSearchList[scenario]->nodeU->next = localSearchList[scenario]->nodeV ;
		localSearchList[scenario]->y->prev = localSearchList[scenario]->x ;

		// et mettre a jour les routes
		localSearchList[scenario]->routeU->updateRouteData();
	}

	local->stopResearch = false ; 
	return true;
	
}

bool Individual::mutation8_indiv() {
	LocalSearch* local = localSearchList[0];
	if  (local->routeU->idx == local->routeV->idx || local->routeU->depot->idx != local->routeV->depot->idx) { return false; }
	double cost = 0.0;

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		LocalSearch* localTemp = localSearchList[scenario];
		double chargeResteU = localTemp->routeU->load - localTemp->nodeU->previousLoad ;
		double chargeResteV = localTemp->routeV->load - localTemp->nodeV->previousLoad ;

		cost += (params->timeCost[localTemp->idxNodeU][localTemp->idxNodeV] 
		+ params->timeCost[localTemp->idxX][localTemp->idxY]
		- params->timeCost[localTemp->idxNodeU][localTemp->idxX] 
		- params->timeCost[localTemp->idxNodeV][localTemp->idxY]
		+ localTemp->routeU->excedentCharge(localTemp->nodeU->previousLoad + localTemp->nodeV->previousLoad)*params->penalityCapa[scenario]
		+ localTemp->routeV->excedentCharge(chargeResteV + chargeResteU)*params->penalityCapa[scenario]
		- localTemp->routeU->excedentCharge(localTemp->routeU->load)*params->penalityCapa[scenario]
		- localTemp->routeV->excedentCharge(localTemp->routeV->load)*params->penalityCapa[scenario]) / (double) nbScenario;
	}

	if ( cost > -0.0001 ) { return false; } 

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
		Node * vv = localTemp->nodeV ;

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
		localTemp->nodeU->next = localTemp->nodeV ;
		localTemp->nodeV->prev = localTemp->nodeU ;
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
		else if ( localTemp->nodeV->isADepot )
		{
			depotV->next = depotUFin->prev ;
			depotV->next->prev = depotV ;
			depotV->prev = depotVFin ;
			depotUFin->prev = localTemp->nodeU ;
			localTemp->nodeU->next = depotUFin ;
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

	local->stopResearch = false ; 
	return true;
	
}

bool Individual::mutation9_indiv() {
	LocalSearch* local = localSearchList[0];
	if  (local->routeU->idx == local->routeV->idx || local->routeU->depot->idx != local->routeV->depot->idx) { return false; }

	double cost = 0.0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		LocalSearch* localTemp = localSearchList[scenario];
		double chargeResteU = localTemp->routeU->load - localTemp->nodeU->previousLoad ;
		double chargeResteV = localTemp->routeV->load - localTemp->nodeV->previousLoad ;

		cost += (params->timeCost[localTemp->idxNodeU][localTemp->idxY] 
		+ params->timeCost[localTemp->idxNodeV][localTemp->idxX]
		- params->timeCost[localTemp->idxNodeU][localTemp->idxX] 
		- params->timeCost[localTemp->idxNodeV][localTemp->idxY]
		+ (localTemp->routeU->excedentCharge(localTemp->nodeU->previousLoad + chargeResteV)
		+ localTemp->routeV->excedentCharge(localTemp->nodeV->previousLoad + chargeResteU)
		- localTemp->routeU->excedentCharge(localTemp->routeU->load)
		- localTemp->routeV->excedentCharge(localTemp->routeV->load))*params->penalityCapa[scenario]) / (double) nbScenario;
	}

	if (cost > -0.0001) {return false;} 

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
		while (!count->isADepot) {
			count->route = localTemp->routeU ;
			count = count->next ;
		}

		count = localTemp->x ;
		while (!count->isADepot) {
			count->route = localTemp->routeV ;
			count = count->next ;
		}

		// mettre a jour les noeuds
		localTemp->nodeU->next = localTemp->y ;
		localTemp->y->prev = localTemp->nodeU ;
		localTemp->nodeV->next = localTemp->x ;
		localTemp->x->prev = localTemp->nodeV ;

		// mettre � jour les extr�mit�s
		if (localTemp->x->isADepot)
		{
			depotUFin->prev = depotVFin->prev ;
			depotUFin->prev->next = depotUFin ;
			localTemp->nodeV->next = depotVFin ;
			depotVFin->prev = localTemp->nodeV ;
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

	local->stopResearch = false ;
	return true;
}
