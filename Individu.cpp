#include "Individu.h"

// constructeur d'un Individu comme simple conteneur

Individu::Individu(vector<Params*> pl) : paramsList(pl)
{
	nbScenario = pl.size();
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
	for (int i = paramsTemp->nbDepots; i < paramsTemp->nbClients + paramsTemp->nbDepots; i++) {
		double startInit = paramsTemp->cli[i].startingInventory;
		vector<double> scenariosInventory(nbScenario);
		if (paramsTemp->rng->genrand64_real1() < 0.5) {
			dailyDelivery = paramsTemp->cli[i].maxInventory - startInit;
			if (dailyDelivery > 0) {
				chromT[1].push_back(i);
			}
			for (int scenario = 0; scenario < nbScenario; scenario++) {
				chromL[1 + scenario * paramsTemp->nbDays][i] = dailyDelivery;
				scenariosInventory[scenario] = paramsTemp->cli[i].maxInventory - paramsList[scenario]->cli[i].dailyDemand[1];
			}
		} else {
			for (int scenario = 0; scenario < nbScenario; scenario++) {
				scenariosInventory[scenario] = std::max(startInit - paramsList[scenario]->cli[i].dailyDemand[1], 0.0);
			}
		}
		startInventory.push_back(scenariosInventory);
	}

	//OTHER DAYS
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		paramsTemp = paramsList[scenario];
		int startIndex = 2 + scenario * (paramsTemp->nbDays);
		for (int i = paramsTemp->nbDepots; i < paramsTemp->nbClients + paramsTemp->nbDepots; i++)
		{
			double tempInventory = startInventory[i - paramsTemp->nbDepots][scenario];
			for (int k = 2; k <= paramsTemp->nbDays; k++) {
				if (tempInventory >= paramsTemp->cli[i].dailyDemand[k] || paramsTemp->rng->genrand64_real1() < 0.5)
				{
					tempInventory = std::max<double>(0., tempInventory - paramsTemp->cli[i].dailyDemand[k]);
				}
				else
				{
					dailyDelivery = paramsTemp->cli[i].maxInventory - tempInventory; 
					tempInventory = paramsTemp->cli[i].maxInventory - paramsTemp->cli[i].dailyDemand[k];
					chromL[k + scenario * paramsTemp->nbDays][i] = dailyDelivery;
					if (dailyDelivery > 0) {
						chromT[k + scenario * (paramsTemp->nbDays - 1)].push_back(i);
					}

				}
			}
		}
	}
    // TODO
	// And shuffle the whole solution
	for (int k = 1; k <= paramsList[0]->nbDays; k++)
	{
		for (int i = 0; i <= (int)chromT[k].size() - 1; i++)
		{
			int j = i + paramsList[0]->rng->genrand64_int64() % ((int)chromT[k].size() - i); 
			// swap i,j
			int temp = chromT[k][i];
			chromT[k][i] = chromT[k][j];
			chromT[k][j] = temp;
		}
	}

	// initialisation of the other structures
	pred = vector<vector<vector<int>>>(nbScenario * (paramsList[0]->nbDays - 1) + 1 + 1, vector<vector<int>>(paramsTemp->nombreVehicules[1] + 1, vector<int>(paramsTemp->nbClients + paramsTemp->nbDepots + 1, 0)));
	potentiels = vector<vector<double>>(paramsList[0]->nombreVehicules[1] + 1, vector<double>((int)paramsList[0]->nbClients + paramsList[0]->nbDepots + 1, 1.e30));
	potentiels[0][0] = 0;

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
void Individu::generalSplit_scenario()
{
	coutSol.evaluation = 0;

	// lancement de la procedure split pour chaque jour
	// on essaye deja le split simple, si c'est sans succes , le split LF
	int consideredScenario;
	if (chromT[1].size() != 0) {
		if (!splitSimpleDay1()) {
			// TODO
			std::cout << "ouUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUps" << std::endl;
		}
	}
	Params *paramsTemp;
	for (int k = 2; k < chromT.size(); k++)
	{
		consideredScenario = (k-2) / (paramsList[0]->nbDays - 1);
		paramsTemp = paramsList[consideredScenario];
		if (chromT[k].size() != 0) {
			if (!splitSimple_scenario(k, paramsTemp, consideredScenario)) {
				splitLF_scenario(k, paramsTemp, consideredScenario);
			}
		}
	}

	// After Split
	// we call a function that fills all other data structures and computes the cost
	measureSol_scenario();

	if (coutSol.evaluation > 1.e20) {
		for (int scenario = 0; scenario < nbScenario; scenario++) {
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

bool Individu::splitSimpleDay1()
{
	// on va utiliser la ligne 1 des potentiels et structures pred
	double distance, averageCost;
	vector<double> loadScenario;
	int s0, s1, sb;
	potentiels[1][0] = 0;
	s0 = paramsList[0]->ordreVehicules[1][0].depotNumber;
	for (int i = 0; i < (int)chromT[1].size(); i++)
	{
		loadScenario = std::vector<double>(nbScenario, 0.0);
		distance = 0;
		bool isCapacityOk = true;
		for (int j = i; j < (int)chromT[1].size() && isCapacityOk; j++)
		{
			s1 = chromT[1][j];
			sb = (i == j) ? s0 : chromT[1][j - 1];
			distance += paramsList[0]->timeCost[sb][s1];
			for (int scenario = 0; scenario < nbScenario; scenario++) {

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
	for (int jj = 0; jj < paramsList[0]->nombreVehicules[1]; jj++)
	{
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
	for (int cam = 0; cam < paramsTemp->nombreVehicules[day]; cam++)
	{
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

void Individu::measureSol_scenario()
{
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

	for (int scenario = 0; scenario < nbScenario; scenario++) {
		Params *paramsTemp = paramsList[scenario];
		vector<vector<double>> I_end(paramsTemp->nbDays + 1, vector<double>(paramsTemp->nbDepots + paramsTemp->nbClients));
		for (int l = paramsTemp->nbDepots; l < paramsTemp->nbDepots + paramsTemp->nbClients; l++){
			I_end[0][l] = paramsTemp->cli[l].startingInventory;
		}

		int startIndex = scenario * (paramsTemp->nbDays);

		vector<int> dayIndexL(paramsTemp->nbDays + 1, 0);
		for (int k = 1; k <= paramsTemp->nbDays; k++){
			dayIndexL[k] = startIndex + k;
		}

		int startIndexT = scenario * (paramsTemp->nbDays - 1);

		vector<int> dayIndexT(paramsTemp->nbDays + 1, 0);
		dayIndexT[1] = 1;
		for (int k = 2; k <= paramsTemp->nbDays; k++){
			dayIndexT[k] = startIndexT + k;
		}

		for (int k = 1; k <= paramsTemp->nbDays; k++){
			int day = dayIndexL[k];
			for (int cus = paramsTemp->nbDepots; cus < paramsTemp->nbDepots + paramsTemp->nbClients; cus++){
				inventoryCost[scenario] += paramsTemp->cli[cus].inventoryCost* std::max<double>(0,I_end[k-1][cus]+chromL[day][cus]-paramsTemp->cli[cus].dailyDemand[k]);
				inventoryCost[scenario] += paramsTemp->cli[cus].stockoutCost* std::max<double>(0,paramsTemp->cli[cus].dailyDemand[k]-I_end[k-1][cus]-chromL[day][cus]);
				inventoryCost[scenario] -= chromL[day][cus] * (paramsTemp->ancienNbDays + 1 - k) * paramsTemp->inventoryCostSupplier;

				I_end[k][cus] = std::max<double>(0,I_end[k-1][cus] + chromL[day][cus] - paramsTemp->cli[cus].dailyDemand[k]);
			}
		}
		for (int kk = 1; kk <= paramsTemp->nbDays; kk++)
		{
			int dayT = dayIndexT[kk];
			int dayL = dayIndexL[kk];

			j = (int)chromT[dayT].size();
			for (int jj = 0; jj < paramsTemp->nombreVehicules[kk]; jj++)
			{
				depot = paramsTemp->ordreVehicules[kk][paramsTemp->nombreVehicules[kk] - jj - 1].depotNumber;
				distance = 0;
				load = 0;
				i = (int)pred[kk][paramsTemp->nombreVehicules[kk] - jj][j];

				if (j == i)
				{
					distance = 0;
					load = 0;
				}
				else if (j == i + 1)
				{
					distance = paramsTemp->timeCost[depot][chromT[dayT][i]] + paramsTemp->timeCost[chromT[dayT][i]][depot];
					load = chromL[dayL][chromT[dayT][i]];
				}
				else
				{
					distance = paramsTemp->timeCost[depot][chromT[dayT][i]];
					load = 0;

					// infos sommets milieu
					for (int k = i; k <= j - 2; k++)
					{
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
	for (int scenario = 0; scenario < nbScenario; scenario++) {
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

void Individu::updateLS_scenario()
{
	int depot;
	int i, j;
	bool traces = false;
	vector<Noeud *> myDepot(nbScenario);
	vector<Noeud *> myDepotFin(nbScenario);
	vector<Noeud *> myClient(nbScenario);
	vector<Route *> myRoute(nbScenario);

	for (int scenario = 0; scenario < nbScenario; scenario++)
	{
		Noeud * tempDepot = myDepot[scenario];
		Noeud * tempDepotFin = myDepotFin[scenario];
		Noeud * tempClient = myClient[scenario];
		Route * tempRoute = myRoute[scenario];

		Params *paramsTemp = paramsList[scenario];
		int startIndex = scenario * (paramsTemp->nbDays) + 1;
		int startIndexT = scenario * (paramsTemp->nbDays-1) + 2;
		vector<vector<double> > demandPerDay(paramsTemp->nbDays+1, vector<double>(paramsTemp->nbClients + paramsTemp->nbDepots, 0));

		for (int k = 1; k <= paramsTemp->nbDays; k++) {
			for (int c = 0; c < paramsTemp->nbClients + paramsTemp->nbDepots; c++) {
				demandPerDay[k][c] = chromL[startIndex + k - 1][c];
			}
		}
		localSearchList[scenario]->demandPerDay = demandPerDay;

		for (int kk = 1; kk <= paramsTemp->nbDays; kk++)
		{
			localSearchList[scenario]->ordreParcours[kk].clear();
			for (i = paramsTemp->nbDepots; i < (int)localSearchList[scenario]->clients[kk].size(); i++)
			{
				localSearchList[scenario]->clients[kk][i]->estPresent = false;
			}

			int chromIndex = (kk == 1) ? 1 : startIndexT + kk - 2;
			j = (int)chromT[chromIndex].size();

			for (int jj = 0; jj < paramsTemp->nombreVehicules[kk]; jj++)
			{
				depot = paramsTemp->ordreVehicules[kk][paramsTemp->nombreVehicules[kk] - jj - 1].depotNumber;
				i = (int)pred[kk][paramsTemp->nombreVehicules[kk] - jj][j];

				tempDepot = localSearchList[scenario]->depots[kk][paramsTemp->nombreVehicules[kk] - jj - 1];
				tempDepotFin = localSearchList[scenario]->depotsFin[kk][paramsTemp->nombreVehicules[kk] - jj - 1];
				tempRoute = localSearchList[scenario]->routes[kk][paramsTemp->nombreVehicules[kk] - jj - 1];

				tempDepot->suiv = tempDepotFin;
				tempDepot->pred = tempDepotFin;
				tempDepotFin->suiv = tempDepot;
				tempDepotFin->pred = tempDepot;

				// cas ou on a un seul sommet dans le cycle
				if (j == i + 1)
				{
					tempClient = localSearchList[scenario]->clients[kk][chromT[chromIndex][i]];
					tempClient->pred = tempDepot;
					tempClient->suiv = tempDepotFin;
					tempClient->route = tempRoute;
					tempClient->estPresent = true;
					tempDepot->suiv = tempClient;
					tempDepotFin->pred = tempClient;
					localSearchList[scenario]->ordreParcours[kk].push_back(tempClient->idx);
				}
				else if (j > i + 1)
				{
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
			for (i = 0; i < (int)localSearchList[scenario]->routes[kk].size(); i++)
				localSearchList[scenario]->routes[kk][i]->updateRouteData();
		}
	}
}

void Individu::localSearchRunSearch_scenario() {
	const int GROUP_SIZE = 8;
	vector<thread> threads;

	for (int i = 0; i < nbScenario; i += GROUP_SIZE) {
		int end = min(i + GROUP_SIZE, nbScenario);
		threads.emplace_back([this, i, end]() {
			for (int scenario = i; scenario < end; scenario++) {
				localSearchList[scenario]->runSearchSameDay(false);
			}
		});
	}
	
	for (auto& t : threads) {
		t.join();
	}
	muterDifferentScenarioDP();
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
  	// while(!rechercheTerminee && nbMoves < 100){
    // 	rechercheTerminee = true;
	for (int client : randomClients) {
		nbMoves += mutationDP(client, rechercheTerminee);
	}
//   }
}

int Individu::mutationDP(int client, bool &rechercheTerminee) {
	Noeud *noeudTravail;
	double currentCost;
	// First, make sure all insertion costs are computed
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		for (int k = 1; k <= paramsList[0]->ancienNbDays; k++){
			noeudTravail = localSearchList[scenario]->clients[k][client]; //node* day k client
			localSearchList[scenario]->computeCoutInsertion(noeudTravail); // detour,place (dominated) for each route
		}
	}
	// Compute the current lot sizing solution cost (from the model point of view)
	//before optimizatio  currentCost = evaluateCurrentCost(client);
	for (int scenario = 0; scenario < nbScenario; scenario++) {
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
	
	// unique_ptr<LotSizingSolver> lotsizingSolver(
	// 		make_unique<LotSizingSolver>(paramsList, insertions, client));

		
	// lotsizingSolver->solveStockoutBackward();
		
	// vector<vector<double>> quantities = vector<vector<double>>(nbScenario, vector<double>(paramsList[0]->nbDays));
	// vector<vector<int>> breakpoints = vector<vector<int>>(nbScenario, vector<int>(paramsList[0]->nbDays));
	// vector<double> objectiveScenarios = lotsizingSolver->objective;
	// double objective = 0.0;
	// for (int scenario = 0; scenario < nbScenario; scenario++){
	// 	objective += objectiveScenarios[scenario];
	// }
	// objective /= nbScenario;
	// quantities = lotsizingSolver->quantities;
	
	// if(lt(currentCost,objective-0.01)) {
	// 	// std::cout << "VERIF: " << currentCost << " " << objective << std::endl;
	// 	return 0;
	// }

	// /* APPLYING THE MOVEMENT */
	// // Later on we will verify whether it's an improving move or not to trigger a
	// // good termination.

	// // First, removing all occurences of the node.
	// for (int scenario = 0; scenario < nbScenario; scenario++) {
	// 	for (int k = 1; k <= paramsList[scenario]->ancienNbDays; k++)
	// 	{
	// 		noeudTravail = localSearchList[scenario]->clients[k][client];
	// 		if (noeudTravail->estPresent){
	// 			localSearchList[scenario]->removeNoeud(noeudTravail);
	// 		}
	// 		localSearchList[scenario]->demandPerDay[k][client] = 0.;

	// 	}
	// 	// Then looking at the solution of the model and inserting in the good place
	// 	for (int k = 1; k <= paramsList[scenario]->ancienNbDays; k++)
	// 	{
	// 		if (quantities[scenario][k - 1] > 0.0001 || (lotsizingSolver->breakpoints[scenario][k - 1]&&gt(0,lotsizingSolver->breakpoints[scenario][k - 1]->detour) )) // don't forget that in the model the index      // goes from 0 to t-1
	// 		{
	// 		localSearchList[scenario]->demandPerDay[k][client] = round(quantities[scenario][k - 1]);
			
	// 		localSearchList[scenario]->clients[k][client]->placeInsertion = lotsizingSolver->breakpoints[scenario][k - 1]->place;
		
	// 		localSearchList[scenario]->addNoeud(localSearchList[scenario]->clients[k][client]);
	// 		}
	// 	}

	// 	double realCost = localSearchList[scenario]->evaluateCurrentCost_stockout(client);
	// 	if (fabs(realCost- objectiveScenarios[scenario])>0.01) {
	// 		std::cout << "The solution doesn't give the expected cost for scenario " << scenario << std::endl;
	// 		std::cout << "Cost: " << realCost << "; Expected cost: " << objectiveScenarios[scenario] << std::endl;
	// 		for (int scenario1 = 0; scenario1 < nbScenario; scenario1++) {
	// 			std::cout << round(quantities[scenario][0]) << std::endl;
	// 			std::cout << quantities[scenario][0] << std::endl;
	// 		}
	// 		throw string("Cost error");
	// 		return 0;
	// 	}
	// }

	// if (currentCost-objective >=0.01 )// An improving move has been found, the search is not finished.
	// {
	// 	rechercheTerminee = false;
	// 	return 1;
	// }
	// else
		return 0;
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

// mise a jour du chromT suite aux modification de localSearch
void Individu::updateIndiv_scenario()
{
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		Params *paramsTemp = paramsList[scenario];
		int startIndex = scenario * (paramsTemp->nbDays) + 1;
		int endIndex = startIndex + (paramsTemp->nbDays);
		for (int i = startIndex; i < endIndex; i++) {

			if (localSearchList[scenario] == nullptr) {
			} else if (i - startIndex + 1 >= localSearchList[scenario]->demandPerDay.size()) {
			} else {
				if (localSearchList[scenario]->demandPerDay[i - startIndex + 1].size() == 0) {
					chromL[i] = vector<double>(paramsTemp->nbClients + paramsTemp->nbDepots, 0);
				} else {
					chromL[i] = localSearchList[scenario]->demandPerDay[i - startIndex + 1];
				}
			}
		}
		vector<Route *> ordreRoutesAngle;
		Route *temp;
		Noeud *node;
		int startIndexT = scenario * (paramsTemp->nbDays - 1) + 2;
		for (int kk = 1; kk <= paramsTemp->nbDays; kk++)
		{
			int chromIndex = kk == 1 ? 1 : startIndexT + kk - 2;
			ordreRoutesAngle = localSearchList[scenario]->routes[kk];
			for (int r = 0; r < (int)ordreRoutesAngle.size(); r++)
				ordreRoutesAngle[r]->updateCentroidCoord();
			chromT[chromIndex].clear();
			for (int r = 0; r < (int)ordreRoutesAngle.size(); r++)
			{
				node = ordreRoutesAngle[r]->depot->suiv;
				while (!node->estUnDepot)
				{
					chromT[chromIndex].push_back(node->idx);
					node = node->suiv;
				}
			}
		}
	}
	generalSplit_scenario();
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
