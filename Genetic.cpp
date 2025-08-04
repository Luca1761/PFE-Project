#include "Genetic.h"
#include <algorithm> 
#include <unistd.h>

// inline double GetMemoryUsage(int pid, int flag=true) {
// 	std::ifstream status("/proc/" + std::to_string(pid) + "/status");
// 	std::string line;
// 	long long vm_size_kb = 0, vm_rss_kb = 0;

// 	while (std::getline(status, line)) {
// 		if (line.find("VmSize:") != std::string::npos) {
// 			vm_size_kb = std::stoll(line.substr(line.find(":") + 1));
// 		} else if (line.find("VmRSS:") != std::string::npos) {
// 			vm_rss_kb = std::stoll(line.substr(line.find(":") + 1));
// 		}
// 	}

// 	double vm_size_mb = (double) vm_size_kb / 1024.0;
// 	double vm_rss_mb = vm_rss_kb / 1024.0;

// 	if (flag) {
// 		std::cout << "Virtual Memory Size: " << vm_size_mb << " MB\n";
// 		std::cout << "Resident Set Size (RSS): " << vm_rss_mb << " MB\n";
// 	}

// 	return vm_size_mb;
// }

void Genetic::evolve(unsigned int maxIter, unsigned int maxIterNonProd, unsigned int maxTime) {
	unsigned int place;
	double rangRelatif;
	unsigned int oldValue = 1;
	nbIter = 0;
	nbIterNonProd = 1;
	while (nbIter < maxIter && nbIterNonProd < maxIterNonProd && (float) (clock() - population->totalTime) / CLOCKS_PER_SEC <=  (float)maxTime ) {
		// on demande deux individus a la population
		population->recopieIndividu(rejeton, population->getIndividuBinT(rangRelatif));
		population->recopieIndividu(rejeton2, population->getIndividuBinT(rangRelatif));
		// on choisit le crossover en fonction du probleme
		if (nbScenario == 1)
			crossPOX2();
		else
			crossPOX_scenario();

		muter_scenario();
		
		// REPAIR IF NEEDED
		if (!rejeton->estValide) {
			place = population->addIndividu(rejeton);
			reparer_scenario();
		}
		if (rejeton->estValide) {
			place = population->addIndividu(rejeton);
		}
		
		nbIterNonProd++;
		if (place == 0 && rejeton->estValide) {
			nbIterNonProd = 1;
			if (params->traces) std::cout << "NEW BEST FEASIBLE ";
			if (params->traces) std::cout << rejeton->coutSol.evaluation;	
			if (params->traces) std::cout << " Cost : " << rejeton->coutSol.fitness
				 << " capacity Violation : " << rejeton->coutSol.capacityViol;
			if (params->traces) std::cout << endl;
			if (params->traces) std::cout << endl;
		}
		if (nbIterNonProd / 500 == oldValue) {
			if (params->traces) std::cout << "Diversification" << endl;
			population->diversify();
			oldValue = (nbIterNonProd / 500) + 1;
		}
	
		// MANAGEMENT OF PARAMETERS
		if (nbIter % 200 == 0) gererPenalites_scenario();
	
		// TRACES
		if (nbIter % 50 == 0 && params->traces) population->afficheEtat(nbIter);
		
		nbIter++;
	}
	
	// fin de l'algorithme , diverses informations affich�es
	population->timeBest = population->timeBest - population->totalTime;
	population->totalTime = clock() - population->totalTime;
	if (params->traces) std::cout << "time passes : " << (float) (population->totalTime)/ CLOCKS_PER_SEC << endl;
	if (params->traces) std::cout << "time to find best solution : " << (float) (population->timeBest)/ CLOCKS_PER_SEC << endl;
	if (params->traces) std::cout << "number of iterations : " << nbIter << endl;
}

void Genetic::muter_scenario() {
	rejeton->generalSplit_scenario();
	rejeton->updateLS_scenario();
	rejeton->localSearchRunSearch_scenario();
	rejeton->updateIndiv_scenario();
	population->updateNbValides(rejeton);
}

// eventuellement effectue une reparation de la solution
void Genetic::reparer_scenario() {
	vector<double> savePenalities(nbScenario, 0.0);

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		savePenalities[scenario] = params->penalityCapa[scenario];
	}

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (rejeton->coutSol_scenario.capacityViol[scenario] > 0.0001)
			params->penalityCapa[scenario] *= 10;
	}
	if (params->rng->genrand64_real1() < params->pRep) {
		rejeton->generalSplit_scenario();
		rejeton->updateLS_scenario();
		rejeton->localSearchRunSearch_scenario();
		rejeton->updateIndiv_scenario();
		if (!rejeton->estValide) {
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				if (rejeton->coutSol_scenario.capacityViol[scenario] > 0.0001)
					params->penalityCapa[scenario] *= 500;
			}
			rejeton->generalSplit_scenario();
			rejeton->updateLS_scenario();
			rejeton->localSearchRunSearch_scenario();
			rejeton->updateIndiv_scenario();
		}
	}

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		params->penalityCapa[scenario] = savePenalities[scenario];
	}
}

// gestion des penalites
//TO CHECK
void Genetic::gererPenalites_scenario() {
	double fractionCharge = population->fractionValidesCharge();
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (fractionCharge < params->minValides && params->penalityCapa[scenario] < 1000)
			params->penalityCapa[scenario] *= 1.25;
		else if (fractionCharge > params->maxValides && params->penalityCapa[scenario] > 0.01)
			params->penalityCapa[scenario] *= 0.8;
	}
	population->validatePen();
}

Genetic::Genetic(Params* _params, clock_t _ticks, Population* _population) : params(_params), ticks(_ticks), population(_population)
{
	nbScenario = params->nbScenarios;
	rejeton = new Individu(params);
	rejeton2 = new Individu(params);
}

int Genetic::crossPOX_scenario() {
	unsigned int debut, fin;
    unsigned int j1, j2, l;
	unsigned int begin, end;
	unsigned int dayT, dayL, realDay;	
	// Keeping track of the chromT of the parent
	vector<vector<unsigned int>> chromTParent1 = rejeton->chromT;
	vector<vector<double>> chromLParent1 = rejeton->chromL;
	vector<vector<unsigned int>> chromTParent2 = rejeton2->chromT;
	vector<vector<double>> chromLParent2 = rejeton2->chromL;

	for (unsigned int k = 1; k < rejeton->chromL.size(); k++) {
		for (unsigned int i = 0; i < rejeton->chromL[k].size(); i++) {
			rejeton->chromL[k][i] = 0.;
		}
	}

	int firstDayInheritance = params->rng->genrand64_int64() % 4;
	rejeton->chromT[1].clear();
	if (firstDayInheritance == 0) {
		for (unsigned int ii : chromTParent1[1]) {
			rejeton->chromT[1].push_back(ii);
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				rejeton->chromL[1 + scenario * params->nbDays][ii] = chromLParent1[1 + scenario * params->nbDays][ii];
			}
		}
	} else if (firstDayInheritance == 1) {
		for (unsigned int ii : chromTParent2[1]) {
			rejeton->chromT[1].push_back(ii);
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				rejeton->chromL[1 + scenario * params->nbDays][ii] = chromLParent2[1 + scenario * params->nbDays][ii];
			}
		}
	} else {
		if (!chromTParent1[1].empty()) {
			begin = (unsigned int) (params->rng->genrand64_int64() % chromTParent1[1].size());
			end = (unsigned int) (params->rng->genrand64_int64() % chromTParent1[1].size());
			l = begin;
			while (l != end) {
				unsigned int ii = chromTParent1[1][l];
				rejeton->chromT[1].push_back(ii);
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					rejeton->chromL[1 + scenario * params->nbDays][ii] = chromLParent1[1 + scenario * params->nbDays][ii];
				}
				l = (unsigned int) ((l + 1) % chromTParent1[1].size());
			}
		}
		
		if (!chromTParent2[1].empty()) {
			for (unsigned int ii : chromTParent2[1]) {
				if (std::find(rejeton->chromT[1].begin(), rejeton->chromT[1].end(), ii) == rejeton->chromT[1].end()) {
					rejeton->chromT[1].push_back(ii);
					for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
						rejeton->chromL[1 + scenario * params->nbDays][ii] = chromLParent2[1 + scenario * params->nbDays][ii];
					}
				}
			}
		}
	}
	if (params->nbDays > 1) {
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			vector<unsigned int> garder;
			vector<int> tableauFin;
			vector<vector<bool>> hasBeenInserted = vector<vector<bool>>(params->nbDays + 1, vector<bool>(params->nbClients + params->nbDepots, false));
			unsigned int scenarioBegin = scenario * (params->nbDays - 1) + 2;
			vector<unsigned int> randomDays = vector<unsigned int>(params->nbDays - 1, 0);
			for (unsigned int i = 0; i < randomDays.size(); i++) {
				randomDays[i] = scenarioBegin + i;
			}

			std::mt19937 g((unsigned int) params->seed);
			std::shuffle(randomDays.begin(), randomDays.end(), g);

			j1 = (unsigned int) (params->rng->genrand64_int64() % randomDays.size());
			j2 = (unsigned int) (params->rng->genrand64_int64() % randomDays.size());
			if (j1 > j2) std::swap(j1, j2);
			for (unsigned int k = 0; k < randomDays.size(); k++) {
				dayT = randomDays[k];
				dayL = scenario + dayT;
				realDay = dayT - scenario * (params->nbDays - 1);
				// on recopie un segment
				if (k < j1 && !rejeton->chromT[dayT].empty()) {
					debut = (unsigned int) (params->rng->genrand64_int64() % rejeton->chromT[dayT].size());
					fin = (unsigned int) (params->rng->genrand64_int64() % rejeton->chromT[dayT].size());
					tableauFin.push_back((int) fin);
					unsigned int j = debut;
					garder.clear();
					while (j != fin) {
						unsigned int ii = rejeton->chromT[dayT][j]; // getting the index to be inherited
						double quantity = chromLParent1[dayL][ii];
						if (quantity > 0.0001) {
							garder.push_back(ii);
							rejeton->chromL[dayL][ii] = quantity;
							hasBeenInserted[realDay][ii] = true;
						}
						j = (unsigned int) ((j + 1) % rejeton->chromT[dayT].size());
					}
					rejeton->chromT[dayT].clear();
					for (unsigned int keptClient : garder) {
						rejeton->chromT[dayT].push_back(keptClient);
					}
				}
				else if (k < j2) // on recopie rien
				{
					rejeton->chromT[dayT].clear();
					tableauFin.push_back(-1);
				}
				else // on recopie tout
				{
					tableauFin.push_back(0);
					for (unsigned int ii : rejeton->chromT[dayT]) {
						garder.push_back(ii);
						rejeton->chromL[dayL][ii] = chromLParent1[dayL][ii];
						hasBeenInserted[realDay][ii] = true;
					}
				}
			}
		
			// completing with rejeton 2
			for (unsigned int k = 0; k < j2; k++) {
				dayT = randomDays[k];
				dayL = dayT + scenario;
				realDay = dayT - scenario * (params->nbDays - 1);
				for (unsigned int i = 0; i < rejeton2->chromT[dayT].size(); i++) {
					fin = (unsigned int) ((int) i + tableauFin[k] + 1);
					unsigned int ii = rejeton2->chromT[dayT][fin % rejeton2->chromT[dayT].size()];
					if (!hasBeenInserted[realDay][ii]) {
						// computing maximum possible delivery quantity
						double quantity = chromLParent2[dayL][ii];
						if (quantity > 0.0001) {
							rejeton->chromT[dayT].push_back(ii);
							rejeton->chromL[dayL][ii] = quantity;
							hasBeenInserted[realDay][ii] = true;
						}
					}
				}
			}
		}
	}
	rejeton->generalSplit_scenario();
	return 0;
}

void Genetic::crossPOX2() {
	vector<int> tableauFin;
	vector<unsigned int> garder;
	vector<unsigned int> joursPerturb;
	unsigned int debut, fin, day;
	unsigned int j1, j2;
	double quantity;

	// Keeping track of the chromL of the parent
	vector<vector<double>> chromLParent1 = rejeton->chromL;
	vector<vector<double>> chromLParent2 = rejeton2->chromL;

	// Reinitializing the chromL of the rejeton (will become the child)
	// Keeping for each day and each customer the total sum of delivered load and initial inventory
	// (when inserting a customer, need to make sure that we are not exceeding this)
	for (unsigned int k = 1; k <= params->nbDays; k++)
		for (unsigned int i = params->nbDepots; i < params->nbDepots + params->nbClients; i++)
			rejeton->chromL[k][i] = 0.;

	// Keeping a vector to remember if a delivery has already been inserted for on day k for customer i
	vector<vector<bool>> hasBeenInserted = vector<vector<bool>>(params->nbDays + 1, vector<bool>(params->nbClients + params->nbDepots, false));

	// choosing the type of inheritance for each day (nothing, all, or mixed)
	for (unsigned int k = 1; k <= params->nbDays; k++)
		joursPerturb.push_back(k);

	std::random_device rd;
	std::mt19937 g((unsigned int) params->seed);
	std::shuffle(joursPerturb.begin(), joursPerturb.end(), g);

	// Picking j1 et j2
	j1 = (unsigned int) (params->rng->genrand64_int64() % params->nbDays);
	j2 = (unsigned int) (params->rng->genrand64_int64() % params->nbDays);
	if (j1 > j2) std::swap(j1, j2);

	// Inheriting the data from rejeton1.
	// For each day, we will keep a sequence going from debut to fin
	for (unsigned int k = 0; k < params->nbDays; k++) {
		day = joursPerturb[k];
		// on recopie un segment
		if (k < j1 && !rejeton->chromT[day].empty()) {
			debut = (unsigned int)(params->rng->genrand64_int64() % rejeton->chromT[day].size());
			fin = (unsigned int)(params->rng->genrand64_int64() % rejeton->chromT[day].size());
			if (debut > fin) std::swap(debut, fin);
			tableauFin.push_back((int)fin);
			unsigned int j = debut;
			garder.clear();
			while (j != ((fin + 1) % rejeton->chromT[day].size())) {
				unsigned int ii = rejeton->chromT[day][j]; // getting the index to be inherited
				quantity = chromLParent1[day][ii];
				if (quantity > 0.0001) {
					garder.push_back(ii);
					rejeton->chromL[day][ii] = quantity;
					hasBeenInserted[day][ii] = true;
				}
				j = (unsigned int) ((j + 1) % rejeton->chromT[day].size());
			}
			rejeton->chromT[day].clear();
			for (unsigned int keptClient : garder) {
				rejeton->chromT[day].push_back(keptClient);
			}
		}
		else if (k < j2) // on recopie rien
		{
			rejeton->chromT[day].clear();
			tableauFin.push_back(-1);
		}
		else // on recopie tout
		{
			tableauFin.push_back(0);
			for (unsigned int ii : rejeton->chromT[day]) {
				garder.push_back(ii);
				rejeton->chromL[day][ii] = chromLParent1[day][ii];
				hasBeenInserted[day][ii] = true;
			}
		}
	}
	
	// completing with rejeton 2
	for (unsigned int k = 0; k < j2; k++) {
		day = joursPerturb[k];
		for (unsigned int i = 0; i < rejeton2->chromT[day].size(); i++) {
			fin = (unsigned int) ((int) i + tableauFin[k] + 1);
			unsigned int ii = rejeton2->chromT[day][fin % rejeton2->chromT[day].size()];
			if (!hasBeenInserted[day][ii]) {
				// computing maximum possible delivery quantity
				quantity = chromLParent2[day][ii];
				if (quantity > 0.0001) {
					rejeton->chromT[day].push_back(ii);
					rejeton->chromL[day][ii] = quantity;
					hasBeenInserted[day][ii] = true;
				}
			}
		}
	}
	rejeton->generalSplit_scenario();
}

Genetic::~Genetic(void)
{
	delete rejeton;
	delete rejeton2;
}
