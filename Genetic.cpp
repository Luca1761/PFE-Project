#include "Genetic.h"
#include <algorithm> 
#include <unistd.h>

inline double GetMemoryUsage(int pid, int flag=true) {
	std::ifstream status("/proc/" + std::to_string(pid) + "/status");
	std::string line;
	long long vm_size_kb = 0, vm_rss_kb = 0;

	while (std::getline(status, line)) {
		if (line.find("VmSize:") != std::string::npos) {
			vm_size_kb = std::stoll(line.substr(line.find(":") + 1));
		} else if (line.find("VmRSS:") != std::string::npos) {
			vm_rss_kb = std::stoll(line.substr(line.find(":") + 1));
		}
	}

	double vm_size_mb = vm_size_kb / 1024.0;
	double vm_rss_mb = vm_rss_kb / 1024.0;

	if (flag) {
		std::cout << "Virtual Memory Size: " << vm_size_mb << " MB\n";
		std::cout << "Resident Set Size (RSS): " << vm_rss_mb << " MB\n";
	}

	return vm_size_mb;
}

void Genetic::evolve(int maxIter, int maxIterNonProd) {
	int place;
	double rangRelatif;
	nbIterNonProd = 1;
	nbIter = 0;
	int oldValue = 1;
	while (nbIter < maxIter && nbIterNonProd < maxIterNonProd && clock() <=  3600 * CLOCKS_PER_SEC) {
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
			cout << "NEW BEST FEASIBLE ";
			cout << rejeton->coutSol.evaluation;	
			cout << " Cost : " << rejeton->coutSol.fitness
				 << " capacity Violation : " << rejeton->coutSol.capacityViol;
			cout << endl;
			cout << endl;
		}
		if (nbIterNonProd / 750 == oldValue) {
			cout << "Diversification" << endl;
			population->diversify();
			oldValue = (nbIterNonProd / 750) + 1;
		}
	
		// MANAGEMENT OF PARAMETERS
		if (nbIter % 100 == 0) gererPenalites_scenario();
	
		// TRACES
		if (nbIter % 50 == 0) population->afficheEtat(nbIter);
		
		nbIter++;
	}
	
	// fin de l'algorithme , diverses informations affich�es
	cout << "time passes : " << clock() << endl;
	cout << "number of iterations : " << nbIter << endl;
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
		savePenalities[scenario] = paramsList[scenario]->penalityCapa;
	}

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (rejeton->coutSol_scenario.capacityViol[scenario] > 0.0001)
			paramsList[scenario]->penalityCapa *= 10;
	}
	if (paramsList[0]->rng->genrand64_real1() < paramsList[0]->pRep) {
		rejeton->generalSplit_scenario();
		rejeton->updateLS_scenario();
		rejeton->localSearchRunSearch_scenario();
		rejeton->updateIndiv_scenario();
		if (!rejeton->estValide) {
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				if (rejeton->coutSol_scenario.capacityViol[scenario] > 0.0001)
					paramsList[scenario]->penalityCapa *= 500;
			}
			rejeton->generalSplit_scenario();
			rejeton->updateLS_scenario();
			rejeton->localSearchRunSearch_scenario();
			rejeton->updateIndiv_scenario();
		}
	}

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		paramsList[scenario]->penalityCapa = savePenalities[scenario];
	}
}

// gestion des penalites
//TO CHECK
void Genetic::gererPenalites_scenario() {
	double fractionCharge = population->fractionValidesCharge();
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (fractionCharge < paramsList[scenario]->minValides && paramsList[scenario]->penalityCapa < 1000)
			paramsList[scenario]->penalityCapa *= 1.25;
		else if (fractionCharge > paramsList[scenario]->maxValides && paramsList[scenario]->penalityCapa > 0.01)
			paramsList[scenario]->penalityCapa *= 0.8;
	}
	population->validatePen();
}

Genetic::Genetic(vector<Params*> pl, Population *population, clock_t ticks, bool traces) : paramsList(pl), population(population), ticks(ticks), traces(traces)
{
	nbScenario = pl.size();
	rejeton = new Individu(paramsList);
	rejeton2 = new Individu(paramsList);
}

int Genetic::crossPOX_scenario() {
    vector<int> vide, tableauEtat;
    vector<double> charge;
	int debut, fin;
    int j1, j2;
	// Keeping track of the chromT of the parent
	vector<vector<int>> chromTParent1 = rejeton->chromT;
	vector<vector<double>> chromLParent1 = rejeton->chromL;
	vector<vector<int>> chromTParent2 = rejeton2->chromT;
	vector<vector<double>> chromLParent2 = rejeton2->chromL;

	for (unsigned int k = 1; k < rejeton->chromL.size(); k++) {
		for (unsigned int i = 0; i < rejeton->chromL[k].size(); i++) {
			rejeton->chromL[k][i] = 0.;
		}
	}

	int firstDayInheritance = paramsList[0]->rng->genrand64_int64() % 4;
	rejeton->chromT[1].clear();
	if (firstDayInheritance == 0) {
		for (int ii : chromTParent1[1]) {
			rejeton->chromT[1].push_back(ii);
			for (int scenario = 0; scenario < nbScenario; scenario++) {
				rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = chromLParent1[1 + scenario * paramsList[0]->nbDays][ii];
			}
		}
	} else if (firstDayInheritance == 1) {
		for (int ii : chromTParent2[1]) {
			rejeton->chromT[1].push_back(ii);
			for (int scenario = 0; scenario < nbScenario; scenario++) {
				rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = chromLParent2[1 + scenario * paramsList[0]->nbDays][ii];
			}
		}
	} else {
		if (!chromTParent1[1].empty()) {
			int begin = paramsList[0]->rng->genrand64_int64() % chromTParent1[1].size();
			int end = paramsList[0]->rng->genrand64_int64() % chromTParent1[1].size();
			int i = begin;
			while (i != end) {
				int ii = chromTParent1[1][i];
				rejeton->chromT[1].push_back(ii);
				for (int scenario = 0; scenario < nbScenario; scenario++) {
					rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = chromLParent1[1 + scenario * paramsList[0]->nbDays][ii];
				}
				i = (i + 1) % chromTParent1[1].size();
			}
		}
		
		if (!chromTParent2[1].empty()) {
			for (int ii : chromTParent2[1]) {
				if (std::find(rejeton->chromT[1].begin(), rejeton->chromT[1].end(), ii) == rejeton->chromT[1].end()) {
					rejeton->chromT[1].push_back(ii);
					for (int scenario = 0; scenario < nbScenario; scenario++) {
						rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = chromLParent2[1 + scenario * paramsList[0]->nbDays][ii];
					}
				}
			}
		}
	}
	if (paramsList[0]->nbDays > 1) {
		for (int scenario = 0; scenario < nbScenario; scenario++) {
			vector<int> garder;
			vector<int> tableauFin;
			vector<vector<bool>> hasBeenInserted = vector<vector<bool>>(paramsList[0]->nbDays + 1, vector<bool>(paramsList[0]->nbClients + paramsList[0]->nbDepots, false));
			int scenarioBegin = scenario * (paramsList[0]->nbDays - 1) + 2;
			vector<int> randomDays = vector<int>(paramsList[0]->nbDays - 1, 0);
			for (int i = 0; i < randomDays.size(); i++) {
				randomDays[i] = scenarioBegin + i;
			}

			std::mt19937 g(paramsList[0]->seed);
			std::shuffle(randomDays.begin(), randomDays.end(), g);

			j1 = paramsList[scenario]->rng->genrand64_int64() % randomDays.size();
			j2 = paramsList[scenario]->rng->genrand64_int64() % randomDays.size();
			if (j1 > j2) std::swap(j1, j2);
			for (unsigned int k = 0; k < randomDays.size(); k++) {
				int day = randomDays[k];
				int dayL = scenario + day;
				int realDay = day - scenario * (paramsList[0]->nbDays - 1);
				// on recopie un segment
				if (k < j1 && !rejeton->chromT[day].empty()) {
					debut = (int)(paramsList[scenario]->rng->genrand64_int64() % rejeton->chromT[day].size());
					fin = (int)(paramsList[scenario]->rng->genrand64_int64() % rejeton->chromT[day].size());
					tableauFin.push_back(fin);
					int j = debut;
					garder.clear();
					while (j != fin) {
						int ii = rejeton->chromT[day][j]; // getting the index to be inherited
						double quantity = chromLParent1[dayL][ii];
						if (quantity > 0.0001) {
							garder.push_back(ii);
							rejeton->chromL[dayL][ii] = quantity;
							hasBeenInserted[realDay][ii] = true;
						}
						j = (j + 1) % rejeton->chromT[day].size();
					}
					rejeton->chromT[day].clear();
					for (int g : garder) {
						rejeton->chromT[day].push_back(g);
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
					for (int ii : rejeton->chromT[day]) {
						garder.push_back(ii);
						rejeton->chromL[dayL][ii] = chromLParent1[dayL][ii];
						hasBeenInserted[realDay][ii] = true;
					}
				}
			}
		
			// completing with rejeton 2
			for (unsigned int k = 0; k < j2; k++) {
				int day = randomDays[k];
				int dayL = day + scenario;
				int realDay = day - scenario * (paramsList[0]->nbDays - 1);
				fin = tableauFin[k];
				for (unsigned int i = 0; i < (int)rejeton2->chromT[day].size(); i++) {
					int ii = rejeton2->chromT[day][(i + fin + 1) % (int)rejeton2->chromT[day].size()];
					if (!hasBeenInserted[realDay][ii]) {
						// computing maximum possible delivery quantity
						double quantity = chromLParent2[dayL][ii];
						if (quantity > 0.0001) {
							rejeton->chromT[day].push_back(ii);
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
	vector<int> garder, joursPerturb, tableauFin;
	int debut, fin, day;
	int j1, j2;
	double quantity;

	// Keeping track of the chromL of the parent
	vector<vector<double>> chromLParent1 = rejeton->chromL;
	vector<vector<double>> chromLParent2 = rejeton2->chromL;

	// Reinitializing the chromL of the rejeton (will become the child)
	// Keeping for each day and each customer the total sum of delivered load and initial inventory
	// (when inserting a customer, need to make sure that we are not exceeding this)
	for (unsigned int k = 1; k <= paramsList[0]->nbDays; k++)
		for (unsigned int i = paramsList[0]->nbDepots; i < paramsList[0]->nbDepots + paramsList[0]->nbClients; i++)
			rejeton->chromL[k][i] = 0.;

	// Keeping a vector to remember if a delivery has already been inserted for on day k for customer i
	vector<vector<bool>> hasBeenInserted = vector<vector<bool>>(paramsList[0]->nbDays + 1, vector<bool>(paramsList[0]->nbClients + paramsList[0]->nbDepots, false));

	// choosing the type of inheritance for each day (nothing, all, or mixed)
	for (int k = 1; k <= paramsList[0]->nbDays; k++)
		joursPerturb.push_back(k);

	std::random_device rd;
	std::mt19937 g(paramsList[0]->seed);
	std::shuffle(joursPerturb.begin(), joursPerturb.end(), g);

	// Picking j1 et j2
	j1 = paramsList[0]->rng->genrand64_int64() % paramsList[0]->nbDays;
	j2 = paramsList[0]->rng->genrand64_int64() % paramsList[0]->nbDays;
	if (j1 > j2) std::swap(j1, j2);

	// Inheriting the data from rejeton1.
	// For each day, we will keep a sequence going from debut to fin
	for (unsigned int k = 0; k < paramsList[0]->nbDays; k++) {
		day = joursPerturb[k];
		// on recopie un segment
		if (k < j1 && !rejeton->chromT[day].empty()) {
			debut = (int)(paramsList[0]->rng->genrand64_int64() % rejeton->chromT[day].size());
			fin = (int)(paramsList[0]->rng->genrand64_int64() % rejeton->chromT[day].size());
			if (debut > fin) std::swap(debut, fin);
			tableauFin.push_back(fin);
			int j = debut;
			garder.clear();
			while (j != ((fin + 1) % rejeton->chromT[day].size())) {
				int ii = rejeton->chromT[day][j]; // getting the index to be inherited
				quantity = chromLParent1[day][ii];
				if (quantity > 0.0001) {
					garder.push_back(ii);
					rejeton->chromL[day][ii] = quantity;
					hasBeenInserted[day][ii] = true;
				}
				j = (j + 1) % rejeton->chromT[day].size();
			}
			rejeton->chromT[day].clear();
			for (int g : garder) {
				rejeton->chromT[day].push_back(g);
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
			for (int ii : rejeton->chromT[day]) {
				garder.push_back(ii);
				rejeton->chromL[day][ii] = chromLParent1[day][ii];
				hasBeenInserted[day][ii] = true;
			}
		}
	}
	
	// completing with rejeton 2
	for (unsigned int k = 0; k < j2; k++) {
		day = joursPerturb[k];
		fin = tableauFin[k];
		for (unsigned int i = 0; i < (int)rejeton2->chromT[day].size(); i++) {
			int ii = rejeton2->chromT[day][(i + fin + 1) % (int)rejeton2->chromT[day].size()];
			if (!hasBeenInserted[day][ii]) {
				// computing maximum possible delivery quantity
				double quantity = chromLParent2[day][ii];
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
