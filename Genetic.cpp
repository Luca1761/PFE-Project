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

void Genetic::evolve(int maxIter, int maxIterNonProd, int nbRec)
{
	int place;
	double rangRelatif;
	nbIterNonProd = 1;
	nbIter = 0;
	
	while (nbIter < maxIter && nbIterNonProd < maxIterNonProd && clock() <=  2400 * CLOCKS_PER_SEC)
	{
		// bool time1 = clock();
		// on demande deux individus a la population
		population->recopieIndividu(rejeton, population->getIndividuBinT(rangRelatif));
		population->recopieIndividu(rejeton2, population->getIndividuBinT(rangRelatif));
		
		// on choisit le crossover en fonction du probleme
		crossPOX2();
		// crossPOX_scenario();
		// bool time2 = clock() - time1;
		// std::cout << time2 << std::endl;

		muter_scenario();
		
		// REPAIR IF NEEDED
		if (!rejeton->estValide)
		{
			place = population->addIndividu(rejeton);
			reparer_scenario();
		}
		// ADD IN THE POPULATION IF IS FEASIBLE
		if (rejeton->estValide)
			place = population->addIndividu(rejeton);
		if (place == 0 && rejeton->estValide)
			nbIterNonProd = 1;
		else
			nbIterNonProd++;
		if (place == 0 && rejeton->estValide)
		{
			cout << "NEW BEST FEASIBLE ";
			cout << rejeton->coutSol.evaluation;	
			cout << " Cost : " << rejeton->coutSol.fitness
				 << " capacity Violation : " << rejeton->coutSol.capacityViol;
			cout << endl;
			cout << endl;
		}
		if (nbRec > 0 && nbIterNonProd % (maxIterNonProd / 3 + 1) == maxIterNonProd / 3)
		{
			if (traces)
				cout << "Diversification" << endl;
			population->diversify();
		}
	
		// MANAGEMENT OF PARAMETERS
		if (nbIter % 100 == 0)
			gererPenalites_scenario();
	
		// TRACES
		if (traces && nbIter % 100 == 0)
			population->afficheEtat(nbIter);
		
		nbIter++;

	
	}
	
	// fin de l'algorithme , diverses informations affich�es
	if (traces)
	{
		cout << "time passes : " << clock() << endl;
		cout << "number of iterations : " << nbIter << endl;
	}
	
}

void Genetic::muter_scenario()
{
	rejeton->generalSplit_scenario();
	rejeton->updateLS_scenario();
	rejeton->localSearchRunSearch_scenario();
	rejeton->updateIndiv_scenario();
	population->updateNbValides(rejeton);
}

// eventuellement effectue une reparation de la solution
void Genetic::reparer_scenario()
{
	vector<double> savePenalities(nbScenario, 0.0);

	for (int i = 0; i < nbScenario; i++) {
		savePenalities[i] = paramsList[i]->penalityCapa;
	}

	for (int i = 0; i < nbScenario; i++) {
		if (rejeton->coutSol_scenario.capacityViol[i] > 0.0001)
			paramsList[i]->penalityCapa *= 10;
	}
	if (paramsList[0]->rng->genrand64_real1() < paramsList[0]->pRep)
	{
		rejeton->generalSplit_scenario();
		rejeton->updateLS_scenario();
		rejeton->localSearchRunSearch_scenario();
		rejeton->updateIndiv_scenario();
		if (!rejeton->estValide) {
			for (int i = 0; i < nbScenario; i++) {
				if (rejeton->coutSol_scenario.capacityViol[i] > 0.0001)
					paramsList[i]->penalityCapa *= 500;
			}
			rejeton->generalSplit_scenario();
			rejeton->updateLS_scenario();
			rejeton->localSearchRunSearch_scenario();
			rejeton->updateIndiv_scenario();
		}
	}

	for (int i = 0; i < nbScenario; i++) {
		paramsList[i]->penalityCapa = savePenalities[i];
	}
}

// gestion des penalites
void Genetic::gererPenalites_scenario()
{
	double fractionCharge = population->fractionValidesCharge();

	for (int i = 0; i < nbScenario; i++) {
		if (fractionCharge < paramsList[i]->minValides && paramsList[i]->penalityCapa < 1000)
			paramsList[i]->penalityCapa = paramsList[i]->penalityCapa * 1.2;
		else if (fractionCharge > paramsList[i]->maxValides && paramsList[i]->penalityCapa > 0.01)
			paramsList[i]->penalityCapa = paramsList[i]->penalityCapa * 0.85;
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
	bool traces = false;
    vector<int> vide, garder, tableauFin, tableauEtat;
    vector<double> charge;
    int debut, fin, day;
    int j1, j2;
    double quantity;

	// Keeping track of the chromL of the parent
	vector<vector<int>> chromTParent1 = rejeton->chromT;
	vector<vector<int>> chromTParent2 = rejeton2->chromT;

	for (int k = 1; k < rejeton->chromL.size(); k++) {
		for (int i = 0; i < rejeton->chromL[k].size(); i++) {
			rejeton->chromL[k][i] = 0.;
		}
	}

	// 接着对第一天进行判断，继承1还是2，还是各自部分继承，选一个随机数0-2取整
	int firstDayInheritance = paramsList[0]->rng->genrand64_int64() % 4;
	if (firstDayInheritance == 0) {
		rejeton->chromT[1].clear();
		if (chromTParent1[1].size() != 0) {
			for (int i = 0; i < chromTParent1[1].size(); i++) {
				if (chromTParent1[1][i]) {
					int ii = chromTParent1[1][i];
					rejeton->chromT[1].push_back(ii);
					for (int scenario = 0; scenario < nbScenario; scenario++) {
						rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = 1.;
					}
				}
			}
		}
	}
	else if (firstDayInheritance == 1) {
		rejeton->chromT[1].clear();
		if (chromTParent2[1].size() != 0) {
			for (int i = 0; i < (int)chromTParent2[1].size(); i++) {
				if (chromTParent2[1][i]) {
					int ii = chromTParent2[1][i];
					rejeton->chromT[1].push_back(ii);
					for (int scenario = 0; scenario < nbScenario; scenario++) {
						rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = 1.;
					}
				}
			}
		}
	}
	else {
		rejeton->chromT[1].clear();
		if (chromTParent1[1].size() != 0) {
			int begin = paramsList[0]->rng->genrand64_int64() % (chromTParent1[1].size());
			int end = paramsList[0]->rng->genrand64_int64() % chromTParent1[1].size();
			int i = begin;
			while (i % chromTParent1[1].size() != end) {	
				if (chromTParent1[1][i]) {
					int ii = chromTParent1[1][i];
					rejeton->chromT[1].push_back(ii);
					for (int scenario = 0; scenario < nbScenario; scenario++) {
						rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = 1.;
					}
				}
				i = (i + 1) % chromTParent1[1].size();
			}
		}
		
		// cout << "ok2" << endl;
		// 之后继承2的部分，遍历rejeton2的chromT，如果未被插入，则插入rejeton
		if (chromTParent2[1].size() != 0) {
			for (int i = 0; i < chromTParent2[1].size(); i++) {
				if (chromTParent2[1][i]) {
					int ii = chromTParent2[1][i];
					if (std::find(rejeton->chromT[1].begin(), rejeton->chromT[1].end(), ii) == rejeton->chromT[1].end()) {
						rejeton->chromT[1].push_back(ii);
						for (int scenario = 0; scenario < nbScenario; scenario++) {
							rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = 1.;
						}
					}
				}
			}
		}
	}

	for (int scenario = 0; scenario < nbScenario; scenario++) {
		int scenarioBegin = scenario * (paramsList[0]->nbDays-1) + 2;
		int scenarioEnd = scenarioBegin + paramsList[0]->nbDays - 1;
		vector<int> randomDays = vector<int>(paramsList[0]->nbDays-1, 0);
		for (int i = 0; i < randomDays.size(); i++) {
			randomDays[i] = scenarioBegin + i;
		}

		std::mt19937 g(paramsList[0]->seed);
		std::shuffle(randomDays.begin(), randomDays.end(), g);
		// 此时randomDays中是对应scenario的天数的随机序列，之后对于每一天
		j1 = paramsList[0]->rng->genrand64_int64() % randomDays.size();
		j2 = paramsList[0]->rng->genrand64_int64() % randomDays.size();
		if (j1 > j2) std::swap(j1, j2);
		for (int k = 0; k < randomDays.size(); k++) {
			int day = randomDays[k];
			// 我们会跳过所有需要继承2的部分，在后面用一个循环统一插入
			if (k < j1 && !chromTParent1[day].empty()) {
				// 完全继承1
				if(traces) cout << "start inherit 1!" << endl;
				rejeton->chromT[day].clear();
				if (chromTParent1[day].size() != 0) {
					for (int i = 0; i < chromTParent1[day].size(); i++) {	
						if (chromTParent1[day][i]) {
							int ii = chromTParent1[day][i];
							rejeton->chromT[day].push_back(ii);
							for (int scenario = 0; scenario < nbScenario; scenario++) {
								rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = 1.;
							}
						}
					}
				}
			}
			else if (k < j2) {
				continue;
			}
			else {
				if (chromTParent1[day].empty()) {
					continue;
				} else {
					}
				
				int begin = paramsList[0]->rng->genrand64_int64() % chromTParent1[day].size();
				int end = paramsList[0]->rng->genrand64_int64() % chromTParent1[day].size();
				int i = begin;
				rejeton->chromT[day].clear();
				while (i % chromTParent1[day].size() != end) {
					if (chromTParent1[day][i]) {
						int ii = chromTParent1[day][i];
						if (std::find(rejeton->chromT[day].begin(), rejeton->chromT[day].end(), ii) == rejeton->chromT[day].end()) {
							rejeton->chromT[day].push_back(ii);
							for (int scenario = 0; scenario < nbScenario; scenario++) {
								rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = 1.;
							}
						}
					}
					i = (i + 1) % chromTParent1[day].size();
				}
			}
		}
		// 输出完整的chromT
		if(traces) cout << "First inherit done!scenario = " << scenario << endl;
		// 之后对当前scenario所有天的2继承进行统一插入
		for (int k = 0; k < randomDays.size(); k++) {
			int day = randomDays[k];
			if (rejeton2->chromT[day].size() != 0) {
				for (int i = 0; i < rejeton2->chromT[day].size(); i++) {	
					if (rejeton2->chromT[day][i]) {
						int ii = rejeton2->chromT[day][i];
						if (std::find(rejeton->chromT[day].begin(), rejeton->chromT[day].end(), ii) == rejeton->chromT[day].end()) {
							rejeton->chromT[day].push_back(ii);
							for (int scenario = 0; scenario < nbScenario; scenario++) {
								rejeton->chromL[1 + scenario * paramsList[0]->nbDays][ii] = 1.;
							}
						}
					}
				}
			}
		}
	}
	// 之后根据库存推导出每一天的默认配送量（应该遵循OU Policy）
	// 对于任何一个scenario 都需要经过第一天，之后找到自己的index范围进行推算
	for (int scenario = 0; scenario < nbScenario; scenario++) {
		vector<vector<double>> I_end(paramsList[0]->nbDays+2, vector<double>(paramsList[0]->nbDepots + paramsList[0]->nbClients));
		vector<int> dayIndex = vector<int>(paramsList[0]->nbDays+1, 0);
		int startIndex = scenario * (paramsList[0]->nbDays) + 1;

		vector<int> dayIndexL(paramsList[0]->nbDays + 1, 0);
		for (int k = 1; k <= paramsList[0]->nbDays; k++){
			dayIndex[k] = startIndex + k - 1;
		}



		for (int i = paramsList[0]->nbDepots; i < paramsList[0]->nbDepots + paramsList[0]->nbClients; i++){
			I_end[0][i] = paramsList[0]->cli[i].startingInventory;
		}
		for (int i = 1; i <= paramsList[0]->nbDays; i++) {
			int day = dayIndex[i];
			for (int cus = paramsList[0]->nbDepots; cus < paramsList[0]->nbDepots + paramsList[0]->nbClients; cus++) {
				if (rejeton->chromL[day][cus] == 1) {
					rejeton->chromL[day][cus] = paramsList[0]->cli[cus].maxInventory - I_end[i-1][cus];
					I_end[i][cus] = std::max<double>(0., paramsList[0]->cli[cus].maxInventory - paramsList[0]->cli[cus].dailyDemand[day]);
				}
				else {
					rejeton->chromL[day][cus] = 0.;
					I_end[i][cus] = std::max<double>(0., I_end[i-1][cus] - paramsList[0]->cli[cus].dailyDemand[day]);
				}
			}
		}
	}
	rejeton->generalSplit_scenario();
	return 0;
}

void Genetic::crossPOX2()
{
	vector<int> garder, joursPerturb, tableauFin;
	int debut, fin, day;
	int j1, j2;
	double quantity;

	// Keeping track of the chromL of the parent
	vector<vector<double>> chromLParent1 = rejeton->chromL;

	// Reinitializing the chromL of the rejeton (will become the child)
	// Keeping for each day and each customer the total sum of delivered load and initial inventory
	// (when inserting a customer, need to make sure that we are not exceeding this)
	for (int k = 1; k <= paramsList[0]->nbDays; k++)
		for (int i = paramsList[0]->nbDepots; i < paramsList[0]->nbDepots + paramsList[0]->nbClients; i++)
			rejeton->chromL[k][i] = 0.;

	// Keeping a vector to remember if a delivery has already been inserted for on day k for customer i
	vector<vector<bool>> hasBeenInserted = vector<vector<bool>>(paramsList[0]->nbDays + 1, vector<bool>(paramsList[0]->nbClients + paramsList[0]->nbDepots, false));

	// choosing the type of inheritance for each day (nothing, all, or mixed)
	for (int k = 1; k <= paramsList[0]->nbDays; k++)
		joursPerturb.push_back(k);

	// std::random_shuffle(joursPerturb.begin(), joursPerturb.end());
	std::random_device rd;
	std::mt19937 g(paramsList[0]->seed);
	std::shuffle(joursPerturb.begin(), joursPerturb.end(), g);

	// Picking j1 et j2
	j1 = paramsList[0]->rng->genrand64_int64() % paramsList[0]->nbDays;
	j2 = paramsList[0]->rng->genrand64_int64() % paramsList[0]->nbDays;
	if (j1 > j2)
	{
		int temp = j2;
		j2 = j1;
		j1 = temp;
	}

	// Inheriting the data from rejeton1.
	// For each day, we will keep a sequence going from debut to fin
	for (int k = 0; k < paramsList[0]->nbDays; k++)
	{
		day = joursPerturb[k];
		// on recopie un segment
		if (k < j1 && !rejeton->chromT[day].empty())
		{
			debut = (int)(paramsList[0]->rng->genrand64_int64() % rejeton->chromT[day].size());
			fin = (int)(paramsList[0]->rng->genrand64_int64() % rejeton->chromT[day].size());
			tableauFin.push_back(fin);
			int j = debut;
			garder.clear();
			while (j != ((fin + 1) % rejeton->chromT[day].size()))
			{
				int ii = rejeton->chromT[day][j]; // getting the index to be inherited
				garder.push_back(ii);
				rejeton->chromL[day][ii] = chromLParent1[day][ii];
				hasBeenInserted[day][ii] = true;
				j = (j + 1) % rejeton->chromT[day].size();
			}
			rejeton->chromT[day].clear();
			for (int g = 0; g < (int)garder.size(); g++)
				rejeton->chromT[day].push_back(garder[g]);
		}
		else if (k < j2) // on recopie rien
		{
			rejeton->chromT[day].clear();
			tableauFin.push_back(-1);
		}
		else // on recopie tout
		{
			tableauFin.push_back(0);
			for (int j = 0; j < (int)rejeton->chromT[day].size(); j++)
			{
				int ii = rejeton->chromT[day][j]; // getting the index to be inherited
				garder.push_back(ii);
				rejeton->chromL[day][ii] = chromLParent1[day][ii];
				hasBeenInserted[day][ii] = true;
				j = (j + 1) % rejeton->chromT[day].size();
			}
		}
	}
	
	// completing with rejeton 2
	for (int k = 0; k < paramsList[0]->nbDays; k++)
	{
		day = joursPerturb[k];
		fin = tableauFin[k];
		if (k < j2)
		{
			for (int i = 0; i < (int)rejeton2->chromT[day].size(); i++)
			{
				int ii = rejeton2->chromT[day][(i + fin + 1) % (int)rejeton2->chromT[day].size()];
				if (!hasBeenInserted[day][ii]) // it has not been inserted yet
				{
					// computing maximum possible delivery quantity
					quantity = std::min<double>(rejeton->maxFeasibleDeliveryQuantity(day, ii), rejeton2->chromL[day][ii]);
					if (quantity > 0.0001)
					{
						rejeton->chromT[day].push_back(ii);
						rejeton->chromL[day][ii] = quantity;
						hasBeenInserted[day][ii] = true;
					}
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
