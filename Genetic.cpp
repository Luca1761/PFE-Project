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
		// on demande deux individus a la population
		population->recopieIndividu(rejeton, population->getIndividuBinT(rangRelatif));
		population->recopieIndividu(rejeton2, population->getIndividuBinT(rangRelatif));
		
		// on choisit le crossover en fonction du probleme
		crossPOX2();

		muter();
		
		// REPAIR IF NEEDED
		if (!rejeton->estValide)
		{
			place = population->addIndividu(rejeton);
			reparer();
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
			gererPenalites();
	
		// TRACES
		if (traces && nbIter % 500 == 0)
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

// effectue la mutation
void Genetic::muter()
{
	rejeton->updateLS_scenario();
	rejeton->localSearchList[0]->runSearchTotal(false);
	
	rejeton->updateIndiv_scenario();
	population->updateNbValides(rejeton);
}

// eventuellement effectue une reparation de la solution
void Genetic::reparer()
{
	double temp, temp2;
	bool continuer = false;

	temp = paramsList[0]->penalityCapa;

	/*First tentative*/
	paramsList[0]->penalityCapa *= 10;
	if (paramsList[0]->rng->genrand64_real1() < paramsList[0]->pRep)
	{
		rejeton->updateLS_scenario();
		rejeton->localSearchList[0]->runSearchTotal(true);
		rejeton->updateIndiv_scenario();

		/* Second tentative*/
		if (!rejeton->estValide)
		{
			paramsList[0]->penalityCapa *= 500;
			rejeton->generalSplit_scenario();
			rejeton->updateLS_scenario();
			rejeton->localSearchList[0]->runSearchTotal(true);
			rejeton->updateIndiv_scenario();
		}
	}
	paramsList[0]->penalityCapa = temp;
}

// gestion des penalites
void Genetic::gererPenalites()
{
	double fractionCharge = population->fractionValidesCharge();

	if (fractionCharge < paramsList[0]->minValides && paramsList[0]->penalityCapa < 1000)
		paramsList[0]->penalityCapa = paramsList[0]->penalityCapa * 1.2;
	else if (fractionCharge > paramsList[0]->maxValides && paramsList[0]->penalityCapa > 0.01)
		paramsList[0]->penalityCapa = paramsList[0]->penalityCapa * 0.85;

	population->validatePen();
}

Genetic::Genetic(vector<Params*> pl, Population *population, clock_t ticks, bool traces) : paramsList(pl), population(population), ticks(ticks), traces(traces)
{
	rejeton = new Individu(paramsList);
	rejeton2 = new Individu(paramsList);
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
