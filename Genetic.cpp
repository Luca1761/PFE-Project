#include "Genetic.h"
#include <algorithm> 
#include <unistd.h>

void Genetic::evolve(unsigned int maxIter, unsigned int maxIterNonProd, unsigned int maxTime) {
	unsigned int place;
	unsigned int diversificationIdx = 1;
	nbIter = 0;
	nbIterNonProd = 1;
	while (nbIter < maxIter && nbIterNonProd < maxIterNonProd && (round((float) (clock() - population->totalTime) / CLOCKS_PER_SEC) <  maxTime)) {
		// we take two individuals from the population
		population->copyIndividual(child, population->getIndividualBinT());
		population->copyIndividual(child2, population->getIndividualBinT());

		// we choose the crossover
		if (nbScenario == 1) crossPOX();
		else crossPOX_scenario();

		// apply the local search
		mutate();
		
		// repair if needed
		if (!child->isFeasible) {
			place = population->addIndividual(child);
			repair();
		}

		// add feasible solutions (repaired or not)
		if (child->isFeasible) { 
			place = population->addIndividual(child);
		}
		
		nbIterNonProd++;
		if (place == 0 && child->isFeasible) {
			nbIterNonProd = 1;
			if (child->solution_cost.capacityViol > 0.0001) {
				std::cout << "SOLUTION IS NOT FEASIBLE" << std::endl;
				throw std::string("SOLUTION IS NOT FEASIBLE");
			}
			if (params->traces) {
				std::cout << "NEW BEST FEASIBLE - Cost: " << child->solution_cost.evaluation << std::endl;
				std::cout << endl;
			}
		}
		if (nbIterNonProd / 500 == diversificationIdx) {
			if (params->traces) std::cout << "Diversification" << endl;
			population->diversify();
			diversificationIdx++;
		}
	
		// MANAGEMENT OF PARAMETERS
		if (nbIter % 200 == 0) managePenalties();
	
		// TRACES
		if (nbIter % 50 == 0 && params->traces) population->displayState(nbIter);
		
		nbIter++;
	}
	
	// end of the algorithm, we display total time, time to find the best solution and total number of iterations
	population->timeBest = population->timeBest - population->totalTime;
	population->totalTime = clock() - population->totalTime;
	if (params->traces) {
		std::cout << "time passes : " << (float) (population->totalTime)/ CLOCKS_PER_SEC << endl;
		std::cout << "time to find best solution : " << (float) (population->timeBest)/ CLOCKS_PER_SEC << endl;
		std::cout << "number of iterations : " << nbIter << endl;
	}
}

void Genetic::mutate() {
	child->split();                     // split the big tour in chromT 
	child->updateLocalSearch();			// fill the localsearch structure with tour information
	child->runLocalSearch();			// launch local search
	child->updateIndividual();			// update the big tour using result from local search
	population->updateNbValid(child); 	// update the number of valid individuals
}

void Genetic::repair() {
	vector<double> savePenalities = params->penalityCapa; // we save penalities to restore them just after

	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (child->solutionCost_scenario.capacityViol[scenario] > 0.0001)
			params->penalityCapa[scenario] *= 10; // increase penalties to force feasibility
	}
	if (params->rng->genrand64_real1() < params->pRep) {
		child->split();
		child->updateLocalSearch();
		child->runLocalSearch();
		child->updateIndividual();
		if (!child->isFeasible) {
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				if (child->solutionCost_scenario.capacityViol[scenario] > 0.0001)
					params->penalityCapa[scenario] *= 500;
			}
			child->split();
			child->updateLocalSearch();
			child->runLocalSearch();
			child->updateIndividual();
		}
	}
	params->penalityCapa = savePenalities;
}

void Genetic::managePenalties() {
	double chargePart = population->validChargePart();
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		if (chargePart < params->minFeasibles && params->penalityCapa[scenario] < 1000)
			params->penalityCapa[scenario] *= 1.25;
		else if (chargePart > params->maxFeasibles && params->penalityCapa[scenario] > 0.01)
			params->penalityCapa[scenario] *= 0.8;
	}
	population->validatePen();
}

Genetic::Genetic(Params* _params, Population* _population) : params(_params), population(_population), nbScenario(params->nbScenarios) {
	child = new Individual(params);
	child2 = new Individual(params);
}

void Genetic::crossPOX_scenario() {
    unsigned int j1, j2, l;
	unsigned int begin, end;
	unsigned int dayT, dayL, realDay;

	// Keeping track of the chromT of the parent
	vector<vector<unsigned int>> chromTParent1 = child->chromT;
	vector<vector<double>> chromLParent1 = child->chromL;
	vector<vector<unsigned int>> chromTParent2 = child2->chromT;
	vector<vector<double>> chromLParent2 = child2->chromL;

	for (unsigned int k = 1; k < child->chromL.size(); k++) {
		for (unsigned int i = 0; i < child->chromL[k].size(); i++) {
			child->chromL[k][i] = 0.;
		}
	}

	int firstDayInheritance = params->rng->genrand64_int64() % 4;
	child->chromT[1].clear();
	if (firstDayInheritance == 0) {
		for (unsigned int ii : chromTParent1[1]) {
			child->chromT[1].push_back(ii);
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				child->chromL[1 + scenario * params->nbDays][ii] = chromLParent1[1 + scenario * params->nbDays][ii];
			}
		}
	} else if (firstDayInheritance == 1) {
		for (unsigned int ii : chromTParent2[1]) {
			child->chromT[1].push_back(ii);
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				child->chromL[1 + scenario * params->nbDays][ii] = chromLParent2[1 + scenario * params->nbDays][ii];
			}
		}
	} else {
		if (!chromTParent1[1].empty()) {
			begin = (unsigned int) (params->rng->genrand64_int64() % chromTParent1[1].size());
			end = (unsigned int) (params->rng->genrand64_int64() % chromTParent1[1].size());
			l = begin;
			while (l != end) {
				unsigned int ii = chromTParent1[1][l];
				child->chromT[1].push_back(ii);
				for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
					child->chromL[1 + scenario * params->nbDays][ii] = chromLParent1[1 + scenario * params->nbDays][ii];
				}
				l = (unsigned int) ((l + 1) % chromTParent1[1].size());
			}
		}
		
		if (!chromTParent2[1].empty()) {
			for (unsigned int ii : chromTParent2[1]) {
				if (std::find(child->chromT[1].begin(), child->chromT[1].end(), ii) == child->chromT[1].end()) {
					child->chromT[1].push_back(ii);
					for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
						child->chromL[1 + scenario * params->nbDays][ii] = chromLParent2[1 + scenario * params->nbDays][ii];
					}
				}
			}
		}
	}
	if (params->nbDays > 1) {
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			vector<unsigned int> keep;
			vector<int> endTable;
			vector<vector<bool>> hasBeenInserted = vector<vector<bool>>(params->nbDays + 1, vector<bool>(params->nbClients + params->nbDepots, false));
			vector<unsigned int> randomDays = vector<unsigned int>(params->nbDays - 1, 0);

			for (unsigned int i = 0; i < randomDays.size(); i++) {
				randomDays[i] = scenario * (params->nbDays - 1) + 2 + i;
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

				// we copy a part of child1 solution
				if (k < j1 && !child->chromT[dayT].empty()) {
					begin = (unsigned int) (params->rng->genrand64_int64() % child->chromT[dayT].size());
					end = (unsigned int) (params->rng->genrand64_int64() % child->chromT[dayT].size());
					endTable.push_back((int) end);
					unsigned int j = begin;
					keep.clear();
					while (j != end) {
						unsigned int ii = child->chromT[dayT][j]; // getting the index to be inherited
						double quantity = chromLParent1[dayL][ii];
						if (quantity > 0.0001) {
							keep.push_back(ii);
							child->chromL[dayL][ii] = quantity;
							hasBeenInserted[realDay][ii] = true;
						}
						j = (unsigned int) ((j + 1) % child->chromT[dayT].size());
					}
					child->chromT[dayT].clear();
					for (unsigned int keptClient : keep) {
						child->chromT[dayT].push_back(keptClient);
					}
				} else if (k < j2) { // copy nothing
					child->chromT[dayT].clear();
					endTable.push_back(-1);
				} else { // copy everything
					endTable.push_back(0);
					for (unsigned int ii : child->chromT[dayT]) {
						keep.push_back(ii);
						child->chromL[dayL][ii] = chromLParent1[dayL][ii];
						hasBeenInserted[realDay][ii] = true;
					}
				}
			}
		
			// completing with child2
			for (unsigned int k = 0; k < j2; k++) {
				dayT = randomDays[k];
				dayL = dayT + scenario;
				realDay = dayT - scenario * (params->nbDays - 1);
				for (unsigned int i = 0; i < child2->chromT[dayT].size(); i++) {
					end = (unsigned int) ((int) i + endTable[k] + 1);
					unsigned int ii = child2->chromT[dayT][end % child2->chromT[dayT].size()];
					if (!hasBeenInserted[realDay][ii]) {
						double quantity = chromLParent2[dayL][ii];
						if (quantity > 0.0001) {
							child->chromT[dayT].push_back(ii);
							child->chromL[dayL][ii] = quantity;
							hasBeenInserted[realDay][ii] = true;
						}
					}
				}
			}
		}
	}
	child->split();
	return;
}

void Genetic::crossPOX() {
	vector<int> endTable;
	vector<unsigned int> keep;
	vector<unsigned int> randomDays;
	unsigned int begin, end, day;
	unsigned int j1, j2;
	double quantity;

	// Keeping track of the chromL of the parent
	vector<vector<double>> chromLParent1 = child->chromL;
	vector<vector<double>> chromLParent2 = child2->chromL;

	// Reinitializing the chromL of the child (will become the new child)
	// Keeping for each day and each customer the total sum of delivered load and initial inventory
	// (when inserting a customer, need to make sure that we are not exceeding this)
	for (unsigned int k = 1; k <= params->nbDays; k++)
		for (unsigned int i = params->nbDepots; i < params->nbDepots + params->nbClients; i++)
			child->chromL[k][i] = 0.;

	// Keeping a vector to remember if a delivery has already been inserted for on day k for customer i
	vector<vector<bool>> hasBeenInserted = vector<vector<bool>>(params->nbDays + 1, vector<bool>(params->nbClients + params->nbDepots, false));

	for (unsigned int k = 1; k <= params->nbDays; k++)
		randomDays.push_back(k);

	std::random_device rd;
	std::mt19937 g((unsigned int) params->seed);
	std::shuffle(randomDays.begin(), randomDays.end(), g);

	// Picking j1 et j2
	j1 = (unsigned int) (params->rng->genrand64_int64() % params->nbDays);
	j2 = (unsigned int) (params->rng->genrand64_int64() % params->nbDays);
	if (j1 > j2) std::swap(j1, j2);

	// Inheriting the data from child1.
	// For each day, we will keep a sequence going from begin to end
	for (unsigned int k = 0; k < params->nbDays; k++) {
		day = randomDays[k];
		// we copy a part
		if (k < j1 && !child->chromT[day].empty()) {
			begin = (unsigned int)(params->rng->genrand64_int64() % child->chromT[day].size());
			end = (unsigned int)(params->rng->genrand64_int64() % child->chromT[day].size());
			if (begin > end) std::swap(begin, end);
			endTable.push_back((int)end);
			unsigned int j = begin;
			keep.clear();
			while (j != ((end + 1) % child->chromT[day].size())) {
				unsigned int ii = child->chromT[day][j]; // getting the index to be inherited
				quantity = chromLParent1[day][ii];
				if (quantity > 0.0001) {
					keep.push_back(ii);
					child->chromL[day][ii] = quantity;
					hasBeenInserted[day][ii] = true;
				}
				j = (unsigned int) ((j + 1) % child->chromT[day].size());
			}
			child->chromT[day].clear();
			for (unsigned int keptClient : keep) {
				child->chromT[day].push_back(keptClient);
			}
		}
		else if (k < j2) { // copy nothing
			child->chromT[day].clear();
			endTable.push_back(-1);
		} else { // copy everything
			endTable.push_back(0);
			for (unsigned int ii : child->chromT[day]) {
				keep.push_back(ii);
				child->chromL[day][ii] = chromLParent1[day][ii];
				hasBeenInserted[day][ii] = true;
			}
		}
	}
	
	// completing with child2
	for (unsigned int k = 0; k < j2; k++) {
		day = randomDays[k];
		for (unsigned int i = 0; i < child2->chromT[day].size(); i++) {
			end = (unsigned int) ((int) i + endTable[k] + 1);
			unsigned int ii = child2->chromT[day][end % child2->chromT[day].size()];
			if (!hasBeenInserted[day][ii]) {
				quantity = chromLParent2[day][ii];
				if (quantity > 0.0001) {
					child->chromT[day].push_back(ii);
					child->chromL[day][ii] = quantity;
					hasBeenInserted[day][ii] = true;
				}
			}
		}
	}
	child->split();
}

Genetic::~Genetic(void) {
	delete child;
	delete child2;
}
