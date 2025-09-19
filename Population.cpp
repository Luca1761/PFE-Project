#include "Population.h"
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>

Population::Population(Params* _params) : params(_params), nbScenario(params->nbScenarios) {
	trainer = new Individual(params);
	valid = new SubPopulation(); 
	invalid = new SubPopulation();

	valid->nbIndiv = 0;
	invalid->nbIndiv = 0;
	capaValidityList = list<bool> (100, true); // arbitrary
	bool count = true;
	
	totalTime = clock();
	
	vector<double> savePenalities = params->penalityCapa;
	for (unsigned int i = 0; i < params->mu * 2; i++) {
		if (i == params->mu) {
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				params->penalityCapa[scenario] *= 50;
			}
			count = false;
		}
		Individual *randomIndiv = new Individual(params);  
		
		train(randomIndiv);
		if (count) updateNbValid(randomIndiv);
		if (addIndividual(randomIndiv) == 0 && randomIndiv->isFeasible) {
			if (params->traces) std::cout << "NEW BEST FEASIBLE FROM INITIALIZATION - Cost: " << randomIndiv->solutionCost.evaluation << std::endl;
		}
		delete randomIndiv;
	}

	params->penalityCapa = savePenalities;
}

Population::~Population() {  // we delete every possible pointers
	for (unsigned int i = 0; i < valid->individuals.size(); i++)
		delete valid->individuals[i];
	for (unsigned int i = 0; i < invalid->individuals.size(); i++)
		delete invalid->individuals[i];
	delete trainer;
}

void Population::evalExtFit(SubPopulation *pop) {
	vector<unsigned int> ranking;
	vector<double> distances;
	for (unsigned int indiv = 0; indiv < pop->nbIndiv; indiv++) {
		ranking.push_back(indiv);
		distances.push_back(pop->individuals[indiv]->distNearest(params->nbCountDistMeasure));
	}

	// ranking of individuals according to their distance measure
	if (pop->nbIndiv > 1) {
		for (unsigned int n = 0; n < pop->nbIndiv; n++) {
			for (unsigned int i = 0; i < pop->nbIndiv - 1 - n; i++) { //i in 0 and 0 or 0 and nbIndiv - 2
				if (distances[ranking[i]] < distances[ranking[i + 1]] - 0.000001)	{
					std::swap(ranking[i], ranking[i + 1]);
				}
			}
		}
	}

	for (unsigned int i = 0; i < pop->nbIndiv; i++) {
		pop->individuals[ranking[i]]->divRank = (float)i / (float)(pop->nbIndiv - 1);
		pop->individuals[ranking[i]]->fitRank = (float)ranking[i] / (float)(pop->nbIndiv - 1);
		pop->individuals[ranking[i]]->extendedFitness = pop->individuals[ranking[i]]->fitRank + ((float)1.0 - (float)params->el / (float)pop->nbIndiv) * pop->individuals[ranking[i]]->divRank;
	}
}

unsigned int Population::addIndividual(Individual *indiv) {
	SubPopulation *subPop = (indiv->isFeasible) ? valid : invalid;
	unsigned int result = placeIndividual(subPop, indiv); // we place the individual in its good subpopulation

	// if two many individuals, remove half of subpopulation
	if (subPop->nbIndiv > params->mu + params->lambda) {
		while (subPop->nbIndiv > params->mu) {
			removeIndividual(subPop, selectCompromise(subPop));
		}
	}
	return result; // return the place of new individual in population (if 0, new best)
}

void Population::updateProximity(SubPopulation *pop, Individual *indiv) {
	for (unsigned int k = 0; k < pop->nbIndiv; k++)
		if (pop->individuals[k] != indiv) {
			pop->individuals[k]->addNearest(indiv); // add new individuals in nearest of others
			indiv->addNearest(pop->individuals[k]); // symetrical
		}
}

bool Population::fitExist(SubPopulation *pop, Individual *indiv) {
	double fitness = indiv->solutionCost.evaluation;
	for (Individual* indiv2 : pop->individuals) {
		if (indiv != indiv2 && indiv2->solutionCost.evaluation >= (fitness - params->delta) && indiv2->solutionCost.evaluation <= (fitness + params->delta))
			return true;
	}
	return false;
}

void Population::diversify() {
	vector<double> savePenalities = params->penalityCapa;

	while (valid->nbIndiv > params->rho * params->mu) { // remove part of feasible individuals
		removeIndividual(valid, valid->nbIndiv - 1);
	}
	while (invalid->nbIndiv > params->rho * params->mu) { // remove part on infeasible individuals
		removeIndividual(invalid, invalid->nbIndiv-1);
	}
	for (unsigned int i = 0; i < params->mu * 2; i++) {
		if (i == params->mu) {
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				params->penalityCapa[scenario] *= 50; // half are not really constrained, the other half is more constrained
			}
		}
		Individual *randomIndiv = new Individual(params);
		train(randomIndiv);
		if (addIndividual(randomIndiv) == 0 && randomIndiv->isFeasible) {
			if (params->traces) std::cout << "NEW BEST FEASIBLE FROM DIVERSIFY - Cost: " << randomIndiv->solutionCost.evaluation << std::endl;
		}
		delete randomIndiv;
	}
	params->penalityCapa = savePenalities;
}

unsigned int Population::placeIndividual(SubPopulation *pop, Individual *indiv) {
	Individual *newIndiv = new Individual(params);
	copyIndividual(newIndiv, indiv);

	unsigned int size = (unsigned int) pop->individuals.size();

	pop->individuals.push_back(newIndiv); // add it at the end (at place size)
	for (unsigned int i = size; i >= 0; i--) { // we try to place it at the good place (pop->individuals is already sorted)
		if (i != 0 && pop->individuals[i - 1]->solutionCost.evaluation >= indiv->solutionCost.evaluation + 0.001) {
			pop->individuals[i] = pop->individuals[i - 1];
		} else {
			pop->individuals[i] = newIndiv;
			pop->nbIndiv++;
			updateProximity(pop, newIndiv);
			if (i == 0 && pop == valid) timeBest = clock();
			return i;
		}
	}
	std::cout << "Failed to place individual" << std::endl;
	throw std::string("Failed to place individual");
	return 1000000; // not supposed to reach this case
}

void Population::removeIndividual(SubPopulation *pop, unsigned int p) {
	Individual *leaving = pop->individuals[p];

	// move the leaving individual at the end
	for (unsigned int i = p + 1; i < pop->individuals.size(); i++) {
		pop->individuals[i - 1] = pop->individuals[i];
	}

	// remove it from the subpopulation
	pop->individuals.pop_back();
	pop->nbIndiv--;

	// remove it from nearest of other individuals of the subpopulation
	for (unsigned int i = 0; i < pop->nbIndiv; i++)
		pop->individuals[i]->removeNearest(leaving);

	// delete the leaving individual
	delete leaving;
}

void Population::validatePen() {
	// penalities have been updated, we also need to update evaluations
	for (Individual* indiv : invalid->individuals) {
		indiv->solutionCost.evaluation = 0.0;
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			indiv->solutionCost_scenario.evaluation[scenario] = indiv->solutionCost_scenario.fitness[scenario] + params->penalityCapa[scenario] * indiv->solutionCost_scenario.capacityViol[scenario]; 
			indiv->solutionCost.evaluation += indiv->solutionCost_scenario.evaluation[scenario];
		}
		indiv->solutionCost.evaluation /= (double) nbScenario;
	}
	for (unsigned int i = 0; i < invalid->nbIndiv; i++) {
		for (unsigned int j = 0; j < invalid->nbIndiv - i - 1; j++) {
			if (invalid->individuals[j]->solutionCost.evaluation > invalid->individuals[j + 1]->solutionCost.evaluation) {
				std::swap(invalid->individuals[j], invalid->individuals[j + 1]);
			}
		}
	}
}

Individual *Population::getIndividualBinT() {
	Individual *individual1;
	Individual *individual2;
	unsigned int place1, place2;

	place1 = (unsigned int) (params->rng->genrand64_int64() % (valid->nbIndiv + invalid->nbIndiv));
	individual1 = (place1 >= valid->nbIndiv) ? invalid->individuals[place1 - valid->nbIndiv] 
											 : valid->individuals[place1];

	place2 = (unsigned int) (params->rng->genrand64_int64() % (valid->nbIndiv + invalid->nbIndiv));
	individual2 = (place2 >= valid->nbIndiv) ? invalid->individuals[place2 - valid->nbIndiv] 
											 : valid->individuals[place2];

	evalExtFit(valid);
	evalExtFit(invalid);

	// the lowest extended fitness wins the binary tournament
	if (individual1->extendedFitness < individual2->extendedFitness) return individual1;
	else return individual2;
}

Individual *Population::getBestFeasIndividual() {
	return (valid->nbIndiv == 0) ? NULL : valid->individuals[0];
}

Individual *Population::getBestInfeasIndividual() {
	return (invalid->nbIndiv == 0) ? NULL : invalid->individuals[0];
}

void Population::copyIndividual(Individual *destination, Individual *source) {
	destination->chromT = source->chromT; // copy of the tour
	destination->chromL = source->chromL; // copy of delivered quantities
	destination->solutionCost = source->solutionCost; // copy of average solution cost
	destination->solutionCost_scenario = source->solutionCost_scenario; // copy of cost, scenario per scenario
	destination->isFeasible = source->isFeasible; //copy of feasibility
}

void Population::ExportPop(string nomFichier, std::vector<double> deliveries, double &totalCost) {
	// exports the current best individual solutions to a folder
	Individual *bestFeasible = getBestFeasIndividual();
	
	if (bestFeasible != NULL) {
		ofstream myfile;
		myfile.precision(10);
		cout.precision(10);
		if (params->jVal == 1) myfile.open(nomFichier.data());
		else myfile.open(nomFichier.data(), std::ios::app); // add on previous (do not erase previous information)
		
		// We will update the local search structure for paths (particularly the split, in case it's not done)
		// We are obliged to set very strong parameters 
		// so that the splitting does not produce a from the best valid solution
		std::vector<double> savePenalities = params->penalityCapa;
		params->penalityCapa = std::vector<double> (params->nbScenarios, 10000.0);
		train(bestFeasible);

		// trainer keeps all the information about bestFeasible routes
		LocalSearch *loc = trainer->localSearchList[0];
		params->penalityCapa = savePenalities;
		
		if (params->jVal == 1) myfile << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
		myfile << "Day "<< params->jVal << endl;
		myfile << endl;
		loc->printInventoryLevels(myfile, deliveries, totalCost);

		// export cost
		myfile << "Total cost: " << totalCost << endl;

		// exporting the number of routes
		int countRoute = 0;
		for (unsigned int i = 0; i < loc->routes[1].size(); i++) {
			countRoute += (!loc->depots[1][i]->next->isADepot);
		}
		myfile << "Route number: " << countRoute << endl;
		myfile << endl;

		// exporting the total CPU time (ms)
		myfile <<"Total Time: " << (float) totalTime / (float)CLOCKS_PER_SEC <<endl; 

		// exporting the time to best solution
		myfile <<"Best Solution Time: " << (float) timeBest / (float)CLOCKS_PER_SEC <<endl; 

		myfile << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
		myfile.close();
		if (params->traces) std::cout << "Successful export" << std::endl;
	} else {
		cout << " Impossible to find a feasible individual " << endl;
	}
}

double Population::validChargePart() {
	int count = 0;
	for (list<bool>::iterator it = capaValidityList.begin(); it != capaValidityList.end(); ++it) {
		count += *it;
	}
	return double(count) / 100.0;
}

double Population::getDiversity(SubPopulation *pop) {
	double total = 0;
	int count = 0;
	for (unsigned int i = 0; i < pop->nbIndiv / 2; i++) {
		for (unsigned int j = i + 1; j < pop->nbIndiv; j++){
			total += pop->individuals[i]->distance(pop->individuals[j]);
			count++;
		}
	}
	return total / (double)count;
}

double Population::getAverageFeasible() {
	double sum = 0.;
	for (Individual* indiv : valid->individuals) {
		sum += indiv->solutionCost.evaluation;
	}
	return sum / (double)valid->nbIndiv;
}

double Population::getAverageInfeasible() {
	double sum = 0.;
	for (Individual* indiv : invalid->individuals) {
		sum += indiv->solutionCost.evaluation;
	}
	return sum / (double)invalid->nbIndiv;
}

unsigned int Population::selectCompromise(SubPopulation *souspop) {
	vector<unsigned int> ranking;

	evalExtFit(souspop);

	// for each individual, we modify the extended fitness
	for (unsigned int i = 0; i < souspop->nbIndiv; i++) {
		ranking.push_back(i);
		if (souspop->individuals[i]->distNearest(1) < params->distMin) souspop->individuals[i]->extendedFitness += 5.0;
		if (fitExist(souspop, souspop->individuals[i])) souspop->individuals[i]->extendedFitness += 5.0;
	}

	// we sort individuals by increasing fitness and we take the worst (highest fitness)
	for (unsigned int n = 0; n < souspop->nbIndiv; n++) {
		for (unsigned int i = 1; i < souspop->nbIndiv - n - 1; i++) {
			if (souspop->individuals[ranking[i]]->extendedFitness > souspop->individuals[ranking[i + 1]]->extendedFitness) {
				std::swap(ranking[i], ranking[i+1]);
			}
		}
	}
	return ranking[souspop->nbIndiv - 1];
}

void Population::train(Individual *indiv) {
	copyIndividual(trainer, indiv); // temporary copy from indiv to trainer
	trainer->split();               // split of big tour in routes + cost measure
	trainer->updateLocalSearch();   // fill local search structures with chromT and chromL
	trainer->runLocalSearch();      // launch local search (including lotsizing)
	trainer->updateIndividual();    // fill back chromT and chromL with local search results
	copyIndividual(indiv, trainer); // copy back from trainer to indiv
}

void Population::updateNbValid(Individual *indiv) {
	capaValidityList.push_back(indiv->solutionCost.capacityViol < 0.0001);
	capaValidityList.pop_front();
}

void Population::measureAndUpdateQuantities(std::vector<double> &deliveries, double &totalCost) {
	if (getBestFeasIndividual() != NULL) {
		Individual* bestIndividual = getBestFeasIndividual();
		bestIndividual->split();
		double dayCost = bestIndividual->measureTrueCost(deliveries); // computes solution cost with true demands
		if (params->traces) std::cout << "Cost of the day " << params->jVal << ": " << dayCost << std::endl;
		totalCost += dayCost;
		if (params->traces) std::cout << "Total cost yet: " << totalCost << std::endl;
    } else {
		std::cout << "NO SOLUTION" << std::endl;
		throw std::string("NO SOLUTION");
    }
    if (params->traces) std::cout << "--------------------------------------------------------------" << std::endl;
}

void Population::displayState(unsigned int nbIter) {
	cout.precision(8);

	cout << "It " << nbIter << " | Sol moy: ";

	if (getBestFeasIndividual() != NULL) {
		cout << getBestFeasIndividual()->solutionCost.evaluation << " Scenarios - ";
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			cout << scenario + 1 << ": " << getBestFeasIndividual()->solutionCost_scenario.evaluation[scenario] << " ";
		}
	} else {
		cout << "NO-VALID ";
	}

	cout << "| ";

	if (getBestInfeasIndividual() != NULL) {
		cout << getBestInfeasIndividual()->solutionCost.evaluation << " Scenarios - ";
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			cout << scenario << ": " << getBestInfeasIndividual()->solutionCost_scenario.evaluation[scenario] << " ";
		}
	} else {
		cout << "NO-INVALID ";
	}

	cout.precision(8);
	double avgPenalityCapa = 0.0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		avgPenalityCapa += params->penalityCapa[scenario];
	}
	avgPenalityCapa /= (double) nbScenario;

	cout << "| Moy " << getAverageFeasible() << " " << getAverageInfeasible()
		 << " | Div " << getDiversity(valid) << " " << getDiversity(invalid)
		 << " | Val " << validChargePart()
		 << " | Pen moy " << avgPenalityCapa << " | Pop " << valid->nbIndiv << " " << invalid->nbIndiv << endl;
}
