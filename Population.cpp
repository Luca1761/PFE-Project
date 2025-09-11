#include "Population.h"
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
// constructeur

Population::Population(Params* _params) : params(_params)
{
	nbScenario = params->nbScenarios;

	trainer = new Individual(params);
	valid = new SousPop(); 
	invalid = new SousPop();

	valid->nbIndiv = 0;
	invalid->nbIndiv = 0;
	listeValiditeCharge = list<bool> (100, true);
	bool compter = true;
	
	vector<double> saveCapa(nbScenario);

	totalTime = clock();
	
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++)
		saveCapa[scenario] = params->penalityCapa[scenario];
	
	for (unsigned int i = 0; i < params->mu * 2; i++) {
		if (i == params->mu) {
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				params->penalityCapa[scenario] *= 50;
			}
			compter = false;
		}
		Individual *randomIndiv = new Individual(params);  
		
		education_scenario(randomIndiv);
		if (compter) updateNbValid(randomIndiv);
		if (addIndividu(randomIndiv) == 0 && randomIndiv->isFeasible) {
			if (params->traces) std::cout << "NEW BEST FEASIBLE FROM INITIALIZATION" << std::endl;
		}
		delete randomIndiv;
	}
	
	// on initialise par defaut a 100, comme si tout etait valide au debut
	// mais c'est arbitraire

	for (unsigned int scenario = 0; scenario < (unsigned int) nbScenario; scenario++) {
		params->penalityCapa[scenario] = saveCapa[scenario];
	}
}

// destructeur
Population::~Population() {
	for (unsigned int i = 0; i < valid->individus.size(); i++)
		delete valid->individus[i];
	for (unsigned int i = 0; i < invalid->individus.size(); i++)
		delete invalid->individus[i];
	delete trainer;
}

void Population::evalExtFit(SousPop *pop) {
	vector<unsigned int> classement;
	vector<double> distances;
	for (unsigned int indiv = 0; indiv < pop->nbIndiv; indiv++) {
		classement.push_back(indiv);
		distances.push_back(pop->individus[indiv]->distNearest(params->nbCountDistMeasure));
	}

	// classement des individus en fonction de leur note de distance
	for (unsigned int n = 0; n < pop->nbIndiv; n++) {
		for (unsigned int i = 0; i < pop->nbIndiv - 1 - n; i++) { //i in 0 and 0 or 0 and nbIndiv - 2
			if (distances[classement[i]] < distances[classement[i + 1]] - 0.000001)	{
				std::swap(classement[i], classement[i + 1]);
			}
		}
	}

	for (unsigned int i = 0; i < pop->nbIndiv; i++) {
		pop->individus[classement[i]]->divRank = (float)i / (float)(pop->nbIndiv - 1);
		pop->individus[classement[i]]->fitRank = (float)classement[i] / (float)(pop->nbIndiv - 1);
		pop->individus[classement[i]]->extendedFitness = pop->individus[classement[i]]->fitRank + ((float)1.0 - (float)params->el / (float)pop->nbIndiv) * pop->individus[classement[i]]->divRank;
	}
}

unsigned int Population::addIndividu(Individual *indiv) {
	SousPop *souspop = (indiv->isFeasible) ? valid : invalid;
	unsigned int result = placeIndividu(souspop, indiv);

	// il faut eventuellement enlever la moitie de la pop
	if (souspop->nbIndiv > params->mu + params->lambda) {
		while (souspop->nbIndiv > params->mu) {
			removeIndividu(souspop, selectCompromis(souspop));
		}
	}
	return result;
}

// met a jour les individus les plus proches d'une population
// en fonction de l'arrivant

void Population::updateProximity(SousPop *pop, Individual *indiv) {
	for (unsigned int k = 0; k < pop->nbIndiv; k++)
		if (pop->individus[k] != indiv) {
			pop->individus[k]->addNearest(indiv);
			indiv->addNearest(pop->individus[k]);
		}
}

// fonction booleenne verifiant si le fitness n'existe pas d�ja

bool Population::fitExist(SousPop *pop, Individual *indiv)
{
	double fitness = indiv->solutionCost.evaluation;
	for (Individual* indiv2 : pop->individus) {
		if (indiv != indiv2 && indiv2->solutionCost.evaluation >= (fitness - params->delta) && indiv2->solutionCost.evaluation <= (fitness + params->delta))
			return true;
	}
	return false;
}
// procede de redemarrage avec remplacement d'une partie de la population
// modele tres simplifie
// on remplace la moitie individus de fitness situes sous la moyenne par des individus aleatoires
void Population::diversify() {
	vector<double> savePenalities(nbScenario, 0.0);

	for (unsigned int scenario = 0; scenario < (unsigned int) nbScenario; scenario++) {
		savePenalities[scenario] = params->penalityCapa[scenario];
	}
	
	while (valid->nbIndiv > params->rho * params->mu) {
		removeIndividu(valid, valid->nbIndiv - 1);
	}
	while (invalid->nbIndiv > params->rho * params->mu) {
		removeIndividu(invalid, invalid->nbIndiv-1);
	}
	for (unsigned int i = 0; i < params->mu * 2; i++) {
		if (i == params->mu) {
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				params->penalityCapa[scenario] *= 50;
			}
		}
		Individual *random = new Individual(params);
		education_scenario(random);
		if (addIndividu(random) == 0 && random->isFeasible) {
			if (params->traces) std::cout << "NEW BEST FEASIBLE FROM DIVERSIFY" << std::endl;
		}
		delete random;
	}
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		params->penalityCapa[scenario] = savePenalities[scenario];
	}
}

//Try to place individual indiv in the sub population and return its position 
//Individuals are increasingly sorted in the subpopulation
unsigned int Population::placeIndividu(SousPop *pop, Individual *indiv) {
	Individual *monIndiv = new Individual(params);
	copyIndividual(monIndiv, indiv);

	// regarde si son fitness est suffisamment espace
	bool placed = false;
	unsigned int i = (unsigned int) pop->individus.size();

	pop->individus.push_back(monIndiv);
	while (i >= 1 && !placed) {
		if (pop->individus[i - 1]->solutionCost.evaluation >= indiv->solutionCost.evaluation + 0.001) {
			pop->individus[i] = pop->individus[i - 1];
			i--;
		} else {
			pop->individus[i] = monIndiv;
			placed = true;
			pop->nbIndiv++;
			updateProximity(pop, monIndiv);
			return i; // reussite
		}
	}
	if (!placed) {
		pop->individus[0] = monIndiv;
		placed = true;
		pop->nbIndiv++;
		updateProximity(pop, monIndiv);
		if (pop == valid) timeBest = clock();
		return 0; // reussite
	}
	std::cout << "Failed to place individual" << std::endl;
	throw std::string("Failed to place individual");
	return 1000000;
}

void Population::removeIndividu(SousPop *pop, unsigned int p) {
	Individual *partant = pop->individus[p];

	// on place notre individu a la fin
	for (unsigned int i = p + 1; i < pop->individus.size(); i++) {
		pop->individus[i - 1] = pop->individus[i];
	}

	// on l'enleve de la population
	pop->individus.pop_back();
	pop->nbIndiv--;

	// on l'enleve des structures de proximit�
	for (unsigned int i = 0; i < pop->nbIndiv; i++)
		pop->individus[i]->removeNearest(partant);

	// et on supprime le partant
	delete partant;
}
// recalcule l'evaluation des individus a partir des violation
// puis effectue un tri a bulles de la population
// operateur de tri peu efficace mais fonction appelee tres rarement
void Population::validatePen() {
	// on met a jour les evaluations
	for (Individual* indiv : invalid->individus) {
		indiv->solutionCost.evaluation = 0.0;
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			indiv->solutionCost_scenario.evaluation[scenario] = indiv->solutionCost_scenario.fitness[scenario] + params->penalityCapa[scenario] * indiv->solutionCost_scenario.capacityViol[scenario]; 
			indiv->solutionCost.evaluation += indiv->solutionCost_scenario.evaluation[scenario];
		}
		indiv->solutionCost.evaluation /= (double)nbScenario;
	}
	for (unsigned int i = 0; i < invalid->nbIndiv; i++) {

		for (unsigned int j = 0; j < invalid->nbIndiv - i - 1; j++) {
			if (invalid->individus[j]->solutionCost.evaluation > invalid->individus[j + 1]->solutionCost.evaluation) {
				std::swap(invalid->individus[j], invalid->individus[j + 1]);
			}
		}
	}
}

Individual *Population::getIndividualBinT() {
	Individual *individu1;
	Individual *individu2;
	unsigned int place1, place2;

	place1 = (unsigned int) (params->rng->genrand64_int64() % (valid->nbIndiv + invalid->nbIndiv));
	if (place1 >= valid->nbIndiv) {
		place1 -= valid->nbIndiv;
		individu1 = invalid->individus[place1];
	} else {
		individu1 = valid->individus[place1];
	}

	place2 = (unsigned int) (params->rng->genrand64_int64() % (valid->nbIndiv + invalid->nbIndiv));
	if (place2 >= valid->nbIndiv) {
		place2 -= valid->nbIndiv;
		individu2 = invalid->individus[place2];
	} else {
		individu2 = valid->individus[place2];
	}

	evalExtFit(valid);
	evalExtFit(invalid);

	if (individu1->extendedFitness < individu2->extendedFitness) {
		return individu1;
	} else {
		return individu2;
	}
}

Individual *Population::getIndividuBestValide() {
	if (valid->nbIndiv != 0) return valid->individus[0];
	return NULL;
}

Individual *Population::getIndividuBestInvalide() {
	if (invalid->nbIndiv != 0) return invalid->individus[0];
	return NULL;
}

// getter simple d'un individu
Individual *Population::getIndividu(unsigned int p) {
	return valid->individus[p];
}
// recopie un Individu dans un autre
// ATTENTION !!! ne recopie que le chromT et les attributs du fitness
void Population::copyIndividual(Individual *destination, Individual *source) {
	destination->chromT = source->chromT;
	destination->chromL = source->chromL;
	destination->solutionCost = source->solutionCost;
	destination->solutionCost_scenario = source->solutionCost_scenario;
	destination->isFeasible = source->isFeasible;
}

void Population::ExportPop(string nomFichier,bool add, std::vector<double> deliveries, double &totalCost) {
	// exporte les solutions actuelles des individus dans un dossier exports current individual solutions to a folder
	int compteur;
	Individual *bestValide = getIndividuBestValide();
	
	if (bestValide != NULL) {
		// TO CHECK
		ofstream myfile;
		clock_t s = totalTime;
		clock_t v = timeBest;
		
		// We will update the local search structure for paths.
		// We are obliged to set very strong parameters so that the splitting does not produce a from the best valid solution
		// so that the splitting does not produce a from the best valid solution
		std::vector<double> temp = params->penalityCapa;
		params->penalityCapa = std::vector<double> (params->nbScenarios, 10000.0);
		education_scenario(bestValide);
		// le trainer a garde les infos des routes de bestValide
		LocalSearch *loc = trainer->localSearchList[0];
		params->penalityCapa = temp;
		
		myfile.precision(10);
		cout.precision(10);
		bool ADD = add && (params->jVal != 1);
		if (ADD) myfile.open(nomFichier.data(), std::ios::app);//add on previous
		else myfile.open(nomFichier.data()); 
		if (params->jVal == 1) myfile << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
		myfile << "Day "<< params->jVal << endl;
		myfile << endl;
		loc->printInventoryLevels(myfile, deliveries, totalCost);
		// export cost
		myfile << "Total cost: " << totalCost << endl;

		// exporting the number of routes
		compteur = 0;
		for (unsigned int i = 0; i < loc->routes[1].size(); i++)
			if (!loc->depots[1][i]->next->isADepot)
				compteur++;
		myfile << "Route number: " << compteur << endl;
		myfile << endl;

		// exporting the total CPU time (ms)
		myfile <<"Total Time: " << (float) s / (float)CLOCKS_PER_SEC <<endl; 

		// exporting the time to best solution
		myfile <<"Best Solution Time: " << (float) v / (float)CLOCKS_PER_SEC <<endl; 

		myfile << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
		myfile.close();
		if (params->traces) std::cout << "Successful export" << std::endl;
	}
	else
	{
		cout << " impossible to find a valid individual " << endl;
	}
}

// retourne la fraction d'individus valides en terme de charge
double Population::validChargePart() {
	int count = 0;
	for (list<bool>::iterator it = listeValiditeCharge.begin(); it != listeValiditeCharge.end(); ++it) {
		count += *it;
	}
	return double(count) / (double)(100);
}

double Population::getDiversity(SousPop *pop) {
	double total = 0;
	int count = 0;
	for (unsigned int i = 0; i < pop->nbIndiv / 2; i++) {
		for (unsigned int j = i + 1; j < pop->nbIndiv; j++){
			total += pop->individus[i]->distance(pop->individus[j]);
			count++;
		}
	}
	return total / (double)count;
}

double Population::getMoyenneValides() {
	double moyenne = 0;
	for (Individual* indiv : valid->individus)
		moyenne += indiv->solutionCost.evaluation;
	return moyenne / (double)valid->nbIndiv;
}
double Population::getMoyenneInvalides() {
	double moyenne = 0;
	for (Individual* indiv : invalid->individus) {
		moyenne += indiv->solutionCost.evaluation;
	}
	return moyenne / (double)invalid->nbIndiv;
}

unsigned int Population::selectCompromis(SousPop *souspop) {
	vector<unsigned int> classement;

	evalExtFit(souspop);

	// pour chaque individu on modifie le fitness etendu
	for (unsigned int i = 0; i < souspop->nbIndiv; i++) {
		classement.push_back(i);
		if (souspop->individus[i]->distNearest(1) < params->distMin) souspop->individus[i]->extendedFitness += 5.0;
		// for the CVRP instances, we need to allow duplicates with the same fitness since in the Golden instances
		// there is a lot of symmetry.
		if (fitExist(souspop, souspop->individus[i])) souspop->individus[i]->extendedFitness += 5.0;
	}

	// on classe les elements par fitness etendu et on prend le plus mauvais
	for (unsigned int n = 0; n < souspop->nbIndiv; n++) {
		for (unsigned int i = 1; i < souspop->nbIndiv - n - 1; i++) {
			if (souspop->individus[classement[i]]->extendedFitness > souspop->individus[classement[i + 1]]->extendedFitness) {
				std::swap(classement[i], classement[i+1]);
			}
		}
	}
	return classement[souspop->nbIndiv - 1];
}

void Population::education_scenario(Individual *indiv) {
	copyIndividual(trainer, indiv); //copie provisoire de indiv dans trainer
	trainer->split(); //split du grand tour en différentes routes + mesure coût
	trainer->updateLocalSearch(); //on remplit les structures de recherches locales grâce à chromL et chromT
	trainer->runLocalSearch(); //phase de recherche locale
	trainer->updateIndividual();  //on remplit chromT et chromL avec les résultats de LS
	copyIndividual(indiv, trainer); //on recopie dans indiv le trainer
}


void Population::updateNbValid(Individual *indiv)
{
	listeValiditeCharge.push_back(indiv->solutionCost.capacityViol < 0.0001);
	listeValiditeCharge.pop_front();
}

void Population::measureAndUpdateQuantities(std::vector<double> &deliveries, double &totalCost) {
	if (getIndividuBestValide() != NULL) {
      Individual * bestIndividual = getIndividuBestValide();
      bestIndividual->split();
      double val = bestIndividual->measureTrueCost(deliveries);
      if (params->traces) std::cout << "Cost of the day " << params->jVal << ": " << val << std::endl;
      totalCost += val;
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

	if (getIndividuBestValide() != NULL) {
		cout << getIndividuBestValide()->solutionCost.evaluation << " Scenarios - ";
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			cout << scenario + 1 << ": " << getIndividuBestValide()->solutionCost_scenario.evaluation[scenario] << " ";
		}
	} else {
		cout << "NO-VALID ";
	}

	cout << " | ";

	if (getIndividuBestInvalide() != NULL) {
		cout << getIndividuBestInvalide()->solutionCost.evaluation << " Scenarios - ";
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			cout << scenario << ": " << getIndividuBestInvalide()->solutionCost_scenario.evaluation[scenario] << " ";
		}
	} else {
		cout << "NO-INVALID";
	}

	cout.precision(8);
	double avgPenalityCapa = 0.0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		avgPenalityCapa += params->penalityCapa[scenario];
	}
	avgPenalityCapa /= (double) nbScenario;

	cout << " | Moy " << getMoyenneValides() << " " << getMoyenneInvalides()
		 << " | Div " << getDiversity(valid) << " " << getDiversity(invalid)
		 << " | Val " << validChargePart()
		 << " | Pen moy " << avgPenalityCapa << " | Pop " << valid->nbIndiv << " " << invalid->nbIndiv << endl;
}
