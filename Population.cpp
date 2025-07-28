#include "Population.h"
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
// constructeur

Population::Population(vector<Params*> pl) : paramsList(pl)
{
	nbScenario = (int) paramsList.size();

	trainer = new Individu(paramsList);
	valides = new SousPop(); 
	invalides = new SousPop();

	valides->nbIndiv = 0;
	invalides->nbIndiv = 0;
	listeValiditeCharge = list<bool> (100, true);
	bool compter = true;
	
	vector<double> saveCapa(nbScenario);
	
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++)
	saveCapa[scenario] = paramsList[scenario]->penalityCapa;
	
	for (int i = 0; i < paramsList[0]->mu * 2; i++) {
		if (i == paramsList[0]->mu) {
			for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
				paramsList[scenario]->penalityCapa *= 50;
			}
			compter = false;
		}
		Individu *randomIndiv = new Individu(paramsList);  
		
		education_scenario(randomIndiv);
		if (compter) updateNbValides(randomIndiv);
		addIndividu(randomIndiv);
		delete randomIndiv;
	}
	
	// on initialise par defaut a 100, comme si tout etait valide au debut
	// mais c'est arbitraire
	
	//TO CHECK
	for (unsigned int scenario = 0; scenario < (unsigned int) nbScenario; scenario++) {
		paramsList[scenario]->penalityCapa = saveCapa[scenario];
	}
}

// destructeur
Population::~Population()
{
	for (int i = 0; i < (int)valides->individus.size(); i++)
		delete valides->individus[i];

	for (int i = 0; i < (int)invalides->individus.size(); i++)
		delete invalides->individus[i];

	delete trainer;
}

void Population::evalExtFit(SousPop *pop) {
	int temp;
	vector<int> classement;
	vector<double> distances;

	for (unsigned int indiv = 0; indiv < pop->nbIndiv; indiv++) {
		classement.push_back(indiv);
		distances.push_back(pop->individus[indiv]->distPlusProche(paramsList[0]->nbCountDistMeasure));
	}

	// classement des individus en fonction de leur note de distance
	for (int n = 0; n < pop->nbIndiv; n++) {
		for (int i = 0; i < pop->nbIndiv - n - 1; i++) {
			if (distances[classement[i]] < distances[classement[i + 1]] - 0.000001)	{
				temp = classement[i + 1];
				classement[i + 1] = classement[i];
				classement[i] = temp;
			}
		}
	}

	for (unsigned int i = 0; i < pop->nbIndiv; i++) {
		pop->individus[classement[i]]->divRank = (float)i / (float)(pop->nbIndiv - 1);
		pop->individus[classement[i]]->fitRank = (float)classement[i] / (float)(pop->nbIndiv - 1);
		pop->individus[classement[i]]->fitnessEtendu = pop->individus[classement[i]]->fitRank + ((float)1.0 - (float)paramsList[0]->el / (float)pop->nbIndiv) * pop->individus[classement[i]]->divRank;
	}
}

int Population::addIndividu(Individu *indiv) {
	SousPop *souspop = (indiv->estValide) ? valides : invalides;
	int result = placeIndividu(souspop, indiv);

	// il faut eventuellement enlever la moitie de la pop
	if (souspop->nbIndiv > paramsList[0]->mu + paramsList[0]->lambda) {
		while (souspop->nbIndiv > paramsList[0]->mu) {
			removeIndividu(souspop, selectCompromis(souspop));
		}
	}
	return result;
}

// met a jour les individus les plus proches d'une population
// en fonction de l'arrivant

void Population::updateProximity(SousPop *pop, Individu *indiv) {
	for (unsigned int k = 0; k < pop->nbIndiv; k++)
		if (pop->individus[k] != indiv) {
			pop->individus[k]->addProche(indiv);
			indiv->addProche(pop->individus[k]);
		}
}

// fonction booleenne verifiant si le fitness n'existe pas d�ja

bool Population::fitExist(SousPop *pop, Individu *indiv)
{
	double fitness = indiv->coutSol.evaluation;
	for (int i = 0; i < (int)pop->nbIndiv; i++) {
		if (indiv != pop->individus[i] && pop->individus[i]->coutSol.evaluation >= (fitness - paramsList[0]->delta) && pop->individus[i]->coutSol.evaluation <= (fitness + paramsList[0]->delta))
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
		savePenalities[scenario] = paramsList[scenario]->penalityCapa;
	}

	while (valides->nbIndiv > (int)(paramsList[0]->rho * (double)paramsList[0]->mu)) {
		delete valides->individus[valides->nbIndiv - 1];
		valides->individus.pop_back();
		valides->nbIndiv--;
	}
	while (invalides->nbIndiv > (int)(paramsList[0]->rho * (double)paramsList[0]->mu)) {
		delete invalides->individus[invalides->nbIndiv - 1];
		invalides->individus.pop_back();
		invalides->nbIndiv--;
	}
	for (int i = 0; i < paramsList[0]->mu * 2; i++) {
		if (i == paramsList[0]->mu) {
			for (int j = 0; j < nbScenario; j++) {
				paramsList[j]->penalityCapa *= 50;
			}
		}
		Individu *random = new Individu(paramsList);
		education_scenario(random);
		addIndividu(random);
	}
	for (unsigned int scenario = 0; scenario < (unsigned int) nbScenario; scenario++) {
		paramsList[scenario]->penalityCapa = savePenalities[scenario];
	}
}

//Try to place individual indiv in the sub population and return its position 
//Individuals are increasingly sorted in the subpopulation
int Population::placeIndividu(SousPop *pop, Individu *indiv) {
	Individu *monIndiv = new Individu(paramsList);
	recopieIndividu(monIndiv, indiv);

	// regarde si son fitness est suffisamment espace
	bool placed = false;
	int i = (int)pop->individus.size() - 1;

	pop->individus.push_back(monIndiv);
	while (i >= 0 && !placed) {
		if (pop->individus[i]->coutSol.evaluation >= indiv->coutSol.evaluation + 0.001) {
			pop->individus[i + 1] = pop->individus[i];
			i--;
		} else {
			pop->individus[i + 1] = monIndiv;
			placed = true;
			pop->nbIndiv++;
			updateProximity(pop, monIndiv);
			return i + 1; // reussite
		}
	}
	if (!placed)
	{
		pop->individus[0] = monIndiv;
		placed = true;
		pop->nbIndiv++;
		updateProximity(pop, monIndiv);
		if (pop == valides) timeBest = clock();
		return 0; // reussite
	}
	throw std::string("failed to place individual");
	return -3;
}

void Population::removeIndividu(SousPop *pop, int p) {
	Individu *partant = pop->individus[p];

	// on place notre individu a la fin
	for (unsigned int i = p + 1; i < pop->individus.size(); i++) {
		pop->individus[i - 1] = pop->individus[i];
	}

	// on l'enleve de la population
	pop->individus.pop_back();
	pop->nbIndiv--;

	// on l'enleve des structures de proximit�
	for (unsigned int i = 0; i < pop->nbIndiv; i++)
		pop->individus[i]->removeProche(partant);

	// et on supprime le partant
	delete partant;
}
// recalcule l'evaluation des individus a partir des violation
// puis effectue un tri a bulles de la population
// operateur de tri peu efficace mais fonction appelee tres rarement
void Population::validatePen() {
	// on met a jour les evaluations
	for (Individu* indiv : invalides->individus) {
		indiv->coutSol.evaluation = 0.0;
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			indiv->coutSol_scenario.evaluation[scenario] = indiv->coutSol_scenario.fitness[scenario] + paramsList[scenario]->penalityCapa * indiv->coutSol_scenario.capacityViol[scenario]; 
			indiv->coutSol.evaluation += indiv->coutSol_scenario.evaluation[scenario];
		}
		indiv->coutSol.evaluation /= (double)nbScenario;
	}
	for (unsigned int i = 0; i < invalides->nbIndiv; i++) {

		for (unsigned int j = 0; j < invalides->nbIndiv - i - 1; j++) {
			if (invalides->individus[j]->coutSol.evaluation > invalides->individus[j + 1]->coutSol.evaluation) {
				std::swap(invalides->individus[j], invalides->individus[j + 1]);
			}
		}
	}
}

Individu *Population::getIndividuBinT(double &rangRelatif) {
	Individu *individu1;
	Individu *individu2;
	int place1, place2;
	double rangRelatif1, rangRelatif2;

	place1 = paramsList[0]->rng->genrand64_int64() % (valides->nbIndiv + invalides->nbIndiv);
	if (place1 >= valides->nbIndiv) {
		place1 -= valides->nbIndiv;
		individu1 = invalides->individus[place1];
		rangRelatif1 = (double) place1 / (double)invalides->nbIndiv;
	} else {
		individu1 = valides->individus[place1];
		rangRelatif1 = (double) place1 / (double)valides->nbIndiv;
	}

	place2 = paramsList[0]->rng->genrand64_int64() % (valides->nbIndiv + invalides->nbIndiv);
	if (place2 >= valides->nbIndiv) {
		place2 -= valides->nbIndiv;
		individu2 = invalides->individus[place2];
		rangRelatif2 = (double) place2 / (double)invalides->nbIndiv;
	} else {
		individu2 = valides->individus[place2];
		rangRelatif2 = (double) place2 / (double)valides->nbIndiv;
	}

	evalExtFit(valides);
	evalExtFit(invalides);

	if (individu1->fitnessEtendu < individu2->fitnessEtendu) {
		rangRelatif = rangRelatif1;
		return individu1;
	} else {
		rangRelatif = rangRelatif2;
		return individu2;
	}
}

Individu *Population::getIndividuBestValide() {
	if (valides->nbIndiv != 0) return valides->individus[0];
	return NULL;
}

Individu *Population::getIndividuBestInvalide() {
	if (invalides->nbIndiv != 0) return invalides->individus[0];
	return NULL;
}

// getter simple d'un individu
Individu *Population::getIndividu(int p) {
	return valides->individus[p];
}
// recopie un Individu dans un autre
// ATTENTION !!! ne recopie que le chromT et les attributs du fitness
void Population::recopieIndividu(Individu *destination, Individu *source) {
	destination->chromT = source->chromT;
	destination->chromL = source->chromL;
	destination->coutSol = source->coutSol;
	destination->coutSol_scenario = source->coutSol_scenario;
	destination->estValide = source->estValide;
}

void Population::ExportPop(string nomFichier,bool add) {
	// exporte les solutions actuelles des individus dans un dossier exports current individual solutions to a folder
	vector<int> rout;
	vector<double> routTime;
	int compteur;
	Noeud *noeudActuel;
	LocalSearch *loc;
	ofstream myfile;
	double cost;
	double temp, temp2;
	char *myBuff;
	Individu *bestValide = getIndividuBestValide();

	if (bestValide != NULL) {

		// We will update the local search structure for paths.
		// We are obliged to set very strong parameters so that the splitting does not produce a from the best valid solution
		// so that the splitting does not produce a from the best valid solution
		temp = paramsList[0]->penalityCapa;
		paramsList[0]->penalityCapa = 10000;
		education_scenario(bestValide);
		// le trainer a gard� les infos des routes de bestValide
		loc = trainer->localSearchList[0];
		paramsList[0]->penalityCapa = temp;

		myfile.precision(10);
		cout.precision(10);
		ofstream myfile;
		if (add) myfile.open(nomFichier.data(), std::ios::app);//add on previous
		else myfile.open(nomFichier.data()); 
		myfile<<endl<<endl;
		loc->printInventoryLevels(myfile,add);
		// export cost
		myfile << trainer->coutSol.evaluation << endl;

		// exporting the number of routes
		compteur = 0;
		for (int k = 1; k <= paramsList[0]->nbDays; k++)
			for (int i = 0; i < (int)loc->routes[k].size(); i++)
				if (!loc->depots[k][i]->suiv->estUnDepot)
					compteur++;
		myfile << compteur << endl;

		// exporting the total CPU time (ms)
		myBuff = new char[100];
		myfile <<"Total Time: ";sprintf(myBuff, "%d", (int)(clock() / 1000000));
		myfile << myBuff << endl;

		myBuff = new char[100];
		myfile <<"PITime: ";sprintf(myBuff, "%d", (int)(paramsList[0]->debut / 1000000));
		myfile << myBuff << endl;

		// exporting the time to best solution
		myBuff = new char[100];
		myfile <<"Best Solution Time: ";sprintf(myBuff, "%d", (int)(timeBest / 1000000));
		myfile << myBuff << endl;
		for (int k = 1; k <= paramsList[0]->nbDays; k++)
		{
			for (int i = 0; i < (int)loc->routes[k].size(); i++)
			{
				compteur = 1;
				if (!loc->depots[k][i]->suiv->estUnDepot)
				{

					myfile << " " << loc->routes[k][i]->depot->idx << " " << (k - 1) % paramsList[0]->ancienNbDays + 1 << " " << compteur << " " << loc->routes[k][i]->temps << " "
						   << loc->routes[k][i]->charge << " ";

					// on place la route dans un vecteur de noeuds clients
					noeudActuel = loc->routes[k][i]->depot->suiv;
					cost = paramsList[0]->timeCost[loc->routes[k][i]->depot->idx][noeudActuel->idx];

					rout.clear();
					routTime.clear();
					rout.push_back(loc->routes[k][i]->depot->idx);
					routTime.push_back(0);
					rout.push_back(noeudActuel->idx);
					routTime.push_back(cost);

					while (!noeudActuel->estUnDepot)
					{
						cost += paramsList[0]->timeCost[noeudActuel->idx][noeudActuel->suiv->idx];
						noeudActuel = noeudActuel->suiv;
						rout.push_back(noeudActuel->idx);
						routTime.push_back(cost);
					}

					myfile << " " << (int)rout.size();

					for (int j = 0; j < (int)rout.size(); j++)
					{
						myfile << " " << rout[j];
					}
					myfile << endl;
					compteur++;
				}
			}
		}

		myfile.close();
		std::cout << "Successful export" << std::endl;
	}
	else
	{
		cout << " impossible to find a valid individual " << endl;
	}
}

void Population::ExportBKS(string nomFichier)
{
	double fit,tim,pri;
	ifstream fichier;
    std::string line;
	fichier.open(nomFichier.c_str());
	if (fichier.is_open())
	{
		while (getline(fichier, line)) {
            std::size_t found = line.find("COST SUMMARY : OVERALL");
            if (found != std::string::npos) {
                std::stringstream ss(line.substr(found));
                std::string temp;
			
                ss >> temp >> temp>> temp >> temp>>fit;
               
            }
			found = line.find("Total Time:");
            if (found != std::string::npos) {
                std::stringstream ss(line.substr(found));
                std::string temp;
			
                ss >> temp >> temp>>tim;
                break; 
            }
        }
		fichier.close();
		timeBest = clock();
		if (getIndividuBestValide() != NULL && getIndividuBestValide()->coutSol.evaluation < fit - 0.01)
		{
			cout << "!!! new BKS !!! : " << getIndividuBestValide()->coutSol.evaluation << endl;
			ExportPop(nomFichier,false);
		}
		else if (getIndividuBestValide() != NULL && std::fabs(getIndividuBestValide()->coutSol.evaluation - fit) < 0.01 &&  tim > (int)(timeBest / 1000000))
		{
			cout << "!!! new time !!! : " << getIndividuBestValide()->coutSol.evaluation << endl;
			
			ExportPop(nomFichier,false);
		}
		
	}
	else
	{
		cout << " BKS file not found " << endl;
		ExportPop(nomFichier,false);
	}
}
// retourne la fraction d'individus valides en terme de charge
double Population::fractionValidesCharge() {
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
	for (Individu* indiv : valides->individus)
		moyenne += indiv->coutSol.evaluation;
	return moyenne / (double)valides->nbIndiv;
}
double Population::getMoyenneInvalides() {
	double moyenne = 0;
	for (Individu* indiv : invalides->individus) {
		moyenne += indiv->coutSol.evaluation;
	}
	return moyenne / (double)invalides->nbIndiv;
}

int Population::selectCompromis(SousPop *souspop) {
	vector<int> classement;

	evalExtFit(souspop);

	// pour chaque individu on modifie le fitness etendu
	for (unsigned int i = 0; i < (unsigned int) souspop->nbIndiv; i++) {
		classement.push_back(i);
		if (souspop->individus[i]->distPlusProche(1) < paramsList[0]->distMin) souspop->individus[i]->fitnessEtendu += 5;
		// for the CVRP instances, we need to allow duplicates with the same fitness since in the Golden instances
		// there is a lot of symmetry.
		if (fitExist(souspop, souspop->individus[i])) souspop->individus[i]->fitnessEtendu += 5.0;
	}

	// on classe les elements par fitness etendu et on prend le plus mauvais
	for (unsigned int n = 0; n < souspop->nbIndiv; n++) {
		for (unsigned int i = 1; i < souspop->nbIndiv - n - 1; i++) {
			if (souspop->individus[classement[i]]->fitnessEtendu > souspop->individus[classement[i + 1]]->fitnessEtendu) {
				std::swap(classement[i], classement[i+1]);
			}
		}
	}
	return classement[souspop->nbIndiv - 1];
}

void Population::education_scenario(Individu *indiv) {
	recopieIndividu(trainer, indiv); //copie provisoire de indiv dans trainer
	trainer->generalSplit_scenario(); //split du grand tour en différentes routes + mesure coût
	trainer->updateLS_scenario(); //on remplit les structures de recherches locales grâce à chromL et chromT
	trainer->localSearchRunSearch_scenario(); //phase de recherche locale
	trainer->updateIndiv_scenario();  //on remplit chromT et chromL avec les résultats de LS
	recopieIndividu(indiv, trainer); //on recopie dans indiv le trainer
}


// met a jour le compte des valides
void Population::updateNbValides(Individu *indiv)
{
	listeValiditeCharge.push_back(indiv->coutSol.capacityViol < 0.0001);
	listeValiditeCharge.pop_front();
}

void Population::afficheEtat(int nbIter)
{
	cout.precision(8);

	cout << "It " << nbIter << " | Sol moy: ";

	if (getIndividuBestValide() != NULL) {
		cout << getIndividuBestValide()->coutSol.evaluation << " ";
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			cout << scenario << ": " << getIndividuBestValide()->coutSol_scenario.evaluation[scenario] << " ";
		}
	} else {
		cout << "NO-VALID ";
	}

	if (getIndividuBestInvalide() != NULL) {
		cout << getIndividuBestInvalide()->coutSol.evaluation << " ";
		for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
			cout << scenario << ": " << getIndividuBestInvalide()->coutSol_scenario.evaluation[scenario] << " ";
		}
	} else {
		cout << "NO-INVALID";
	}

	cout.precision(8);
	double avgPenalityCapa = 0.0;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
		avgPenalityCapa += paramsList[scenario]->penalityCapa;
	}
	avgPenalityCapa /= (double) nbScenario;

	cout << " | Moy " << getMoyenneValides() << " " << getMoyenneInvalides()
		 << " | Div " << getDiversity(valides) << " " << getDiversity(invalides)
		 << " | Val " << fractionValidesCharge()
		 << " | Pen moy " << avgPenalityCapa << " | Pop " << valides->nbIndiv << " " << invalides->nbIndiv << endl;
}
