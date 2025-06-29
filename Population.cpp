#include "Population.h"
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
// constructeur

Population::Population(vector<Params*> pl) : paramsList(pl)
{

	nbScenario = paramsList.size();
	params = paramsList[0];

	trainer = new Individu(paramsList);

	double temp = params->penalityCapa;
	
	valides = new SousPop(); 
	invalides = new SousPop();

	valides->nbIndiv = 0;
	invalides->nbIndiv = 0;
	valides->nbGenerations = 0;
	invalides->nbGenerations = 0;
	bool compter = true;
	
	vector<double> saveCapa;
	vector<double> saveLength;

	for (int scenario = 0; scenario < nbScenario; scenario++) {
		saveCapa.push_back(paramsList[scenario]->penalityCapa);
	}

	for (int i = 0; i < paramsList[0]->mu * 2; i++)
	{
		if (i == paramsList[0]->mu)
		{
			paramsList[0]->penalityCapa *= 50;
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
	listeValiditeCharge = list<bool> (100, true);
	for (int scenario = 0; scenario < nbScenario; scenario++) {
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
void Population::evalExtFit(SousPop *pop)
{
	int temp;
	vector<int> classement;
	vector<double> distances;

	for (int i = 0; i < pop->nbIndiv; i++)
	{
		classement.push_back(i);
		distances.push_back(pop->individus[i]->distPlusProche(params->nbCountDistMeasure));
	}

	// classement des individus en fonction de leur note de distance
	for (int n = 0; n < pop->nbIndiv; n++)
		for (int i = 0; i < pop->nbIndiv - n - 1; i++)
			if (distances[classement[i]] < distances[classement[i + 1]] - 0.000001)
			{
				temp = classement[i + 1];
				classement[i + 1] = classement[i];
				classement[i] = temp;
			}

	for (int i = 0; i < pop->nbIndiv; i++)
	{
		pop->individus[classement[i]]->divRank = (float)i / (float)(pop->nbIndiv - 1);
		pop->individus[classement[i]]->fitRank = (float)classement[i] / (float)(pop->nbIndiv - 1);
		pop->individus[classement[i]]->fitnessEtendu = pop->individus[classement[i]]->fitRank + ((float)1.0 - (float)params->el / (float)pop->nbIndiv) * pop->individus[classement[i]]->divRank;
	}
}

int Population::addIndividu(Individu *indiv)
{
	SousPop *souspop = (indiv->estValide) ? valides : invalides;
	int k, result;

	result = placeIndividu(souspop, indiv);
	// il faut eventuellement enlever la moitie de la pop
	if (result != -1 && souspop->nbIndiv > params->mu + params->lambda)
	{

		while (souspop->nbIndiv > params->mu)
		{
			k = selectCompromis(souspop);
			removeIndividu(souspop, k);
		}

		souspop->nbGenerations++;
	}
	return result;
}

// met a jour les individus les plus proches d'une population
// en fonction de l'arrivant

void Population::updateProximity(SousPop *pop, Individu *indiv)
{
	for (int k = 0; k < pop->nbIndiv; k++)
		if (pop->individus[k] != indiv)
		{
			pop->individus[k]->addProche(indiv);
			indiv->addProche(pop->individus[k]);
		}
}

// fonction booleenne verifiant si le fitness n'existe pas d�ja

bool Population::fitExist(SousPop *pop, Individu *indiv)
{
	double fitness = indiv->coutSol.evaluation;
	for (int i = 0; i < (int)pop->nbIndiv; i++)
	{
		if (indiv != pop->individus[i] && pop->individus[i]->coutSol.evaluation >= (fitness - params->delta) && pop->individus[i]->coutSol.evaluation <= (fitness + params->delta))
			return true;
	}
	return false;
}
// procede de redemarrage avec remplacement d'une partie de la population
// modele tres simplifie
// on remplace la moitie individus de fitness situes sous la moyenne par des individus aleatoires
void Population::diversify()
{
	Individu *randomIndiv;
	double temp = params->penalityCapa;

	while (valides->nbIndiv > (int)(params->rho * (double)params->mu))
	{
		delete valides->individus[valides->nbIndiv - 1];
		valides->individus.pop_back();
		valides->nbIndiv--;
	}
	while (invalides->nbIndiv > (int)(params->rho * (double)params->mu))
	{
		delete invalides->individus[invalides->nbIndiv - 1];
		invalides->individus.pop_back();
		invalides->nbIndiv--;
	}

	for (int i = 0; i < params->mu * 2; i++)
	{
		if (i == params->mu)
		{
			params->penalityCapa *= 50;
		}
		randomIndiv = new Individu(paramsList);
		education(randomIndiv);
		addIndividu(randomIndiv);
		delete randomIndiv;
	}

	params->penalityCapa = temp;
}

int Population::placeIndividu(SousPop *pop, Individu *indiv)
{

	Individu *monIndiv = new Individu(paramsList);
	recopieIndividu(monIndiv, indiv);

	// regarde si son fitness est suffisamment espace
	bool placed = false;
	int i = (int)pop->individus.size() - 1;

	pop->individus.push_back(monIndiv);
	while (i >= 0 && !placed)
	{
		if (pop->individus[i]->coutSol.evaluation >= indiv->coutSol.evaluation + 0.001)
		{
			pop->individus[i + 1] = pop->individus[i];
			i--;
		}
		else
		{
			pop->individus[i + 1] = monIndiv;
			placed = true;
			pop->nbIndiv++;
			updateProximity(pop, pop->individus[i + 1]);
			return i + 1; // reussite
		}
	}
	if (!placed)
	{
		pop->individus[0] = monIndiv;
		placed = true;
		pop->nbIndiv++;
		updateProximity(pop, pop->individus[0]);
		if (pop == valides)
			timeBest = clock();
		return 0; // reussite
	}
	cout << "erreur placeIndividu" << endl;
	return -3;
}

void Population::removeIndividu(SousPop *pop, int p)
{
	Individu *partant = pop->individus[p];

	// on place notre individu � la fin
	for (int i = p + 1; i < (int)pop->individus.size(); i++)
		pop->individus[i - 1] = pop->individus[i];

	// on l'enleve de la population
	pop->individus.pop_back();
	pop->nbIndiv--;

	// on l'enleve des structures de proximit�
	for (int i = 0; i < pop->nbIndiv; i++)
		pop->individus[i]->removeProche(partant);

	// et on supprime le partant
	delete partant;
}
// recalcule l'evaluation des individus a partir des violation
// puis effectue un tri a bulles de la population
// operateur de tri peu efficace mais fonction appelee tres rarement
void Population::validatePen()
{
	Individu *indiv;
	// on met � jour les evaluations
	for (int i = 0; i < invalides->nbIndiv; i++)
		invalides->individus[i]->coutSol.evaluation = invalides->individus[i]->coutSol.fitness + params->penalityCapa * invalides->individus[i]->coutSol.capacityViol;

	for (int i = 0; i < invalides->nbIndiv; i++)
		for (int j = 0; j < invalides->nbIndiv - i - 1; j++)
		{
			if (invalides->individus[j]->coutSol.evaluation > invalides->individus[j + 1]->coutSol.evaluation)
			{
				indiv = invalides->individus[j];
				invalides->individus[j] = invalides->individus[j + 1];
				invalides->individus[j + 1] = indiv;
			}
		}
}

Individu *Population::getIndividuBinT(double &rangRelatif)
{
	Individu *individu1;
	Individu *individu2;
	int place1, place2;
	double rangRelatif1, rangRelatif2;

	place1 = params->rng->genrand64_int64() % (valides->nbIndiv + invalides->nbIndiv);
	if (place1 >= valides->nbIndiv)
	{
		place1 -= valides->nbIndiv;
		individu1 = invalides->individus[place1];
		rangRelatif1 = (double) place1 / (double)invalides->nbIndiv;
	} else
	{
		individu1 = valides->individus[place1];
		rangRelatif1 = (double) place1 / (double)valides->nbIndiv;
	}

	place2 = params->rng->genrand64_int64() % (valides->nbIndiv + invalides->nbIndiv);
	if (place2 >= valides->nbIndiv)
	{
		place2 -= valides->nbIndiv;
		individu2 = invalides->individus[place2];
		rangRelatif2 = (double) place2 / (double)invalides->nbIndiv;
	} else
	{
		individu2 = valides->individus[place2];
		rangRelatif2 = (double) place2 / (double)valides->nbIndiv;
	}

	evalExtFit(valides);
	evalExtFit(invalides);

	if (individu1->fitnessEtendu < individu2->fitnessEtendu)
	{
		rangRelatif = rangRelatif1;
		return individu1;
	}
	else
	{
		rangRelatif = rangRelatif2;
		return individu2;
	}
}

// rank donne le rang de l'individu choisi
Individu *Population::getIndividuRand(double &rangRelatif)
{
	Individu *individu1;

	int place1 = params->rng->genrand64_int64() % (valides->nbIndiv + invalides->nbIndiv);
	if (place1 >= valides->nbIndiv)
	{
		individu1 = invalides->individus[place1 - valides->nbIndiv];
		rangRelatif = (double)(place1 - valides->nbIndiv) / (double)invalides->nbIndiv;
	}
	else
	{
		individu1 = valides->individus[place1];
		rangRelatif = (double)place1 / (double)valides->nbIndiv;
	}
	return valides->individus[place1];
}
// getter de 1 individu par selection au hasard dans un certain % des meilleurs meilleurs
// rank donne le rang de l'individu choisi
Individu *Population::getIndividuPourc(int pourcentage, double &rangRelatif)
{
	int place;
	if ((valides->nbIndiv * pourcentage) / 100 != 0)
	{
		place = params->rng->genrand64_int64() % ((valides->nbIndiv * pourcentage) / 100);
		rangRelatif = (double)place / (double)params->mu;
		return valides->individus[place];
	}
	else
		return getIndividuBinT(rangRelatif);
}
Individu *Population::getIndividuBestValide()
{
	if (valides->nbIndiv != 0) return valides->individus[0];
	else return NULL;
}

Individu *Population::getIndividuBestInvalide()
{
	if (invalides->nbIndiv != 0) return invalides->individus[0];
	else return NULL;
}

// getter simple d'un individu
Individu *Population::getIndividu(int p)
{
	return valides->individus[p];
}
// recopie un Individu dans un autre
// ATTENTION !!! ne recopie que le chromT et les attributs du fitness
void Population::recopieIndividu(Individu *destination, Individu *source)
{
	destination->chromT = source->chromT;
	destination->chromL = source->chromL;
	destination->coutSol = source->coutSol;
	destination->estValide = source->estValide;
}

void Population::ExportPop(string nomFichier,bool add)
{
	
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

	if (bestValide != NULL)
	{

		// We will update the local search structure for paths.
		// We are obliged to set very strong parameters so that the splitting does not produce a from the best valid solution
		// so that the splitting does not produce a from the best valid solution
		temp = params->penalityCapa;
		params->penalityCapa = 10000;
		education(bestValide);
		// le trainer a gard� les infos des routes de bestValide
		loc = trainer->localSearchList[0];
		params->penalityCapa = temp;

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
		for (int k = 1; k <= params->nbDays; k++)
			for (int i = 0; i < (int)loc->routes[k].size(); i++)
				if (!loc->depots[k][i]->suiv->estUnDepot)
					compteur++;
		myfile << compteur << endl;

		// exporting the total CPU time (ms)
		myBuff = new char[100];
		myfile <<"Total Time: ";sprintf(myBuff, "%d", (int)(clock() / 1000000));
		myfile << myBuff << endl;

		myBuff = new char[100];
		myfile <<"PITime: ";sprintf(myBuff, "%d", (int)(params->debut / 1000000));
		myfile << myBuff << endl;

		// exporting the time to best solution
		myBuff = new char[100];
		myfile <<"Best Solution Time: ";sprintf(myBuff, "%d", (int)(timeBest / 1000000));
		myfile << myBuff << endl;

		for (int k = 1; k <= params->nbDays; k++)
		{
			for (int i = 0; i < (int)loc->routes[k].size(); i++)
			{
				compteur = 1;
				if (!loc->depots[k][i]->suiv->estUnDepot)
				{

					myfile << " " << loc->routes[k][i]->depot->idx << " " << (k - 1) % params->ancienNbDays + 1 << " " << compteur << " " << loc->routes[k][i]->temps << " "
						   << loc->routes[k][i]->charge << " ";

					// on place la route dans un vecteur de noeuds clients
					noeudActuel = loc->routes[k][i]->depot->suiv;
					cost = params->timeCost[loc->routes[k][i]->depot->idx][noeudActuel->idx];

					rout.clear();
					routTime.clear();
					rout.push_back(loc->routes[k][i]->depot->idx);
					routTime.push_back(0);
					rout.push_back(noeudActuel->idx);
					routTime.push_back(cost);

					while (!noeudActuel->estUnDepot)
					{
						cost += params->timeCost[noeudActuel->idx][noeudActuel->suiv->idx];
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
double Population::fractionValidesCharge()
{
	int count = 0;
	for (list<bool>::iterator it = listeValiditeCharge.begin(); it != listeValiditeCharge.end(); ++it)
		if (*it == true)
			count++;

	return double(count) / (double)(100);
}

double Population::getDiversity(SousPop *pop)
{
	double total = 0;
	int count = 0;
	for (int i = 0; i < pop->nbIndiv / 2; i++)
		for (int j = i + 1; j < pop->nbIndiv / 2; j++)
		{
			total += pop->individus[i]->distance(pop->individus[j]);
			count++;
		}
	return total / (double)count;
}

double Population::getMoyenneValides()
{
	double moyenne = 0;
	for (int i = 0; i < valides->nbIndiv / 2; i++)
		moyenne += valides->individus[i]->coutSol.evaluation;
	return moyenne / (valides->nbIndiv / 2);
}
double Population::getMoyenneInvalides()
{
	double moyenne = 0;
	for (int i = 0; i < invalides->nbIndiv / 2; i++)
		moyenne += invalides->individus[i]->coutSol.evaluation;
	return moyenne / (invalides->nbIndiv / 2);
}

int Population::selectCompromis(SousPop *souspop)
{
	double pireFitnessEtendu = 0;
	int mauvais = -1;
	vector<int> classement;
	int temp, sortant;

	evalExtFit(souspop);

	// pour chaque individu on modifie le fitness etendu
	for (int i = 0; i < souspop->nbIndiv; i++)
	{
		classement.push_back(i);
		if (souspop->individus[i]->distPlusProche(1) < params->distMin)
			souspop->individus[i]->fitnessEtendu += 5;
		// for the CVRP instances, we need to allow duplicates with the same fitness since in the Golden instances
		// there is a lot of symmetry.
		if (fitExist(souspop, souspop->individus[i]))
			souspop->individus[i]->fitnessEtendu += 5;
	}

	// on classe les elements par fitness etendu et on prend le plus mauvais
	for (int n = 0; n < souspop->nbIndiv; n++)
		for (int i = 1; i < souspop->nbIndiv - n - 1; i++)
			if (souspop->individus[classement[i]]->fitnessEtendu > souspop->individus[classement[i + 1]]->fitnessEtendu)
			{
				temp = classement[i + 1];
				classement[i + 1] = classement[i];
				classement[i] = temp;
			}

	sortant = classement[souspop->nbIndiv - 1];

	if (params->rng->genrand64_real1() < -1)
		cout << souspop->individus[sortant]->fitRank << " "
			 << souspop->individus[sortant]->divRank << " "
			 << souspop->individus[sortant]->fitnessEtendu << endl;

	return sortant;
}

void Population::education(Individu *indiv)
{
	recopieIndividu(trainer, indiv);
	trainer->generalSplit_scenario();
	trainer->updateLS_scenario();
	trainer->localSearchList[0]->runSearchTotal(false);
	trainer->updateIndiv_scenario();
	recopieIndividu(indiv, trainer);
}

void Population::education_scenario(Individu *indiv)
{
	recopieIndividu(trainer, indiv);
	trainer->generalSplit_scenario();
	trainer->updateLS_scenario();
	trainer->localSearchList[0]->runSearchTotal(false);
	trainer->updateIndiv_scenario();
	recopieIndividu(indiv, trainer);
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

	cout << "It " << nbIter << " | Sol ";

	if (getIndividuBestValide() != NULL)
		cout << getIndividuBestValide()->coutSol.evaluation << " ";
	else
		cout << "NO-VALID ";

	if (getIndividuBestInvalide() != NULL)
		cout << getIndividuBestInvalide()->coutSol.evaluation;
	else
		cout << "NO-INVALID";

	cout.precision(8);

	cout << " | Moy " << getMoyenneValides() << " " << getMoyenneInvalides()
		 << " | Div " << getDiversity(valides) << " " << getDiversity(invalides)
		 << " | Val " << fractionValidesCharge()
		 << " | Pen " << params->penalityCapa << " | Pop " << valides->nbIndiv << " " << invalides->nbIndiv << endl;
}
