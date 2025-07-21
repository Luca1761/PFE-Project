#include "Params.h"

// creating the parameters from the instance file
Params::Params(string nomInstance, string nomSolution, int nbVeh, int seedRNG, int rou, int var, double randomValue, int idxScenario) : 
	nbVehiculesPerDep(nbVeh), seed(seedRNG), ancienNbDays(nbDays), pathToInstance(nomInstance), pathToSolution(nomSolution), idxScenario(idxScenario)
{

	if (seed == 0)
		rng = new Rng((unsigned long long)time(NULL));
	else
		rng = new Rng((unsigned long long)(seed + idxScenario));	

	debut = clock();
	fichier.open(nomInstance.c_str());

	if (fichier.is_open())
		preleveDonnees(nomInstance, 1000000);	
	else
	{
		cout << "Unable to open file : " << nomInstance << endl;
		throw string(" Unable to open file ");
	}	

	adjustDemands(randomValue, var);

	setMethodParams();

	computeDistancierFromCoords();

	calculeStructures();	

	// Compute the constant value in the objective function
	computeConstant_stockout();

	
}

void Params::adjustDemands(double rv, int var) {
    mt19937 gen(seed + static_cast<int>(rv * 10000)); 
	
    for (int i = 0; i < nbClients + nbDepots; i++) {
		if (cli[i].custIdx == 0) continue; 
		std::cout << "Client " << i << ": ";
		std::cout << "original " << cli[i].dailyDemand[1] << " // " << cli[i].maxInventory << " // ";
        for (int k = 1; k <= nbDays; k++) {
			normal_distribution<double> normDist(cli[i].dailyDemand[k], (double) var);    
            double x = normDist(gen);        // x ~ N(0,1)
            x = max<double>(0.0, x);   
			x = min<double>(x, cli[i].maxInventory);
            cli[i].dailyDemand[k] = round(x);
			std::cout << cli[i].dailyDemand[k] << " ";
        }
		std::cout << std::endl;
    }
}


Params::Params(string nomInstance, string nomSolution, int nbVeh, string nomBKS, int seedRNG) : 
	nbVehiculesPerDep(nbVeh)
{
	seed = seedRNG;
	if (seed == 0)
		rng = new Rng((unsigned long long)time(NULL));
	else
		rng = new Rng((unsigned long long)(seed));

	pathToInstance = nomInstance;
	pathToSolution = nomSolution;
	pathToBKS = nomBKS;

	debut = clock();
	nbVehiculesPerDep = nbVeh;

	// ouverture du fichier en lecture
	fichier.open(nomInstance.c_str());

	// parsing des donn�es
	if (fichier.is_open())
		preleveDonnees(nomInstance, 0);
	else
	{
		cout << "Unable to open file : " << nomInstance << endl;
		throw string(" Unable to open file ");
	}
	cout << "read file done" << endl;

	// Setting the parameters
	setMethodParams();

	// Simply compute the distances from the coordinates
	computeDistancierFromCoords();

	// calcul des structures
	calculeStructures();

	// for instances 27-32 of the PVRP, there is a small patch to avoid errors because the objective function
	// is of a very different magnitude, need to give a more adapted starting value for the penalties
	
}


void Params::computeDistancierFromCoords()
{
	double d;
	vector<double> dist;

	// on remplit les distances dans timeCost
	for (int i = 0; i < nbClients + nbDepots; i++)
	{
		dist.clear();
		for (int j = 0; j < nbClients + nbDepots; j++)
		{
			d = sqrt((cli[i].coord.x - cli[j].coord.x) * (cli[i].coord.x - cli[j].coord.x) +
					 (cli[i].coord.y - cli[j].coord.y) * (cli[i].coord.y - cli[j].coord.y));

			// integer rounding
			if (isRoundingInteger) // to be able to deal with other rounding conventions
			{
				d += 0.5;
				d = (double)(int)d;
			}
			else if (isRoundingTwoDigits)
			{
				d = d * 100.0;
				d += 0.5;
				d = (double)(int)d;
				d = d * 0.01;
			}

			dist.push_back((double)d);
		}
		timeCost.push_back(dist);
	}
}

void Params::setMethodParams()
{
	// parameters related to how the problem is treated
	isRoundingInteger = false; // using the rounding (for now set to true, because currently testing on the instances of Uchoa)
	isRoundingTwoDigits = false;

	// parameters of the population
	el = 3;					// ***
	mu = 12;				// *
	lambda = 25;			// *
	nbCountDistMeasure = 5; // o
	rho = 0.30;				// o
	delta = 0.001;			// o

	// parameters of the mutation
	prox = 40;			 // granularity parameter (expressed as a percentage of the problem size -- 35%) // ***
	proxCst = 1000000;	 // granularity parameter (expressed as a fixed maximum)
	pRep = 0.5;			 // probability of repair // o

	// parametres lies aux pEnalites adaptatives
	penalityCapa = 50;
	minValides = 0.2;	// Target range for the number of feasible solutions // **
	maxValides = 0.25;	// Target range for the number of feasible solutions // **
	distMin = 0.01;		// Distance in terms of objective function under which the solutions are considered to be the same // o
	borneSplit = 2.0;	// Split parameter (how far to search beyond the capacity limit) // o
}

Params::~Params(void) {}

// effectue le prelevement des donnees du fichier
void Params::preleveDonnees(string nomInstance, int rou)
{
	// variables temporaires utilisees dans la fonction
	vector<Vehicle> tempI;
	double capacity;
	//C. Archetti, L. Bertazzi, G. Laporte, and M.G. Speranza. A branch-and-cut algorithm for a vendor-managed inventory-routing problem. Transportation Science, 41:382-391, 2007 Instances
	// IRP format of Archetti http://or-brescia.unibs.it/instances 
	cout << "path: " << nomInstance << endl;
	if (nbVehiculesPerDep == -1) {
		cout << "ERROR : Need to specify fleet size" << endl;
		throw string("ERROR : Need to specify fleet size");
	}
	fichier >> nbClients >> nbDays >> capacity;
	bool isDSIRP = false;
	if (!isDSIRP)
		nbClients--; // the given number of nodes also counts the supplier/depot
	nbDepots = 1;
	ancienNbDays = nbDays;

	ordreVehicules.push_back(tempI);
	nombreVehicules.push_back(0);
	for (int kk = 1; kk <= nbDays; kk++) {
		ordreVehicules.push_back(tempI);
		nombreVehicules.push_back(nbDepots * nbVehiculesPerDep);
		for (int i = 0; i < nbDepots; i++)
			for (int j = 0; j < nbVehiculesPerDep; j++)
				ordreVehicules[kk].push_back(Vehicle(i, capacity));
	}
	// Liste des clients
	for (int i = 0; i < nbClients + nbDepots; i++) {
		if (isDSIRP)
			cli.push_back(getClientDSIRP(i,rou));
		else
			cli.push_back(getClient(i,rou));
	}
}

// calcule les autres structures du programme
void Params::calculeStructures()
{
	int temp;

	// initializing the correlation matrix
	isCorrelated1 = vector<vector<bool>>(nbClients + nbDepots, vector<bool>(nbClients + nbDepots, false));

	for (int i = 0; i < nbClients + nbDepots; i++)
	{
		cli[i].ordreProximite.clear();
		cli[i].sommetsVoisins.clear();
	}

	// on remplit la liste des plus proches pour chaque client
	for (int i = nbDepots; i < nbClients + nbDepots; i++)
	{
		for (int j = 0; j < nbClients + nbDepots; j++)
			if (i != j)
				cli[i].ordreProximite.push_back(j); 

		// et on la classe
		for (int a1 = 0; a1 < nbClients + nbDepots - 1; a1++)
			for (int a2 = 0; a2 < nbClients + nbDepots - a1 - 2; a2++)
				if (timeCost[i][cli[i].ordreProximite[a2]] > timeCost[i][cli[i].ordreProximite[a2 + 1]])
				{
					temp = cli[i].ordreProximite[a2 + 1];
					cli[i].ordreProximite[a2 + 1] = cli[i].ordreProximite[a2];
					cli[i].ordreProximite[a2] = temp;
				}

		// on remplit les x% plus proches
		for (int j = 0; j < (prox * (int)cli[i].ordreProximite.size()) / 100 && j < proxCst; j++)
		{
			cli[i].sommetsVoisins.push_back(cli[i].ordreProximite[j]);
			isCorrelated1[i][cli[i].ordreProximite[j]] = true;
		}
	}

	// on melange les proches
	shuffleProches();
}

// sous routine du prelevement de donnees
Client Params::getClient(int i,int rou)
{
	struct couple coordonnees;
	Client client;

	// file format of Cordeau et al.
	fichier >> client.custIdx;
	client.custIdx--;
	fichier >> coordonnees.x >> coordonnees.y;
	client.coord = coordonnees;

	if (client.custIdx == 0) // information of the supplier
	{
		double initInventory;
		double dailyProduction;
		fichier >> initInventory;
		fichier >> dailyProduction;
		availableSupply = vector<double>(nbDays + 1, 0.); // days are indexed from 1 ... t
		for (int t = 1; t <= nbDays; t++)
			availableSupply[t] = dailyProduction;
		availableSupply[1] += initInventory;
		fichier >> inventoryCostSupplier;
	} else //information of each customer
	{
		fichier >> client.startingInventory;
		fichier >> client.maxInventory;
		fichier >> client.minInventory;
		double myDailyDemand;
		fichier >> myDailyDemand;
		client.dailyDemand = vector<double>(nbDays + 1, myDailyDemand);
		fichier >> client.inventoryCost;
		client.stockoutCost = client.inventoryCost * rou;
	}

	return client;
}

Client Params::getClientDSIRP(int i,int rou)
{
	struct couple coordonnees;
	Client client;

	// file format of Cordeau et al.
	fichier >> client.custIdx;
	fichier >> coordonnees.x >> coordonnees.y;
	client.coord = coordonnees;

	if (client.custIdx == 0) // information of the supplier
	{
		double initInventory;
		double dailyProduction;
		fichier >> initInventory;
		availableSupply = vector<double>(nbDays + 1, 0.); // days are indexed from 1 ... t
		for (int t = 1; t <= nbDays; t++) {
			fichier >> dailyProduction;
			availableSupply[t] = dailyProduction;
		}
		availableSupply[1] += initInventory;
		fichier >> inventoryCostSupplier;
	} else //information of each customer
	{
		fichier >> client.startingInventory;
		client.minInventory = 0.0;
		// fichier >> client.minInventory;
		int a;
		for (unsigned int i = 0; i < 50; i++) {
			fichier >> a;
		}
		double myDailyDemand;
		client.dailyDemand = vector<double>(nbDays + 1, 0.0);
		for (unsigned int t = 1; t <= nbDays; t++) {
			fichier >> client.dailyDemand[t];
		}
		fichier >> client.maxInventory;
		fichier >> client.inventoryCost;
		fichier >> client.stockoutCost;
	}

	return client;
}

void Params::shuffleProches()
{
	int temp, temp2;

	// on introduit du desordre dans la liste des plus proches pour chaque client
	for (int i = nbDepots; i < nbClients + nbDepots; i++) {
		// on introduit du desordre
		for (int a1 = 0; a1 < (int)cli[i].sommetsVoisins.size() - 1; a1++) {
			temp2 = a1 + rng->genrand64_int64() % ((int)cli[i].sommetsVoisins.size() - a1);
			temp = cli[i].sommetsVoisins[a1];
			cli[i].sommetsVoisins[a1] = cli[i].sommetsVoisins[temp2];
			cli[i].sommetsVoisins[temp2] = temp;
		}
	}
}

Params::Params(Params *params, int decom, int *serieVisites, Vehicle **serieVehicles, int *affectDepots, int *affectPatterns, int depot, int jour, int nbVisites, int nbVeh)
{
	debut = clock();
	rng = params->rng;
	seed = params->seed;

	/* For now I just kept the CVRP decomposition */
	if (decom == 0)
		decomposeRoutes(params, serieVisites, serieVehicles, nbVisites, nbVeh);
	else
		cout << "Error : Transformation not available actually" << endl;

	// affectation des autres parametres
	setMethodParams();
	penalityCapa = params->penalityCapa;

	mu = params->mu;
	lambda = params->lambda;
	el = params->el;
	nbCountDistMeasure = params->nbCountDistMeasure;

	// calcul des distances
	// on remplit les distances dans timeCost
	for (int i = 0; i < nbClients + nbDepots; i++)
	{
		timeCost.push_back(vector<double>());
		for (int j = 0; j < nbClients + nbDepots; j++)
			timeCost[i].push_back(params->timeCost[correspondanceTable[i]][correspondanceTable[j]]);
	}

	// calcul des structures
	calculeStructures();
}

void Params::decomposeRoutes(Params *params, int *serieVisites, Vehicle **serieVehicles, int nbVisites, int nbVeh)
{
	vector<Vehicle> temp;

	correspondanceTable2.clear();
	for (int i = 0; i < params->nbClients + params->nbDepots; i++)
		correspondanceTable2.push_back(-1);

	correspondanceTable.push_back(0);

	for (int i = 0; i < nbVisites; i++)
	{
		correspondanceTable.push_back(serieVisites[i]);
		correspondanceTable2[serieVisites[i]] = (int)correspondanceTable.size() - 1;
	}

	nbVehiculesPerDep = nbVeh;
	nbClients = correspondanceTable.size() - 1;
	nbDepots = 1;
	nbDays = 1;
	ancienNbDays = 1;
	borneSplit = 1.5;

	// on place les donn�es sur les v�hicules
	ordreVehicules.clear();
	nombreVehicules.clear();
	ordreVehicules.push_back(temp);
	nombreVehicules.push_back(0);
	ordreVehicules.push_back(temp);
	nombreVehicules.push_back(nbVeh);
	for (int v = 0; v < nbVeh; v++)
	{
		ordreVehicules[1].push_back(*serieVehicles[v]);
	}

	// on met les bons clients
	// ils ont toujours les anciens num�ros
	for (int i = 0; i < nbDepots + nbClients; i++)
		cli.push_back(params->cli[correspondanceTable[i]]);
}

void Params::computeConstant_stockout() {
	objectiveConstant_stockout = 0.;
	// Adding the total cost for supplier inventory over the planning horizon (CONSTANT IN OBJECTIVE)
	for (unsigned int k = 1; k <= ancienNbDays; k++) {
		objectiveConstant_stockout += availableSupply[k] * (ancienNbDays + 1 - k) * inventoryCostSupplier;
	}
}