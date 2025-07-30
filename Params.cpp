#include "Params.h"

double calculateMean(const std::vector<double>& data) {
    if (data.empty()) {
        return 0.0; // Or handle error for empty data
    }
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / (double) data.size();
}

double calculateStandardDeviation(const std::vector<double>& data, bool isSample = true) {
    if (data.empty()) {
        return 0.0; // Or handle error for empty data
    }

    double mean = calculateMean(data); // Reusing the mean function

    double sumSquaredDifferences = 0.0;
    for (double value : data) {
        sumSquaredDifferences += std::pow(value - mean, 2);
    }

    if (isSample && data.size() > 1) {
        return std::sqrt(sumSquaredDifferences / ((double)data.size() - 1.0)); // Sample standard deviation
    } else if (!isSample) {
        return std::sqrt(sumSquaredDifferences / (double) data.size()); // Population standard deviation
    } else {
        return 0.0; // Handle case of single element for sample std dev
    }
}



// creating the parameters from the instance file
Params::Params(string nomInstance, int seedRNG, unsigned int nbScenario, unsigned int nbVeh, bool trace) : 
	seed(seedRNG), nbScenarios(nbScenario), nbVehiculesPerDep(nbVeh), traces(trace)
{
	if (seed == 0)
		rng = new Rng((unsigned long long)time(NULL));
	else
		rng = new Rng((unsigned long long)(seed));	

	fichier.open(nomInstance.c_str());

	if (fichier.is_open())
		preleveDonnees(nomInstance);	
	else
	{
		cout << "Unable to open file : " << nomInstance << endl;
		throw string(" Unable to open file ");
	}
	setMethodParams();
	computeDistancierFromCoords();
}

void Params::adjustDemands() {
	for (unsigned int i = nbDepots; i < nbClients + nbDepots; i++) {
		cli[i].dailyDemand = vector<vector<double>>(nbScenarios, vector<double>(nbDays + 1, 0.0));
	}
	availableSupply = vector<double>(nbDays + 1 + 1, 0.);
	for (unsigned int t = 1; t <= nbDays; t++) {
		availableSupply[t] = allSupply[jVal + t - 1];
	}
	availableSupply[1] += cli[0].startingInventory;

	for (unsigned int d = 0; d < nbDepots; d++) {
		if (traces) std::cout << "Depot " << d << ": " << std::endl;
		for (unsigned int k = 1; k <= nbDays; k++) {
			if (traces) std::cout << availableSupply[k] << " ";
		}
		if (traces) std::cout << std::endl;
	}
	mt19937 gen((unsigned int) seed + jVal); 
	for (unsigned int scenario = 0; scenario < nbScenarios; scenario++) {
		if (traces) std::cout << "Scenario " << scenario + 1 << " :" << std::endl;
		for (unsigned int i = nbDepots; i < nbClients + nbDepots; i++) {	
			if (traces) std::cout << "Client " << i << ": ";
			for (unsigned int day = 1; day <= nbDays; day++) {
				normal_distribution<double> normDist((int) meanDemands[i - nbDepots], (int) stdDemands[i - nbDepots]);    
				double x = normDist(gen);        // x ~ N(0,1)
				x = max<double>(0.0, x);   
				x = min<double>(x, cli[i].maxInventory);
				cli[i].dailyDemand[scenario][day] = (double) round(x);
				if (traces) std::cout << cli[i].dailyDemand[scenario][day] << " ";
        	}
			if (traces) std::cout << std::endl;
		}
	}
	if (traces) std::cout << std::endl;
}

void Params::computeDistancierFromCoords()
{
	double d;
	vector<double> distances;

	// on remplit les distances dans timeCost
	for (unsigned int i = 0; i < nbClients + nbDepots; i++) {
		distances.clear();
		for (unsigned int j = 0; j < nbClients + nbDepots; j++) {
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

			distances.push_back((double)d);
		}
		timeCost.push_back(distances);
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
	penalityCapa = std::vector<double> (nbScenarios, 50.0);
	minValides = 0.2;	// Target range for the number of feasible solutions // **
	maxValides = 0.25;	// Target range for the number of feasible solutions // **
	distMin = 0.01;		// Distance in terms of objective function under which the solutions are considered to be the same // o
	borneSplit = std::vector<double> (nbScenarios, 1.0);	// Split parameter (how far to search beyond the capacity limit) // o
}

Params::~Params(void) {}

// effectue le prelevement des donnees du fichier
void Params::preleveDonnees(string nomInstance)
{	
	jVal = 1;
	// variables temporaires utilisees dans la fonction
	vector<Vehicle> tempI;
	double capacity;
	//C. Archetti, L. Bertazzi, G. Laporte, and M.G. Speranza. A branch-and-cut algorithm for a vendor-managed inventory-routing problem. Transportation Science, 41:382-391, 2007 Instances
	// IRP format of Archetti http://or-brescia.unibs.it/instances 
	if (nbVehiculesPerDep == 0) {
		cout << "ERROR : Need to specify fleet size" << endl;
		throw string("ERROR : Need to specify fleet size");
	}
	fichier >> nbClients >> pHorizon >> capacity;

	nbDepots = 1;

	ordreVehicules.push_back(tempI);
	nombreVehicules.push_back(0);
	for (unsigned int kk = 1; kk <= pHorizon; kk++) {
		ordreVehicules.push_back(tempI);
		nombreVehicules.push_back(nbDepots * nbVehiculesPerDep);
		for (unsigned int i = 0; i < nbDepots; i++)
			for (unsigned int j = 0; j < nbVehiculesPerDep; j++)
				ordreVehicules[kk].push_back(Vehicle(i, capacity));
	}
	// Liste des clients
	for (unsigned int i = 0; i < nbClients + nbDepots; i++) {
		cli.push_back(getClientDSIRP());
	}
}

// calcule les autres structures du programme
void Params::calculeStructures() {
	// initializing the correlation matrix
	isCorrelated1 = vector<vector<bool>>(nbClients + nbDepots, vector<bool>(nbClients + nbDepots, false));

	for (unsigned int i = 0; i < nbClients + nbDepots; i++)
	{
		cli[i].ordreProximite.clear();
		cli[i].sommetsVoisins.clear();
	}

	// on remplit la liste des plus proches pour chaque client
	for (unsigned int i = nbDepots; i < nbClients + nbDepots; i++) {
		for (unsigned int j = 0; j < nbClients + nbDepots; j++)
			if (i != j)
				cli[i].ordreProximite.push_back(j); 

		// et on la classe
		for (unsigned int a1 = 0; a1 < nbClients + nbDepots - 1; a1++)
			for (unsigned int a2 = 0; a2 < nbClients + nbDepots - a1 - 2; a2++)
				if (timeCost[i][cli[i].ordreProximite[a2]] > timeCost[i][cli[i].ordreProximite[a2 + 1]]) {
					std::swap(cli[i].ordreProximite[a2 + 1], cli[i].ordreProximite[a2]);
				}

		// on remplit les x% plus proches
		for (unsigned int j = 0; j < (prox * cli[i].ordreProximite.size()) / 100 && j < proxCst; j++) {
			cli[i].sommetsVoisins.push_back(cli[i].ordreProximite[j]);
			isCorrelated1[i][cli[i].ordreProximite[j]] = true;
		}
	}

	// on melange les proches
	shuffleProches();
}

Client Params::getClientDSIRP()
{
	struct couple coordonnees;
	Client client;

	// file format of Cordeau et al.
	fichier >> client.custIdx;
	fichier >> coordonnees.x >> coordonnees.y;
	client.coord = coordonnees;

	if (client.custIdx == 0) // information of the supplier
	{
		fichier >> client.startingInventory;
		allSupply = vector<double>(pHorizon + 1 + 1, 0.);
		for (unsigned int t = 1; t <= pHorizon; t++) {
			fichier >> allSupply[t];
		}
		fichier >> inventoryCostSupplier;
	} else //information of each customer
	{
		std::vector<double> oldDemand = vector<double>(50, 0.);
		client.testDemand = vector<double>(pHorizon + 1, 0.0);
		fichier >> client.startingInventory;
		client.minInventory = 0.0;
		for (unsigned int i = 0; i < 50; i++) {
			fichier >> oldDemand[i];
		}
		for (unsigned int t = 1; t <= pHorizon; t++) {
			fichier >> client.testDemand[t];
		}
		fichier >> client.maxInventory;
		fichier >> client.inventoryCost;
		fichier >> client.stockoutCost;
		// TO CHECK
		meanDemands.push_back(calculateMean(oldDemand));
		stdDemands.push_back(calculateStandardDeviation(oldDemand));
		oldDemands.push_back(oldDemand);
		
	}

	return client;
}

void Params::shuffleProches() {
	// on introduit du desordre dans la liste des plus proches pour chaque client
	for (unsigned int i = nbDepots; i < nbClients + nbDepots; i++) {
		// on introduit du desordre
		for (int a1 = 0; a1 < (int)cli[i].sommetsVoisins.size() - 1; a1++) {
			int temp = a1 + rng->genrand64_int64() % ((int)cli[i].sommetsVoisins.size() - a1);
			std::swap(cli[i].sommetsVoisins[a1], cli[i].sommetsVoisins[temp]);
		}
	}
}

void Params::computeConstant_stockout() {
	objectiveConstant_stockout = 0.;
	// Adding the total cost for supplier inventory over the planning horizon (CONSTANT IN OBJECTIVE)
	for (unsigned int k = 1; k <= ancienNbDays; k++) {
		objectiveConstant_stockout += availableSupply[k] * (ancienNbDays + 1 - k) * inventoryCostSupplier;
	}
}

void Params::updateToDay(unsigned int j, std::vector<double> deliveries) {
    setJval(j);
	if (j!=1) cli[0].startingInventory = availableSupply[1] - std::accumulate(deliveries.begin(), deliveries.end(), 0.0);
	if (traces) std::cout << "Init | MAX | Previous delivery | True demand" << std::endl;
    for (unsigned int c = nbDepots; c < nbDepots + nbClients; c++) {  
      double deliver = (j == 1) ? -1 : deliveries[c - nbDepots];
      if (j!=1) cli[c].startingInventory = std::max(0.0, cli[c].startingInventory + (deliveries[c - nbDepots] - cli[c].testDemand[j - 1]));
      if (traces) std::cout << "Client " << c << ": " <<  cli[c].startingInventory << " " << cli[c].maxInventory << " " << deliver << " "  << cli[c].testDemand[j] << std::endl;
    }
	if (traces) std::cout << std::endl;
	if (j != 1) {
		for (unsigned int i = nbDepots; i < nbClients + nbDepots; i++) {
			oldDemands[i - nbDepots].push_back(cli[i].testDemand[j - 1]);
			meanDemands[i - nbDepots] = calculateMean(oldDemands[i - nbDepots]);
			stdDemands[i - nbDepots] = calculateStandardDeviation(oldDemands[i - nbDepots]);
		}
	}
    
    adjustDemands();
    setMethodParams();
    calculeStructures();	
    computeConstant_stockout();
}