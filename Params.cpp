#include "Params.h"

double computeMean(const std::vector<double>& data) {
    if (data.empty()) {
        return 0.0;
    }
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    return sum / (double) data.size();
}

double computeStd(const std::vector<double>& data) {
    if (data.empty()) {
        return 0.0; 
    }
    double mean = computeMean(data); // Reusing the mean function
    double sumSquaredDifferences = 0.0;
    for (double value : data) {
        sumSquaredDifferences += std::pow(value - mean, 2);
    }
    return std::sqrt(sumSquaredDifferences / (double) data.size()); // Population standard deviation

}

Params::Params(string nomInstance, int seedRNG, unsigned int nbCore, unsigned int nbScenario, unsigned int nbVeh, bool endInventories, bool trace, bool determinist, bool trueDemand1) : 
	seed(seedRNG), nbCores(nbCore), nbScenarios(nbScenario), nbVehiculesPerDep(nbVeh), endDayInventories(endInventories), traces(trace), deterministic(determinist), trueDemandDay1(trueDemand1)
{
	jVal = 1; // initial value for jVal
	if (seed == 0)
		rng = new Rng((unsigned long long)time(NULL));
	else
		rng = new Rng((unsigned long long)(seed));	

	file.open(nomInstance.c_str());

	if (file.is_open())
		collectData();	
	else {
		cout << "Unable to open file : " << nomInstance << endl;
		throw string(" Unable to open file ");
	}
	setMethodParams();
	computeDistancesFromCoords();
}

void Params::computeScenarios() {
	// clean the demand and supply
	for (unsigned int i = nbDepots; i < nbClients + nbDepots; i++) {
		cli[i].dailyDemand = vector<vector<double>>(nbScenarios, vector<double>(nbDays + 1, 0.0));
	}
	availableSupply = vector<double>(nbDays + 1 + 1, 0.);

	// supply is deterministic and availableSupply[1] depends on jVal (allSupply contains supply for the whole horizon)
	for (unsigned int day = 1; day <= nbDays; day++) {
		availableSupply[day] = allSupply[jVal + day - 1];
	}
	availableSupply[1] += cli[0].startingInventory; // at day 1, we add the initial inventory

	if (traces) {
		std::cout << "Depot supply: " << std::endl;
		for (unsigned int day = 1; day <= nbDays; day++) {
			std::cout << availableSupply[day] << " ";
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}

	if (traces) std::cout << "Scenarios for customers demands: " << std::endl;
	mt19937 gen((unsigned int) seed + jVal); 
	for (unsigned int scenario = 0; scenario < nbScenarios; scenario++) {
		if (traces) std::cout << "Scenario " << scenario + 1 << " :" << std::endl;
		for (unsigned int i = nbDepots; i < nbClients + nbDepots; i++) {	
			if (traces) std::cout << "Client " << i << ": ";
			normal_distribution<double> normDist(meanDemands[i - nbDepots], stdDemands[i - nbDepots]); // normal distribution with mean and std 
			for (unsigned int day = 1; day <= nbDays; day++) {
				double x = max<double>(0.0, normDist(gen)); // make sure the demand is >= 0
				if ((trueDemandDay1 && day == 1) || deterministic) cli[i].dailyDemand[scenario][day] = cli[i].trueDemand[jVal + day - 1]; // use the true demand (for deterministic or for day 1)
				else cli[i].dailyDemand[scenario][day] = (double) round(x); 
				if (traces) std::cout << cli[i].dailyDemand[scenario][day] << " ";
        	}
			if (traces) std::cout << std::endl;
		}
	}
	if (traces) std::cout << std::endl;
}

void Params::computeDistancesFromCoords() {
	double cij;
	vector<double> distances;

	// fille timeCost distances
	for (unsigned int i = 0; i < nbClients + nbDepots; i++) {
		distances.clear();
		for (unsigned int j = 0; j < nbClients + nbDepots; j++) {
			cij = sqrt((cli[i].coord.x - cli[j].coord.x) * (cli[i].coord.x - cli[j].coord.x) +
					 (cli[i].coord.y - cli[j].coord.y) * (cli[i].coord.y - cli[j].coord.y));

			// integer rounding
			if (isRoundingInteger) // to be able to deal with other rounding conventions (as in coelho instances)
			{
				cij += 0.5;				 // even if I'm pretty sure it doesn't respect triangular inequality
				cij = (double)(int)cij;
			}
			else if (isRoundingTwoDigits)
			{
				cij = cij * 100.0;
				cij += 0.5;
				cij = (double)(int)cij;
				cij = cij * 0.01;
			}

			distances.push_back((double) cij);
		}
		timeCost.push_back(distances);
	}
}

void Params::setMethodParams()
{
	// parameters related to how the problem is treated
	isRoundingInteger = true; // using the rounding (for now set to true, because currently testing on the instances of Coelho)
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

	// adaptative penalities related parameters
	penalityCapa = std::vector<double> (nbScenarios, 50.0);
	minFeasibles = 0.2;	// Target range for the number of feasible solutions // **
	maxFeasibles = 0.25;	// Target range for the number of feasible solutions // **
	distMin = 0.01;		// Distance in terms of objective function under which the solutions are considered to be the same // o
	splitBounds = std::vector<double> (nbScenarios, 1.0);	// Split parameter (how far to search beyond the capacity limit) // o
}


void Params::collectData() {	
	// temporary variables
	vector<Vehicle> tempI;
	double vehicleCapacity;
	if (nbVehiculesPerDep == 0) {
		cout << "ERROR : Need to specify fleet size" << endl;
		throw string("ERROR : Need to specify fleet size");
	}
	
	file >> nbClients >> pHorizon >> vehicleCapacity;
	nbDepots = 1; // in the IRP, only one depot
	
	vehicleOrder.push_back(tempI); // (vehicleOrder[0] and vehicleNumber[0] are not used)
	vehicleNumber.push_back(0);
	
	for (unsigned int kk = 1; kk <= pHorizon; kk++) {
		vehicleOrder.push_back(tempI);
		vehicleNumber.push_back(nbDepots * nbVehiculesPerDep);
		for (unsigned int i = 0; i < nbDepots; i++)
			for (unsigned int j = 0; j < nbVehiculesPerDep; j++)
				vehicleOrder[kk].push_back(Vehicle(i, vehicleCapacity));
	}
	// get client (and supplier) information
	for (unsigned int i = 0; i < nbClients + nbDepots; i++)
		cli.push_back(getNextClient());
}

void Params::computeStructures() {
	// initializing the correlation matrix
	isCorrelated = vector<vector<bool>>(nbClients + nbDepots, vector<bool>(nbClients + nbDepots, false));
	
	// clean the old structures (in case of rolling horizon, proximity order is still the same but not neighbors because we shuffle it at the beginning)
	for (unsigned int i = 0; i < nbClients + nbDepots; i++) {
		cli[i].proximityOrder.clear();
		cli[i].neighbors.clear();
	}
	
	for (unsigned int i = nbDepots; i < nbClients + nbDepots; i++) {
		// we fill the vector of nearest for every client
		for (unsigned int j = 0; j < nbClients + nbDepots; j++)
			if (i != j) cli[i].proximityOrder.push_back(j); 

		// then we sort them
		for (unsigned int a1 = 0; a1 < nbClients + nbDepots - 1; a1++)
			for (unsigned int a2 = 0; a2 < nbClients + nbDepots - a1 - 2; a2++)
				if (timeCost[i][cli[i].proximityOrder[a2]] > timeCost[i][cli[i].proximityOrder[a2 + 1]]) {
					std::swap(cli[i].proximityOrder[a2 + 1], cli[i].proximityOrder[a2]);
				}
		
		// we only keep the prox% nearest neighbors
		for (unsigned int j = 0; j < (prox * (unsigned int) cli[i].proximityOrder.size()) / 100 && j < proxCst; j++) {
			cli[i].neighbors.push_back(cli[i].proximityOrder[j]);
			isCorrelated[i][cli[i].proximityOrder[j]] = true;
		}
	}
		
	// we shuffle these neighbors
	shuffleNeighbors();
}
	
Client Params::getNextClient() {
	struct couple coordonnees;
	Client client;
	
	// file format of Coelho et al.
	file >> client.custIdx;
	file >> coordonnees.x >> coordonnees.y;
	client.coord = coordonnees;
	
	if (client.custIdx == 0) // information of the supplier
	{
		file >> client.startingInventory;
		allSupply = vector<double>(pHorizon + 1, 0.);
		for (unsigned int t = 1; t <= pHorizon; t++) {
			file >> allSupply[t];
		}
		file >> inventoryCostSupplier;
	} else //information of each customer
	{
		file >> client.startingInventory;

		std::vector<double> historicalDemands = vector<double>(50, 0.);
		for (unsigned int i = 0; i < 50; i++) {
			file >> historicalDemands[i];
		}
		meanDemands.push_back(computeMean(historicalDemands));
		stdDemands.push_back(computeStd(historicalDemands));
		oldDemands.push_back(historicalDemands);

		client.trueDemand = vector<double>(pHorizon + 1, 0.0);
		for (unsigned int t = 1; t <= pHorizon; t++) {
			file >> client.trueDemand[t];
		}

		file >> client.maxInventory >> client.inventoryCost >> client.stockoutCost;
	}
	return client;
}
		
void Params::shuffleNeighbors() {
	// we shuffle the nearest neighbors vector of each client
	for (unsigned int i = nbDepots; i < nbClients + nbDepots; i++) {
		for (unsigned int a1 = 1; a1 < cli[i].neighbors.size(); a1++) {
			unsigned int temp = a1 - 1 + (unsigned int) (rng->genrand64_int64() % (cli[i].neighbors.size() - a1 + 1));
			std::swap(cli[i].neighbors[a1], cli[i].neighbors[temp]);
		}
	}
}
		
void Params::computeObjectiveConstant() {
	objectiveConstant = 0.;
	// Adding the total cost for supplier inventory over the planning horizon (CONSTANT IN OBJECTIVE)
	for (unsigned int k = 1; k <= nbDays; k++) {
		objectiveConstant += availableSupply[k] * (nbDays + 1 - k) * inventoryCostSupplier;
	}
}

void Params::setJval(unsigned int _jVal) {
	jVal = _jVal;
	nbDays = pHorizon - jVal + 1;
}

void Params::updateToDay(unsigned int j, std::vector<double> deliveries) {
	setJval(j); // set jVal to the actual day of the rolling horizon
	if (rng != nullptr) delete rng;
	rng = new Rng((unsigned long long)(seed + jVal));  // new rng to compute scenarios, etc.
	if (j != 1) cli[0].startingInventory = availableSupply[1] - std::accumulate(deliveries.begin(), deliveries.end(), 0.0); // if not first day of rolling horizon, actualize starting inventory of supplier
	if (traces) std::cout << " Client | Init | MAX | Previous delivery | True demand" << std::endl;
    for (unsigned int c = nbDepots; c < nbDepots + nbClients; c++) {  
		string deliver = (j == 1) ? "?" : to_string((int) deliveries[c - nbDepots]);
		if (j != 1) cli[c].startingInventory = std::max(0.0, cli[c].startingInventory + deliveries[c - nbDepots] - cli[c].trueDemand[j - 1]); // same for customers
      if (traces) std::cout << "Client " << c << ": " <<  cli[c].startingInventory << " " << cli[c].maxInventory << " " << deliver << " "  << cli[c].trueDemand[j] << std::endl;
    }
	if (traces) std::cout << std::endl;
	if (j != 1) {
		for (unsigned int i = nbDepots; i < nbClients + nbDepots; i++) {
			oldDemands[i - nbDepots].push_back(cli[i].trueDemand[j - 1]); // we add the previous demand to historical data and compute again means and std
			meanDemands[i - nbDepots] = computeMean(oldDemands[i - nbDepots]);
			stdDemands[i - nbDepots] = computeStd(oldDemands[i - nbDepots]);
		}
	}
    
    computeScenarios(); // compute new scenarios (and actualize available supply)
    setMethodParams(); // reset method params (especially borneSplit and penalityCapa)
    computeStructures();  // re compute structures, especially neighbors, randomly shuffled again
    computeObjectiveConstant(); // compute the new objective constant, related to new available supply 
}

Params::~Params(void) {}