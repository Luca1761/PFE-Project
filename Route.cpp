#include "Route.h"
#include "LocalSearch.h"

Route::Route(void){}

Route::Route(Params* _params, LocalSearch* _myLS, unsigned int _idx, unsigned int _day, Node * _depot, double _time, double _load, double _capacity) : 
params(_params), myLS(_myLS), idx(_idx), day(_day), depot(_depot), time(_time) , load(_load), capacity(_capacity)
{
	bestInsertion = vector<Insertion>(params->nbClients + params->nbDepots);
	nodeAndRouteTested = vector<bool>(params->nbClients + params->nbDepots, false);
}

Route::~Route(void){}

void Route::printRoute(std::ostream& file) {
	int place = 0 ;

	// we loop across all the route
	Node * node = depot;
	node->place = place;

	while (!node->isADepot || place == 0) {
		file << " " << node->idx << " ->";
		node = node->next;
		place++;
		node->place = place;
	}
	file <<" " << node->idx <<endl;
}

void Route::updateRouteData () {
	int place = 0;
	load = 0;
	time = 0;

	// we loop across all the route
	Node* node = depot;
	node->place = place;
	depot->previousLoad = 0;

	while (!node->isADepot || place == 0) {
		node = node->next;
		place++;
		node->place = place;
		load += myLS->deliveryPerDay[day][node->idx];
		time += params->timeCost[node->prev->idx][node->idx];
		node->previousLoad = load;
	}

	initiateInsertions();
}

void Route::evalInsertClient(Node* U)  {
	Node* currentNode;
	double cost;
	bestInsertion[U->idx].detour = 1.e30;
	bestInsertion[U->idx].place = NULL;
	bestInsertion[U->idx].load = -1.e30;
	
	bool firstIt = true ;
	if (U->route != this || !U->isPresent) {
		bestInsertion[U->idx].load = std::max<double>(0.0, capacity - load);
		currentNode = depot;
		while (!currentNode->isADepot || firstIt == true) {
			firstIt = false;
			cost = params->timeCost[currentNode->idx][U->idx] 
			+ params->timeCost[U->idx][currentNode->next->idx] 
			- params->timeCost[currentNode->idx][currentNode->next->idx];
			
			if (cost < bestInsertion[U->idx].detour - 0.0001) { 
				bestInsertion[U->idx].detour = cost;
				bestInsertion[U->idx].place = currentNode;
			}
			currentNode = currentNode->next;
		}
	} else // U is already in the route
	{
		bestInsertion[U->idx].load = std::max<double>(0.0, capacity + myLS->deliveryPerDay[day][U->idx] - load);
		bestInsertion[U->idx].detour = params->timeCost[U->prev->idx][U->idx] - params->timeCost[U->prev->idx][U->next->idx]   
										+ params->timeCost[U->idx][U->next->idx];
		bestInsertion[U->idx].place = U->prev;

		// however, we'll see if there's a better insertion possible
		// temporarily we'll remove the node from the list (in O(1))
		U->prev->next = U->next;
		U->next->prev = U->prev;
		currentNode = depot;

		// and explore the route, again
		while (!currentNode->isADepot || firstIt == true ) {
			firstIt = false;
			cost = params->timeCost[currentNode->idx][U->idx] 
			+ params->timeCost[U->idx][currentNode->next->idx]
			- params->timeCost[currentNode->idx][currentNode->next->idx];

			// finally, check if we can find a best insertion
			if (cost < bestInsertion[U->idx].detour - 0.0001) { 
				bestInsertion[U->idx].detour = cost;
				bestInsertion[U->idx].place = currentNode;
			}
			currentNode = currentNode->next ;
		}
		
		// we replace the node
		U->prev->next = U;
		U->next->prev = U;
	}
}

void Route::initiateInsertions() {
	for (unsigned int i = 0 ; i < params->nbClients + params->nbDepots ; i++) {
		bestInsertion[i].detour = 1.e30;
		bestInsertion[i].load = -1.e30;
		bestInsertion[i].place = NULL;
	}
}

void Route::reinitSingleDayMoves() {
	nodeAndRouteTested = std::vector<bool>(params->nbClients + params->nbDepots, false);
}
