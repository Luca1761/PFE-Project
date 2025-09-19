//
// Created by pta on 15/08/2018.
//

#include "PLFunction.h"
#include <cmath>
#include <cassert>

PLFunction::PLFunction(Params *_params) : params(_params) {
    nbPieces = 0;
    pieces = vector<shared_ptr<LinearPiece>>();
}

PLFunction::PLFunction(PLFunction *plf) : params(plf->params) {
    nbPieces = plf->nbPieces;
    pieces = plf->pieces;
}

void PLFunction::clear() {
    nbPieces = 0;
    pieces.clear();
}

PLFunction::PLFunction(Params *_params, Insertion insertion, Client client, unsigned int idxScenario) : params(_params){
    double maxInventory = (params->endDayInventories) ? client.maxInventory + client.theoricalMinDemand
                                                      : client.maxInventory;
    nbPieces = 0;
    pieces = vector<shared_ptr<LinearPiece>>();
    double lim = std::min(insertion.load, maxInventory);
    double costNoPenality = insertion.detour - params->inventoryCostSupplier * (double)(params->nbDays) * lim;
    shared_ptr<LinearPiece> tmp(make_shared<LinearPiece>(0, insertion.detour, lim, costNoPenality));
    tmp->fromInst = make_shared<Insertion>(insertion.detour, insertion.load, insertion.place);
    append(tmp);
    if (insertion.load < maxInventory) {
        double costPenality = insertion.detour 
                            + params->penalityCapa[idxScenario] * std::max(0.0, maxInventory - insertion.load) 
                            - params->inventoryCostSupplier * (double)(params->nbDays) * (maxInventory);

        shared_ptr<LinearPiece> tmp2(make_shared<LinearPiece>(insertion.load, costNoPenality, maxInventory, costPenality));
        tmp2->fromInst = make_shared<Insertion>(insertion.detour, insertion.load, insertion.place);
        append(tmp2);
    }
}


PLFunction::PLFunction(Params *_params, vector<Insertion> insertions, unsigned int day, Client client, unsigned int idxScenario) : params(_params) {
    if (insertions.empty()) {
        std::cout << "Cannot initialize a PL function with zero element in insertions!!!" << std::endl;
        throw std::string("Cannot initialize a PL function with zero element in insertions");
    }

    clear();

    double pre_x = 0;
    double pre_y, x, y, pre_load, pre_detour;
    Node *pre_place = nullptr;
    
    // loop through all pieces
    std::vector<Insertion>::iterator index = insertions.begin();
    while (index != insertions.end()) {
        if (index->detour < -1) {
            std::cout << "Negative detour - Triangular inequality?" << std::endl;
            throw std::string("Negative detour - Triangular inequality?");
        }
        if (index == insertions.begin()) {
            pre_x = 0;
            pre_y = calculateF(day, index->detour, pre_x,  pre_x, idxScenario);
            x = index->load;
            y = calculateF(day, index->detour, x, index->load, idxScenario);
        } else {
            double current_cost = calculateF(day, index->detour, index->load, index->load, idxScenario);
            std::vector<Insertion>::iterator pre_index = index - 1;
            double pre_insertion_cost_with_current_demand = calculateF(day, pre_index->detour, index->load, pre_index->load, idxScenario);

            bool is_dominated_vehicle = (index->load <= pre_load) || (pre_insertion_cost_with_current_demand <= current_cost);

            if (!is_dominated_vehicle) {    // make piece
                x = pre_load + (index->detour - pre_detour) / params->penalityCapa[idxScenario];

                y = calculateF(day, pre_detour, x, pre_load, idxScenario);

                shared_ptr<LinearPiece> tmp(make_shared<LinearPiece>(pre_x, pre_y, x, y));
                tmp->fromInst = make_shared<Insertion>(pre_index->detour, pre_index->load, pre_index->place);
                append(tmp);   
                    
                pre_x = x;
                pre_y = y;
            }   

            x = index->load;
            y = current_cost;
        }

        // make piece
        pre_y = calculateF(day, index->detour, pre_x, index->load, idxScenario);
        shared_ptr<LinearPiece> tmp(make_shared<LinearPiece>(pre_x, pre_y, x, y));
        tmp->fromInst = make_shared<Insertion>(index->detour, index->load, index->place);

        append(tmp);

        pre_x = x;
        pre_y = y;
        pre_load = index->load;
        pre_detour = index->detour;
        pre_place = index->place;
        index++;
    }

    // the last piece
    x =  (params->endDayInventories) ? client.maxInventory + client.theoricalMinDemand : client.maxInventory;
    if (x - pre_x <= 0.01 ) {
        return;
    }

    y = calculateF(day, pre_detour, x, pre_load, idxScenario);

    // make last piece
    shared_ptr<LinearPiece> tmp(make_shared<LinearPiece>(pre_x, pre_y, x, y));
    tmp->fromInst = make_shared<Insertion>(pre_detour, pre_load, pre_place);
    append(tmp);
}

std::shared_ptr<LinearPiece> PLFunction::getMinimalPiece(double &minAt, double &minValue){
    if (nbPieces == 0) return nullptr;
    shared_ptr<LinearPiece> minPiece = pieces[0];
    minValue = pieces[0]->p1->y;
    minAt = pieces[0]->p1->x;
    
    if (lt(minPiece->p2->y, minValue)) {
        minValue = minPiece->p2->y;
        minAt = minPiece->p2->x;
    }
    
    for (auto& piece : pieces) {
        if (gt(minValue, piece->p2->y)){ 
            minPiece = piece;
            
            minValue = piece->p2->y;
            minAt = piece->p2->x;
        }
        if (gt(minValue, piece->p1->y)){
            minPiece = piece;
            
            minValue = piece->p1->y;
            minAt = piece->p1->x;
        }
    }
    return minPiece;
}

bool PLFunction::intersect(shared_ptr<LinearPiece> lp1, shared_ptr<LinearPiece> lp2, double &x, double &y) {
    if (lp1) lp1->updateSlope();
    if (lp2) lp2->updateSlope();
    
    // parallel
    if (eq(lp1->slope, lp2->slope)) return false;
    
    // get intersection point
    double del_b = ((lp1->p2->y - lp1->slope * lp1->p2->x ) - (lp2->p2->y - lp2->slope * lp2->p2->x));
    double del_k = lp2->slope - lp1->slope;
    x = del_b / del_k;
    y = lp1->p2->y - lp1->slope * (lp1->p2->x - x);

    return le(lp1->p1->x, x) && le(x, lp1->p2->x) && le(lp2->p1->x, x) && le(x, lp2->p2->x);
}

void PLFunction::append(shared_ptr<LinearPiece> lp){
    if (!lp) return;
    double minNum = std::min(lp->p1->x, lp->p2->x);
    double maxNum = std::max(lp->p1->x, lp->p2->x);
    if (ceil(minNum) > floor(maxNum)) return;
    double exp = 0;
    if(lp->slope<100 && lp->slope>-100) exp = 0.01;
    else if(lp->slope<1000&& lp->slope>-1000) exp = 0.001;
    else exp = 0.0001;
    if (fabs(lp->p1->y*1000-floor(lp->p1->y*1000)) <= exp) lp->p1->y = floor(lp->p1->y*1000)/1000;  
    if (fabs(lp->p1->y*1000-ceil(lp->p1->y*1000)) <= exp) lp->p1->y = ceil(lp->p1->y*1000)/1000;  
    if (fabs(lp->p2->y*1000-floor(lp->p2->y*1000)) <= exp) lp->p2->y = floor(lp->p2->y*1000)/1000;  
    if (fabs(lp->p2->y*1000-ceil(lp->p2->y*1000)) <= exp) lp->p2->y = ceil(lp->p2->y*1000)/1000; 
    
    if (fabs(lp->p1->x*1000-floor(lp->p1->x*1000)) <= exp) lp->p1->x = floor(lp->p1->x*1000)/1000;  
    if (fabs(lp->p1->x*1000-ceil(lp->p1->x*1000)) <= exp) lp->p1->x = ceil(lp->p1->x*1000)/1000;  
    if (fabs(lp->p2->x*1000-floor(lp->p2->x*1000)) <= exp) lp->p2->x = floor(lp->p2->x*1000)/1000;  
    if (fabs(lp->p2->x*1000-ceil(lp->p2->x*1000)) <= exp) lp->p2->x = ceil(lp->p2->x*1000)/1000;

    
    shared_ptr<LinearPiece> newPiece = lp->clone();
    if (newPiece) newPiece->updateSlope();

    if (nbPieces == 0) {
        pieces.push_back(newPiece);
        nbPieces += 1;
    } else{
        shared_ptr<LinearPiece> lastPiece = pieces[nbPieces - 1];
        lastPiece->updateSlope();

        bool isLastPoint = eq(lastPiece->p1->x, lastPiece->p2->x);
        bool isNewPoint = eq(newPiece->p1->x,newPiece->p2->x);
        bool test = eq(newPiece->p1->y,newPiece->p2->y);
        if (isNewPoint && !test) {
            newPiece->print();
            double valMin = std::min(newPiece->p1->y, newPiece->p2->y);
            newPiece->update(newPiece->p1->x, valMin, newPiece->p2->x, valMin);
        }

        if (!lastPiece->eqDeep(newPiece)) {     
            if (isLastPoint && isNewPoint && eq(newPiece->p1->x, lastPiece->p2->x)) {
                if (gt(lastPiece->p2->y, newPiece->p1->y)) {
                    pieces.pop_back();
                    pieces.push_back(newPiece);
                }
            } else if (isLastPoint && eq(newPiece->p1->x, lastPiece->p2->x)) {
                if(ge(lastPiece->p2->y, newPiece->p1->y)){
                    pieces.pop_back();            
                } else{
                    nbPieces += 1;
                }
                pieces.push_back(newPiece);
            } else if (isNewPoint && eq(newPiece->p1->x, lastPiece->p2->x)) {
                if (gt(lastPiece->p2->y, newPiece->p1->y)){
                    pieces.push_back(newPiece);
                    nbPieces += 1;
                }
            } else{
                if(eq(lastPiece->slope, newPiece->slope) && eq(newPiece->p1->y, lastPiece->p2->y) && lastPiece->eqFromC(newPiece) && lastPiece->eqFromF(newPiece) &&lastPiece->eqFromCpre(newPiece) && lastPiece->stockout == newPiece->stockout){
                    newPiece->p1->x = lastPiece->p1->x;
                    newPiece->p1->y = lastPiece->p1->y;
                    pieces.pop_back();
                    pieces.push_back(newPiece);  
                } else{
                    pieces.push_back(newPiece);
                    nbPieces += 1;
                }
            }
        }
    }
}

void PLFunction::shiftRight(double x_axis) {
    for (auto &piece : pieces) {
        piece->p1->x += x_axis;
        piece->p2->x += x_axis;
    }
}

void PLFunction::addHoldingf(double inventoryCost) {
    for (auto &piece : pieces) {
        piece->p1->y += inventoryCost * (piece->p1->x);
        piece->p2->y += inventoryCost * (piece->p2->x); 
        piece->update(piece->p1->x, piece->p1->y, piece->p2->x, piece->p2->y);
    }
}

void PLFunction::addStockoutf(double stockoutCost) {
    for (auto &piece : pieces) {
        piece->p1->y -= stockoutCost*(piece->p1->x);
        piece->p2->y -= stockoutCost*(piece->p2->x);
        piece->update(piece->p1->x, piece->p1->y, piece->p2->x, piece->p2->y);
    }
}

void PLFunction::moveUp(double y_axis) {
    for (auto &piece : pieces) {
        piece->p1->y += y_axis;
        piece->p2->y += y_axis;
    }
}

std::shared_ptr<PLFunction> PLFunction::getInBound(double lb, double ub) {
    std::shared_ptr<PLFunction> plFunction(std::make_shared<PLFunction>(params));

    if (lt(ub, lb)) return plFunction;

    for (auto &piece : pieces){
        shared_ptr<LinearPiece> lp = piece->getInBound(lb, ub);
        if (lp != nullptr)
            plFunction->append(lp);
    }
    return plFunction;
}

std::shared_ptr<PLFunction> PLFunction::clone() {
  std::shared_ptr<PLFunction> destination(std::make_shared<PLFunction>(params));
  for (auto & piece : pieces) {
    destination->append(piece);
  }
  return destination;
}

void PLFunction::print() {
    for (unsigned int i = 0; i < nbPieces; i++) {        
        cout << i <<" :(" << pieces[i]->p1->x << ", " << pieces[i]->p1->y << ", " << pieces[i]->p2->x
             << ", " << pieces[i]->p2->y << "), ";
    }
    cout << endl;
}

double PLFunction::calculateF(unsigned int day, double detour, double replenishment, double freeload, unsigned int idxScenario) {
    // detour
    double cost = detour;

    //cost-depot(holdingcost)
    cost -= params->inventoryCostSupplier * replenishment * (double)(params->nbDays - day); 
    
    // possible excess capacity
    cost += params->penalityCapa[idxScenario] * std::max(0.0, replenishment - freeload);
    return cost;
}

PLFunction::~PLFunction() {}
