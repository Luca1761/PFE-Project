//
// Created by pta on 15/08/2018.
//

#ifndef IRP_PLFUNCTION_H
#define IRP_PLFUNCTION_H

#include "memory"
#include <algorithm>
#include <string>
#include <vector>
#include "LinearPiece.h"
#include "Params.h"

using namespace std;

class PLFunction
{
private:
    // instance and solve parameters
    Params *params;

public:
    // linear pieces
    vector<std::shared_ptr<LinearPiece>> pieces;

    // number of pieces
    unsigned int nbPieces;

    // default constructor
    PLFunction(Params *params);
    
    // constructor for basic F
    PLFunction(Params *params, vector<Insertion> insertions, unsigned int day, Client client, unsigned int idxScenario);
    
    // constructor for Fk (first day)
    PLFunction(Params *params, Insertion insertion, Client client, unsigned int idxScenario);
    
    // constructor by copy
    PLFunction(PLFunction *plf);

    // for a given insertion (detour, and freeload)
    // F(q) = detour + penalityCapa * (replenishment - freeload) - q * (horizon - day) 
    double calculateF(unsigned int day, double detour, double replenishment, double freeload, unsigned int idxScenario);

    // total clear of PL function 
    void clear();
    
    // get (x, y) such that pl(x) = y with minimal y
    std::shared_ptr<LinearPiece> getMinimalPiece(double &minAt, double &minValue);

    // check intersection of two pieces and fill (x, y) with the intersection point
    bool intersect(shared_ptr<LinearPiece> lp1, shared_ptr<LinearPiece> lp2, double &x, double &y);

    // append a new linear piece to function
    void append(shared_ptr<LinearPiece> lp);

    // lp(x) = newLp(x + x_axis)
    // so newLp(x) = lp(x - x_axis)
    void shiftRight(double x_axis);

    // newLP(x) = lp(x) + inventoryCost * x
    void addHoldingf(double inventoryCost);

    // newLP(x) = lp(x) - stockoutCost * x
    void addStockoutf(double stockoutCost);

    // newLP(x) = lp(x) + y_axis
    void moveUp(double y_axis);

    // create a copy in [lb, ub] range
    std::shared_ptr<PLFunction> getInBound(double lb, double ub);

    // return a copy of the piecewise linear function
    std::shared_ptr<PLFunction> clone();

    // display the PL function
    void print();

    // destructor
    ~PLFunction();
};

#endif // IRP_PLFUNCTION_H
