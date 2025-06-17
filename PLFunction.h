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
    Params *params;

public:
    vector<std::shared_ptr<LinearPiece>> pieces;
    int nbPieces;

    std::shared_ptr<LinearPiece> minimalPiece;
    double maxValue, minValue;
    PLFunction(Params *params);

    // initialize a PL function from arc profile
    PLFunction(Params *params, vector<Insertion> insertions, int day, int client);

    // initialize a PLFunction from list of pieces
    PLFunction(Params *params, vector<shared_ptr<LinearPiece>> pieces);

    PLFunction(PLFunction *plf);

    double cost(double x);
    void clear();

    // get piece that fits with time t
    std::shared_ptr<LinearPiece> getPiece(double t);
    
    std::shared_ptr<LinearPiece> getMinimalPiece(int client, double &minAt, double &minValue);

    // check intersection of two pieces
    bool intersect(shared_ptr<LinearPiece> lp1, shared_ptr<LinearPiece> lp2, double &x, double &y);

    void append(shared_ptr<LinearPiece> lp);

    void shiftRight(double x_axis);
    void addHoldingf(double InventoryCost);
    void addStockoutf(double stockoutCost);
    void moveUp(double y_axis);

    std::shared_ptr<PLFunction> getInBound(double lb, double ub, bool updateValueAt0 = false);

    double calculateGFunction(int day, int client, double detour, double replenishment, double freeload);

    void print();

    ~PLFunction();
};

#endif // IRP_PLFUNCTION_H
