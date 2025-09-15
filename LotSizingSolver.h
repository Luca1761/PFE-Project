//
// Created by pta on 16/08/2018.
//

#ifndef IRP_LOTSIZINGSOLVER_H
#define IRP_LOTSIZINGSOLVER_H

#include "memory"
#include "PLFunction.h"
#include "LinearPiece.h"
#include "Params.h"

using namespace std;

class LotSizingSolver {
private:
    // instance and solve parameters
    Params* params; 

    // index of client to reinsert
    unsigned int client;

    // horizon for resolution
    unsigned int horizon;

    //number of scenarios
    unsigned int nbScenario;
public:

    vector<vector<vector<Insertion>>> insertions;        // possible insertions for client on day t (insertions[scenario][t])     
    vector<vector<std::shared_ptr<PLFunction>>> C;       // cost-to-go functions (size: nbScenario * nbDays)
    vector<vector<double>> Ck;                           // cost-to-go values for day 1 (one per insertion)
    vector<vector<std::shared_ptr<PLFunction>>> savedCk; // save PLFunction that led to Ck values
    vector<vector<std::shared_ptr<PLFunction>>> Fk;      // F_k = detour_k + capacityPenalty_k - q * supplierCost
    vector<vector<std::shared_ptr<PLFunction>>> F;       // F(q) = min_k(detour_k + capacityPenalty_k) - q * supplierCost
    vector<vector<double>> I;                            // backtracked inventories (after DP)
    vector<vector<double>> quantities;                   // backtracked quantities (after DP)
    vector<vector<std::shared_ptr<Insertion>>> breakpoints; //
    vector<double> cost;                                  // cost[scenario] : DP total cost for this scenario

    // fill _breakpoints with every breakpoint of piecewise linear function 
    // (i.e. points where the slope changes)
    // also fills repeat with points x such that [x, x] is a linear piece of f
    void extractBreakpointsAndRepeat(const std::shared_ptr<PLFunction>& f, vector<double>& _breakpoints, vector<double>& repeat);

    // first combine breakpoints of f1 and f2 (without repetition)
    // then add repeats of f1 and f2 because they are particular points
    vector<double > getBreakpoints_final(std::shared_ptr<PLFunction> f1, std::shared_ptr<PLFunction> f2);

    // min(f1, f2) with f1 and f2 two piecewise linear functions
    std::shared_ptr<PLFunction> min_final(std::shared_ptr<PLFunction> f1, std::shared_ptr<PLFunction> f2);

    // create a copy of chosenPiece (fromC, fromF, etc.) but with (x1, y1) and (x2, y2) points
    std::shared_ptr<LinearPiece> createPieceFromLowerY(std::shared_ptr<LinearPiece> chosenPiece,
        double x1, double y1, double x2, double y2);
   
    // min_q(fromC(x + q) + fromF(q))
    std::shared_ptr<PLFunction> supperposition(std::shared_ptr<PLFunction> fromC, std::shared_ptr<PLFunction> fromF);
    
    // previous supperposition but for 2 pieces (used in global supperposition)
    std::shared_ptr<PLFunction> supperpositionPieces(std::shared_ptr<LinearPiece> fromPieceC, std::shared_ptr<LinearPiece> fromPieceF);
    
    // C[horizon] = 0 (start for backward DP)
    void lastDay(std::shared_ptr<PLFunction> Cf, unsigned int scenario);
    
    // DP + backtrack + fill local search structures
    bool solveStockoutBackward();
    
    // backward DP for a given scenario and for day >= 2 (and fill associated functions C)
    void solveOneScenario(unsigned int scenario);

    // first day is a specific case (route of day 1 is shared by scenarios)
    void day1(unsigned int scenario);

    // equation system to retrieve DP solution with backtracking
    void solveEquationSystemHoldingBackward(std::shared_ptr<LinearPiece> C, std::shared_ptr<LinearPiece> fromC,
        std::shared_ptr<LinearPiece> fromF, double IAtT, double demand, double &nextI, double &quantity);
    
    // backtrack to retrieve the corresponding solution
    bool backtrack(unsigned int scenario, unsigned int idxInsert);

    // constructor
    LotSizingSolver(Params* params, vector<vector<vector<Insertion>>> insertions, unsigned int _client);
    
    // destructor
    ~LotSizingSolver();
};


#endif 
