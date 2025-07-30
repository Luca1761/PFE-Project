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
    Params* params;
     
public:
    vector<vector<vector<Insertion>>> insertions;
    vector<vector<std::shared_ptr<PLFunction>>> C;
    vector<vector<double>> Ck;
    vector<vector<std::shared_ptr<PLFunction>>> savedCk; 
    vector<vector<std::shared_ptr<PLFunction>>> Gk;
    vector<vector<std::shared_ptr<PLFunction>>> F;
    vector<vector<double>> I;
    vector<vector<double>> quantities;
    vector<vector<std::shared_ptr<Insertion>>> breakpoints;

    int client;
    int horizon;
    vector<double> objective;
    int nbScenario;

    LotSizingSolver(Params* params, vector<vector<vector<Insertion>>> insertions, int client);

    std::shared_ptr<PLFunction> copyPLFunction(std::shared_ptr<PLFunction> source, Params* paramsTemp);
    void extractBreakpoints(const std::shared_ptr<PLFunction>& f, vector<double>& breakpoints);
    void extractRepeat(const std::shared_ptr<PLFunction>& f, vector<double>& breakpoints);
    vector<double > getBreakpoints_final(std::shared_ptr<PLFunction> f1, std::shared_ptr<PLFunction> f2);
    std::shared_ptr<PLFunction> min_final(std::shared_ptr<PLFunction> f1, std::shared_ptr<PLFunction> f2, Params* paramsTemp);
    std::shared_ptr<LinearPiece> createPieceFromLowerY(std::shared_ptr<LinearPiece> chosenPiece,
         double x1, double y1, double x2, double y2);
    std::shared_ptr<PLFunction> supperpositionPieces(std::shared_ptr<LinearPiece> fromPieceC, std::shared_ptr<LinearPiece> fromPieceF, Params* paramsTemp);
    std::shared_ptr<PLFunction> supperposition(std::shared_ptr<PLFunction> fromC, std::shared_ptr<PLFunction> fromF, Params* paramsTemp);
    
    void Lastday(vector<std::shared_ptr<PLFunction>> &C, int scenario);
    bool solveStockoutBackward();
    void day1(int scenario);
    void solveOneScenario(int scenario);
    void solveEquationSystemHoldingBackward(std::shared_ptr<LinearPiece> C, std::shared_ptr<LinearPiece> fromC,
        std::shared_ptr<LinearPiece> fromF, double I, double demand, double &fromI, double &quantity);
    bool backtrackingStockoutBackward(unsigned int scenario, int valInsertMin);
    ~LotSizingSolver();
};


#endif //IRP_LOTSIZINGSOLVER_H
