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
    Params * params;
     
public:
    vector<vector<Insertion>> insertions;
    vector<std::shared_ptr<PLFunction>> C;
    vector<std::shared_ptr<PLFunction>> F;
    vector<std::shared_ptr<PLFunction>> F1;
    vector<std::shared_ptr<PLFunction>> F2;
    vector<double> I;
    vector<double> quantities;
    vector<std::shared_ptr<Insertion>> breakpoints;

    int client;
    int horizon;
    double objective;

    LotSizingSolver(Params * params, vector<vector<Insertion>> insertions, int client);

    std::shared_ptr<PLFunction> copyPLFunction(std::shared_ptr<PLFunction> source);
    void extractBreakpoints(const std::shared_ptr<PLFunction>& f, vector<double>& breakpoints);
    void extractRepeat(const std::shared_ptr<PLFunction>& f, vector<double>& breakpoints);
    vector<double > getBreakpoints_final(std::shared_ptr<PLFunction> f1, std::shared_ptr<PLFunction> f2);
    std::shared_ptr<PLFunction> min_final(std::shared_ptr<PLFunction> f1, std::shared_ptr<PLFunction> f2);
    std::shared_ptr<LinearPiece> createPieceFromLowerY(std::shared_ptr<LinearPiece> &chosenPiece,
         double x1, double y1, double x2, double y2);
    std::shared_ptr<PLFunction> supperpositionPieces(std::shared_ptr<LinearPiece> fromPieceC, std::shared_ptr<LinearPiece> fromPieceF);
    std::shared_ptr<PLFunction> supperposition(std::shared_ptr<PLFunction> fromC, std::shared_ptr<PLFunction> fromF);
    
    void Lastday(vector<std::shared_ptr<PLFunction>> &C);
    bool solveStockoutBackward();
    void solveEquationSystemHoldingBackward(std::shared_ptr<LinearPiece> C, std::shared_ptr<LinearPiece> fromC,
        std::shared_ptr<LinearPiece> fromF, double I, double demand, double &fromI, double &quantity,std::shared_ptr<LinearPiece> Cpre,double inventory);
    bool backtrackingStockoutBackward();
    ~LotSizingSolver();
};


#endif //IRP_LOTSIZINGSOLVER_H
