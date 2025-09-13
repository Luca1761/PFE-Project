//
// Created by pta on 15/08/2018.
//

#include "LinearPiece.h"
#include <algorithm>
#include <iomanip>
LinearPiece::LinearPiece() {
    p1 = make_shared<Point>();
    p2 = make_shared<Point>();

    fromC = NULL;
    fromF = NULL;
    fromC_pre = NULL;
    fromInst = NULL;
    stockout = false;
}

LinearPiece::LinearPiece(double left_x, double left_y, double right_x, double right_y) {
    p1 = make_shared<Point>(left_x, left_y);
    p2 = make_shared<Point>(right_x, right_y);
    updateSlope();

    stockout = false;
    fromC = NULL;
    fromC_pre = NULL;
    fromF = NULL;
    fromInst = NULL;
}

LinearPiece::LinearPiece(LinearPiece *lp) {
    p1 = make_shared<Point>(lp->p1->x, lp->p1->y);
    p2 = make_shared<Point>(lp->p2->x, lp->p2->y);
    updateSlope();
    
    slope = lp->slope;
    fromC = lp->fromC;
    fromC_pre = lp->fromC_pre;
    fromF = lp->fromF;
    fromInst = lp->fromInst;
    stockout = lp->stockout;
}

void LinearPiece::updateSlope() {
    if(eq(p1->x, p2->x)) slope = 0.0;
    else slope = (p2->y - p1->y) / (p2->x - p1->x);
}

void LinearPiece::update(double left_x, double left_y, double right_x, double right_y) {
    p1 = make_shared<Point>(left_x, left_y);
    p2 = make_shared<Point>(right_x, right_y);
    updateSlope();
}

double LinearPiece::cost(double x) {
    if (lt(x, p1->x) || gt(x, p2->x)) {
        std::cout << "Linear piece is not defined for this value!!" << std::endl; 
        throw std::string("Linear piece not defined for this value!!");
    }
    if(eq(p1->x, p2->x)){
        return min(p1->y, p2->y);
    }
    updateSlope();
    return slope * (x - p2->x) + p2->y;
}

std::shared_ptr<LinearPiece> LinearPiece::getInBound(double lb, double ub) const{
    if (gt(lb, ub)) return nullptr;
    if (lt(p2->x, lb) || gt(p1->x, ub)) return nullptr; 
    if (le(lb, p1->x) && le(p2->x, ub)) return std::make_shared<LinearPiece>(*this);
    
    double startX = std::max<double>(p1->x, lb);
    double startY = (startX - p1->x) * slope + p1->y;
    
    double endX = std::min<double>(p2->x, ub);
    double endY = (endX - p1->x) * slope + p1->y;
    
    std::shared_ptr<LinearPiece> lp = std::make_shared<LinearPiece>(*this);
    lp->update(startX, startY, endX, endY);    
    return lp;
}

bool LinearPiece::eqlp(const shared_ptr<LinearPiece> rhs) {
    return rhs != NULL && eq(p1->x, rhs->p1->x) && eq(p1->y, rhs->p1->y) && eq(p2->x, rhs->p2->x) && eq(p2->y, rhs->p2->y);
}

bool LinearPiece::eqFromC(const shared_ptr<LinearPiece> rhs) {
    if (rhs == NULL) return false;
    else if (fromC == NULL && rhs->fromC == NULL) return true;
    else if (fromC == NULL || rhs->fromC == NULL) return false;

    return fromC->eqlp(rhs->fromC);
}

bool LinearPiece::eqFromCpre(const shared_ptr<LinearPiece> rhs) {
    if (rhs == NULL) return false;
    else if (fromC_pre == NULL && rhs->fromC_pre == NULL) return true;
    else if (fromC_pre == NULL || rhs->fromC_pre == NULL) return false;

    return fromC_pre->eqlp(rhs->fromC_pre);
}

bool LinearPiece::eqFromF(const shared_ptr<LinearPiece> rhs) {
    if (rhs == NULL) return false;
    else if (fromF == NULL && rhs->fromF == NULL) return true;
    else if (fromF == NULL || rhs->fromF == NULL) return false;

    return fromF->eqlp(rhs->fromF);
}

bool LinearPiece::eqDeep(const shared_ptr<LinearPiece> rhs) {
    return eqlp(rhs) && eqFromC(rhs) && eqFromF(rhs) && stockout == rhs->stockout;
}

std::shared_ptr<LinearPiece> LinearPiece::clone() {
    return make_shared<LinearPiece>(this);
}

void LinearPiece::print() {   
    cout << "(" << p1->x << ", " << p1->y << ", " << p2->x
    << ", " << p2->y << ") ";
    cout << endl;
}

LinearPiece::~LinearPiece() {}
