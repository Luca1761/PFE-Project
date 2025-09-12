//
// Created by pta on 15/08/2018.
//

#ifndef IRP_LINEARPIECE_H
#define IRP_LINEARPIECE_H

#include "memory"
#include <algorithm>
#include <string>
#include <vector>
#include "Params.h"

using namespace std;

constexpr double EPSILON = 1e-5;

inline bool le(const double &x, const double &y) { return x <= y + EPSILON; } // less than or equal to

inline bool lt(const double &x, const double &y) { return x + EPSILON < y; } // less than

inline bool eq(const double &x, const double &y) { return fabs(x - y) < EPSILON; } // equal

inline bool neq(const double &x, const double &y) { return !eq(x, y); } // not equal

inline bool gt(const double &x, const double &y) { return lt(y, x); } // greater than

inline bool ge(const double &x, const double &y) { return le(y, x); } // greater than or equal to

class Params;

struct Point
{
    double x;
    double y;

    Point() {}

    Point(double _x, double _y) : x(_x), y(_y) {}

    virtual ~Point() {}
    
    friend bool operator==(const Point& lhs, const Point& rhs) {
        return eq(lhs.x, rhs.x) && eq(lhs.y, rhs.y);
    }

    inline Point convolve(shared_ptr<Point> p) {
        return Point(x + p->x, y + p->y);
    }
};

class LinearPiece
{

public:
    double slope = 0.;
    shared_ptr<Point> p1;   // starting point of linear piece
    shared_ptr<Point> p2;   // end point of linear piece
    
    shared_ptr<LinearPiece> fromC_pre;      // original C function used to get that linear piece
    shared_ptr<LinearPiece> fromC;  // translated and modified C function to get that linear piece
    shared_ptr<LinearPiece> fromF;      // F function used to get that linear piece
    shared_ptr<Insertion> fromInst;     // insertion used to get that linear piece

    // specify is this linear piece comes from a stockout
    bool stockout;

    // constructors
    LinearPiece();
    LinearPiece(LinearPiece *lp);
    LinearPiece(double left_x, double left_y, double right_x, double right_y);

    // update the slope if values changed
    void updateSlope();
    
    // update linear piece with new values (also update the slope)
    void update(double left_x, double left_y, double right_x, double right_y);
    
    // return value of linear piece applied on x
    double cost(double x);
    
    // get linear piece in [lb, ub]
    std::shared_ptr<LinearPiece> getInBound(double lb, double ub) const;
    
    // check equality of linear pieces (without fromC, fromCPre and fromF)
    bool eqlp(const shared_ptr<LinearPiece> rhs);

    // check equality of fromC
    bool eqFromC(const shared_ptr<LinearPiece> rhs);

    // check equality of fromCPre
    bool eqFromCpre(const shared_ptr<LinearPiece> rhs);

    // check equality of fromF
    bool eqFromF(const shared_ptr<LinearPiece> rhs);
    
    // check equality of linear pieces (with fromC, fromCPre and fromF)
    bool eqDeep(const shared_ptr<LinearPiece> rhs);
    
    // copy the linear piece
    std::shared_ptr<LinearPiece> clone();
    
    // display linear piece
    void print();

    // destructor
    ~LinearPiece();
};

#endif // IRP_LINEARPIECE_H
