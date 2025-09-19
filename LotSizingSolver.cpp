//
// //
// Created by zjy on 13/05/2024 (SOA).
//
//

#include "LotSizingSolver.h"
#include <memory>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iomanip>
#include <thread>

LotSizingSolver::LotSizingSolver(Params* _params, vector<vector<vector<Insertion>>> inst, unsigned int _client)
    : params(_params), insertions(inst), client(_client), nbScenario(params->nbScenarios), horizon(params->nbDays) {

  unsigned int nbInsertionsFirstDay = (unsigned int) insertions[0][0].size();

  F = vector<vector<std::shared_ptr<PLFunction>>>(nbScenario, vector<std::shared_ptr<PLFunction>>(horizon));
  Fk = vector<vector<std::shared_ptr<PLFunction>>>(nbScenario);
  Ck = vector<vector<double>>(nbScenario, vector<double>(nbInsertionsFirstDay + 1,10000000));
  savedCk = vector<vector<std::shared_ptr<PLFunction>>>(nbScenario, vector<std::shared_ptr<PLFunction>>(nbInsertionsFirstDay + 1));
 
  // compute structures for day 1
  for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
    if (nbInsertionsFirstDay != insertions[scenario][0].size()) {
      std::cout << "INSERTIONS SHOULD BE THE SAME FOR DAY 1" << std::endl;
      throw std::string("INSERTIONS SHOULD BE THE SAME FOR DAY 1");
    }
    
    for (unsigned int i = 0; i < nbInsertionsFirstDay; i++) {
      std::shared_ptr<PLFunction> reverseFk = std::make_shared<PLFunction>(params, insertions[scenario][0][i], 
          params->cli[client], scenario);
          
      Fk[scenario].push_back(std::make_shared<PLFunction>(params));
      //we need to "reverse" in a negative sense
      // so we can use C(x+q) + F(q) to compute (C(x-q) + F(q))
      for (auto &piece : reverseFk->pieces) {
        std::shared_ptr<LinearPiece> tmp = piece->clone();
        tmp->update(-tmp->p2->x, tmp->p2->y, -tmp->p1->x, tmp->p1->y);
        Fk[scenario][i]->append(tmp);
      }
    }
  }

  // compute structures for day >= 2
  for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
    for (unsigned int t = 1; t < horizon; t++) {
      std::shared_ptr<PLFunction> reverseF = std::make_shared<PLFunction>(params, insertions[scenario][t], t, params->cli[client], scenario); 
      
      F[scenario][t] = std::make_shared<PLFunction>(params);
      // same reverse
      for (unsigned int i = reverseF->nbPieces; i >= 1 ; i--) {
        std::shared_ptr<LinearPiece> tmp = reverseF->pieces[i - 1]->clone();
        tmp->update(-tmp->p2->x, tmp->p2->y, -tmp->p1->x, tmp->p1->y);
        F[scenario][t]->append(tmp);
      }
    }
  }
}

void LotSizingSolver::extractBreakpointsAndRepeat(const std::shared_ptr<PLFunction> &f, vector<double> &_breakpoints, vector<double> &repeat) {
  for (unsigned int i = 0; i < f->nbPieces; i++) {
      if(fabs(f->pieces[i]->p1->x - f->pieces[i]->p2->x) < 0.00001) {
        repeat.push_back(f->pieces[i]->p1->x);
      } else if((i == 0 || neq(f->pieces[i]->p1->x, f->pieces[i - 1]->p2->x))) {
        _breakpoints.push_back(f->pieces[i]->p1->x);
      }
      _breakpoints.push_back(f->pieces[i]->p2->x);
  }
}

vector<double> LotSizingSolver::getBreakpoints_final(std::shared_ptr<PLFunction> f1, std::shared_ptr<PLFunction> f2){
  vector<double> _breakpoints;
  vector<double> repeat;
  extractBreakpointsAndRepeat(f1, _breakpoints, repeat);
  extractBreakpointsAndRepeat(f2, _breakpoints, repeat);
  std::sort(_breakpoints.begin(), _breakpoints.end());
  std::sort(repeat.begin(), repeat.end());
  auto it = std::unique(_breakpoints.begin(), _breakpoints.end(), 
                      [](double a, double b) { return fabs(a - b) < 0.00001; });
  _breakpoints.erase(it, _breakpoints.end()); // erase repetitions

  auto it1 = std::unique(repeat.begin(), repeat.end(), 
                        [](double a, double b) { return fabs(a - b) < 0.00001; });
  repeat.erase(it1, repeat.end()); // erase repetitions

  _breakpoints.insert(_breakpoints.end(), repeat.begin(), repeat.end()); // add particular points
  
  std::sort(_breakpoints.begin(), _breakpoints.end());
  return _breakpoints;
}

std::shared_ptr<LinearPiece> LotSizingSolver::createPieceFromLowerY(
  std::shared_ptr<LinearPiece> chosenPiece, double x1, double y1, double x2, double y2)  {
    std::shared_ptr<LinearPiece> tmpPiece = chosenPiece->clone();
    tmpPiece->update(x1, y1, x2, y2);
    return tmpPiece;
}

std::shared_ptr<PLFunction> LotSizingSolver::min_final(std::shared_ptr<PLFunction> f1, 
                                                      std::shared_ptr<PLFunction> f2) {
  std::shared_ptr<PLFunction> f = std::make_shared<PLFunction>(params);
  if (f1->nbPieces == 0 && f2->nbPieces == 0) return f;
  if (f1->nbPieces == 0 && f2->nbPieces > 0) return f2->clone();
  if (f1->nbPieces > 0 && f2->nbPieces == 0) return f1->clone();

  vector<double> _breakpoints = getBreakpoints_final(f1, f2);
  std::shared_ptr<LinearPiece> lp1 = f1->pieces[0];
  std::shared_ptr<LinearPiece> lp2 = f2->pieces[0];
  unsigned int n1 = 0, n2 = 0;
  for (unsigned int i = 1; i < _breakpoints.size(); i++) {
    if (lp1) lp1->updateSlope();
    if (lp2) lp2->updateSlope();
    
    std::shared_ptr<LinearPiece> piece1 = nullptr;
    if (lp1) piece1 = lp1->getInBound(_breakpoints[i - 1], _breakpoints[i]); 
    
    std::shared_ptr<LinearPiece> piece2 = nullptr;
    if (lp2) piece2 = lp2->getInBound(_breakpoints[i - 1], _breakpoints[i]); 
    
    if (!piece1 && !piece2);
    else if (!piece1) f->append(piece2);
    else if (!piece2) f->append(piece1);
    else {
      double x, y;
      bool intersects = f->intersect(piece1, piece2, x, y);
      if (eq(piece2->p1->x, piece2->p2->x) && eq(piece1->p2->x, piece2->p1->x)) {
        f->append(piece1);
        if (gt(piece1->p2->y, piece2->p1->y)) 
          f->append(piece2);
      } else if (eq(piece1->p1->x, piece1->p2->x) && eq(piece2->p2->x, piece1->p1->x)) {
        f->append(piece2);
        if (gt(piece2->p2->y, piece1->p1->y)) 
          f->append(piece1);
      } else if (intersects) {
        std::shared_ptr<LinearPiece> firstPiece, secondPiece, chosenPiece1, chosenPiece2;
        chosenPiece1 = (piece1->p1->y <= piece2->p1->y) ? piece1 : piece2;
        if (round(x) == chosenPiece1->p2->x) {
          firstPiece = createPieceFromLowerY(chosenPiece1->clone(), chosenPiece1->p1->x, chosenPiece1->p1->y, chosenPiece1->p2->x, chosenPiece1->p2->y);
          f->append(firstPiece);
          chosenPiece2 = ((piece1->p2->y < piece2->p2->y) )? piece1 : piece2;

          secondPiece = createPieceFromLowerY(chosenPiece2->clone(), chosenPiece2->p2->x, chosenPiece2->p2->y, chosenPiece2->p2->x, chosenPiece2->p2->y);
          f->append(secondPiece);
        } else {
          firstPiece = createPieceFromLowerY(chosenPiece1->clone(), chosenPiece1->p1->x, chosenPiece1->p1->y, x, y);
          f->append(firstPiece);

          chosenPiece2 =((piece1->p2->y < piece2->p2->y) )? piece1 : piece2;
          secondPiece = createPieceFromLowerY(chosenPiece2->clone(), x, y, chosenPiece2->p2->x, chosenPiece2->p2->y);
          f->append(secondPiece);
        }
      } else if (eq(piece1->p1->x ,x) ){
          if(piece1->p2->y < piece2->p2->y) 
              f->append(piece1);
            else 
              f->append(piece2);
      }
      else if (eq(piece1->p2->x ,x) ){
          if(piece1->p1->y <piece2->p1->y) 
              f->append(piece1);
            else 
              f->append(piece2);
      }
      else {
        if ((piece1->p1->y -piece2->p1->y) + (piece1->p2->y -piece2->p2->y) < 0) 
              f->append(piece1);
          else 
              f->append(piece2);
          
      }
      
    }
    if (lp1&& fabs(_breakpoints[i]- lp1->p2->x)<0.00001 ){
        if(n1 + 1 < f1->nbPieces) lp1 = f1->pieces[++n1];
        else lp1 = nullptr;
    }
    if (lp2&& fabs(_breakpoints[i]- lp2->p2->x)<0.00001 ){
      if(n2 + 1 < f2->nbPieces) lp2 = f2->pieces[++n2];
      else lp2 = nullptr;
    }
    if (!lp1 && !lp2) break;
  }
  return f;
}

std::shared_ptr<PLFunction> LotSizingSolver::supperpositionPieces(std::shared_ptr<LinearPiece> fromPieceC, std::shared_ptr<LinearPiece> fromPieceF) {
  std::shared_ptr<PLFunction> f(std::make_shared<PLFunction>(params));
  vector<Point> points;
  shared_ptr<LinearPiece> tmpPiece, newpiece;
  points.push_back(fromPieceC->p1->convolve(fromPieceF->p2));  
  points.push_back(fromPieceC->p1->convolve(fromPieceF->p1));  
  points.push_back(fromPieceC->p2->convolve(fromPieceF->p1));
  points.push_back(fromPieceC->p2->convolve(fromPieceF->p2));
   
  std::vector<Point> sortedPoints = points; 
  std::sort(sortedPoints.begin(), sortedPoints.end(), 
    [](const Point& a, const Point& b) { return lt(a.y, b.y); });
  Point minYPoint1 = sortedPoints[0];
  Point minYPoint2 = sortedPoints[1]; Point minYPoint3 = sortedPoints[2];Point minYPoint4 = sortedPoints[3];

  if(eq(minYPoint1.x,minYPoint2.x) && eq(minYPoint1.y,minYPoint2.y)){
    minYPoint2 = sortedPoints[2];
    if(gt(minYPoint1.x,minYPoint2.x)) std::swap(minYPoint1,minYPoint2);

    tmpPiece = make_shared<LinearPiece>(minYPoint1.x, minYPoint1.y, minYPoint2.x, minYPoint2.y);
    f->append(tmpPiece); 
    f->pieces[0]->fromC = fromPieceC->clone();
    f->pieces[0]->fromF = fromPieceF->clone();
    f->pieces[0]->fromInst = fromPieceF->fromInst; 
    return f;
  }

  unsigned int minYIndex1 = (unsigned int) (std::find(points.begin(), points.end(), minYPoint1) - points.begin());
  unsigned int minYIndex2 = (unsigned int) (std::find(points.begin(), points.end(), minYPoint2) - points.begin());
  unsigned int next1 = (minYIndex1 + 1) % 4,next2 = (minYIndex1+3)%4;

  double slope1 = (eq(points[next1].x, points[minYIndex1].x) || eq(points[next1].y, points[minYIndex1].y))
                    ? 0 : (points[next1].y- points[minYIndex1].y )/(points[next1].x- points[minYIndex1].x);
  double slope2 = (eq(points[next2].x, points[minYIndex1].x) || eq(points[next2].y, points[minYIndex1].y))
                    ? 0 : (points[next2].y- points[minYIndex1].y )/(points[next2].x- points[minYIndex1].x);
  if(eq(slope1, 0) || eq(slope2, 0)){//one 0
    if(minYPoint1.x>minYPoint2.x )std::swap(minYIndex1,minYIndex2);
    next1 = (minYIndex1 + 1) % 4; next2 = (minYIndex1+3)%4;
    unsigned int rightup = 6-minYIndex1-next1-next2;
    unsigned int leftup = 6-minYIndex1-minYIndex2-rightup;

    if(lt(points[rightup].x,points[minYIndex2].x)){
        tmpPiece = make_shared<LinearPiece>(points[leftup].x, points[leftup].y,points[minYIndex1].x, points[minYIndex1].y); 
        newpiece = make_shared<LinearPiece>(points[minYIndex1].x, points[minYIndex1].y, points[minYIndex2].x, points[minYIndex2].y);
    } else{
        tmpPiece = make_shared<LinearPiece>(points[minYIndex1].x, points[minYIndex1].y, points[minYIndex2].x, points[minYIndex2].y); 
        newpiece = make_shared<LinearPiece>(points[minYIndex2].x, points[minYIndex2].y,points[rightup].x, points[rightup].y);
    }
  }  
  else if ((lt(slope1, 0.0) && lt(0.0, slope2)) || (lt(slope2, 0.0) && lt(0.0, slope1))) {//one negative one postive
    if(points[next1].x>points[next2].x)std::swap(next1,next2);
    tmpPiece = make_shared<LinearPiece>(points[next1].x, points[next1].y,points[minYIndex1].x, points[minYIndex1].y); 
    newpiece = make_shared<LinearPiece>(points[minYIndex1].x, points[minYIndex1].y,  points[next2].x, points[next2].y);
    
  }
  else if(lt(slope1, 0.0)){//both negative
    if(eq(minYPoint2.x,minYPoint3.x) && eq(minYPoint2.y,minYPoint3.y) ){
      tmpPiece = make_shared<LinearPiece>(minYPoint4.x, minYPoint4.y,minYPoint2.x, minYPoint2.y); 
      newpiece = make_shared<LinearPiece>(minYPoint2.x, minYPoint2.y,minYPoint1.x, minYPoint1.y);  
    } else{
      unsigned int tmp1 = gt(slope1,slope2) ? next1 : next2;
      unsigned int tmp2 = 6 - minYIndex1 - next1 - next2;
      tmpPiece = make_shared<LinearPiece>(points[tmp2].x, points[tmp2].y,points[tmp1].x, points[tmp1].y); 
      newpiece = make_shared<LinearPiece>(points[tmp1].x, points[tmp1].y,points[minYIndex1].x, points[minYIndex1].y);  
    }
  } else{//both positive
    if(eq(minYPoint2.x,minYPoint3.x) && eq(minYPoint2.y,minYPoint3.y) ){//line
      tmpPiece = make_shared<LinearPiece>(minYPoint1.x, minYPoint1.y,minYPoint2.x, minYPoint2.y); 
      newpiece = make_shared<LinearPiece>(minYPoint2.x, minYPoint2.y,minYPoint4.x, minYPoint4.y);  
    } else{
      unsigned int tmp1 = gt(slope1,slope2) ? next2 : next1;
      unsigned int tmp2 = 6 - minYIndex1 - next1 - next2;
     
      tmpPiece = make_shared<LinearPiece>(points[minYIndex1].x, points[minYIndex1].y,points[tmp1].x, points[tmp1].y); 
      newpiece = make_shared<LinearPiece>(points[tmp1].x, points[tmp1].y,points[tmp2].x, points[tmp2].y);  
    }
  }
  f->append(tmpPiece);
  f->append(newpiece);
  for(auto &piece : f->pieces){
    piece->fromC = fromPieceC->clone();
    piece->fromF = fromPieceF->clone();
    piece->fromInst = fromPieceF->fromInst; 
  }
  return f;
}

std::shared_ptr<PLFunction> LotSizingSolver:: supperposition(std::shared_ptr<PLFunction> fromC, std::shared_ptr<PLFunction> fromF) {
 std::shared_ptr<PLFunction> f(std::make_shared<PLFunction>(params));

  for (auto& pieceC : fromC->pieces){
    for (auto& pieceF : fromF->pieces){
      
      std::shared_ptr<PLFunction> tmpF = std::make_shared<PLFunction>(params);
      std::shared_ptr<LinearPiece> fromPieceC = pieceC->clone();
      
      std::shared_ptr<LinearPiece> fromPieceF = pieceF->clone();

      tmpF = supperpositionPieces(fromPieceC, fromPieceF);
      
      if (f->nbPieces > 0 || tmpF->nbPieces > 0){
        std::shared_ptr<PLFunction> fmin;
        fmin = min_final(f, tmpF);
        f = fmin->clone();
      }

      tmpF.reset();
      fromPieceC.reset();
      fromPieceF.reset();
    }
  }
  return f;
}

void LotSizingSolver::solveEquationSystemHoldingBackward(std::shared_ptr<LinearPiece> Ct,
                                          std::shared_ptr<LinearPiece> fromC,
                                          std::shared_ptr<LinearPiece> fromF,
                                          double IAtT, double demand,
                                          double &nextI, double &quantity){
  if(eq(fromC->p1->x,fromC->p2->x)){
    nextI = std::max(0.0, fromC->p2->x - demand);
    quantity = nextI + demand - IAtT;
    return;
  }
  if (eq(fromF->p1->x, fromF->p2->x)) {
    quantity = -fromF->p1->x;
    nextI = IAtT - demand + quantity;
    return;
  }

  double slopeC = (fromC->p2->y - fromC->p1->y) / (fromC->p2->x - fromC->p1->x);
  double slopeF = -(fromF->p2->y - fromF->p1->y) / (fromF->p2->x - fromF->p1->x);
  
  if(eq(slopeC, -slopeF)){
    double upperbound = std::min<double>(fromC->p2->x - IAtT, -fromF->p1->x);
    double lowerbound = std::max<double>(fromC->p1->x - IAtT, -fromF->p2->x);
    
    if (ge(upperbound, lowerbound) ){
        quantity = upperbound;
        nextI = IAtT - demand + quantity;
    } else{
      std::cout << "Lowerbound is greater than upperbound" << std::endl;
      throw std::string("Lowerbound is greater than upperbound");
    }
    return;
  }
  slopeC*=10000;slopeF*=10000;
  double x1 = fromC->p2->x, y1 = fromC->p2->y,x2 = -fromF->p1->x,y2 = fromF->p1->y;
  double numerator = Ct->cost(std::max<double>(0, IAtT)) * 10000- y1 * 10000 - y2 * 10000;
    numerator -= slopeC * (IAtT-x1);
    numerator += slopeF * x2;
    quantity = numerator / (slopeF + slopeC);
    double left = std::max<double>(fromC->p1->x - IAtT,-fromF->p2->x);
    double right = std::min<double>(fromC->p2->x - IAtT,-fromF->p1->x);
    if(gt(quantity ,right )||lt(quantity,left)){
      if (gt(quantity ,right)) { quantity = right;}
      if (lt(quantity,left)) {quantity = left;}
    }
    
  nextI = IAtT - demand + quantity;
}

bool LotSizingSolver::backtrack(unsigned int scenario, unsigned int idxInsert) {
  // initialization
  for (unsigned int t = 0; t < horizon; t++){
    quantities[scenario][t] = 0.0;
    breakpoints[scenario][t] = nullptr;
    I[scenario][t] = 0;
  }
  I[scenario].push_back(0.0);
  unsigned int day = 0;
  
  if (savedCk[scenario][idxInsert]->nbPieces == 0) {
    std::cout << "ERROR, at least one piece is needed" << std::endl;
    throw std::string("ERROR, at least one piece is needed");
  }
  std::shared_ptr<LinearPiece> tmp = 
      savedCk[scenario][idxInsert]->getMinimalPiece(I[scenario][day], cost[scenario]);

  if (I[scenario][day] != params->cli[client].startingInventory) throw std::string("Wrong initial inventory");
  ///////////////////////////////////////////
  while (tmp != nullptr) {

    // if do not delivery any thing, then inventory at the end of next day
    // equals this day demand
    if (!tmp->stockout && !tmp->fromF ) { //f1 q = 0
      I[scenario][day + 1] = I[scenario][day] - params->cli[client].dailyDemand[scenario][day+1];
    } else if (tmp->stockout && !tmp->fromF) {//f2  q= 0, not enough
      I[scenario][day + 1] = 0.0;
    } else if (!tmp->stockout && tmp->fromF) { //f3
      std::shared_ptr<LinearPiece> fromC = tmp->fromC->clone();
      std::shared_ptr<LinearPiece> fromF = tmp->fromF->clone();
      solveEquationSystemHoldingBackward(tmp, fromC, fromF, I[scenario][day],
                      params->cli[client].dailyDemand[scenario][day + 1], I[scenario][day + 1],
                      quantities[scenario][day]);
        
      shared_ptr<LinearPiece> tmpF = tmp->fromF;
      breakpoints[scenario][day] = tmpF->fromInst;
    } else if (tmp->stockout && tmp->fromF ) {//f4
      std::shared_ptr<LinearPiece> fromC = tmp->fromC->clone();
      std::shared_ptr<LinearPiece> fromF = tmp->fromF->clone();
      solveEquationSystemHoldingBackward(tmp, fromC, fromF, I[scenario][day],
                      params->cli[client].dailyDemand[scenario][day + 1], I[scenario][day + 1],
                      quantities[scenario][day]);
      
      I[scenario][day + 1] = 0.0;

      shared_ptr<LinearPiece> tmpF = tmp->fromF;
      breakpoints[scenario][day] = tmpF->fromInst;
    }
    tmp = tmp->fromC_pre;

    day = day + 1;
    if(day == horizon) break;
  }
  return true;
}

void LotSizingSolver::lastDay(std::shared_ptr<PLFunction> Cf, unsigned int scenario){//for the last day
  double maxInventory = (params->endDayInventories) ? params->cli[client].maxInventory 
                                                    : std::max(0.0, params->cli[client].maxInventory - params->cli[client].dailyDemand[scenario][horizon]);
  std::shared_ptr<LinearPiece> tmp = std::make_shared<LinearPiece>(0, 0, maxInventory, 0);
  Cf->append(tmp);
}

bool LotSizingSolver::solveStockoutBackward() {
  C = vector<vector<std::shared_ptr<PLFunction>>>(nbScenario, vector<std::shared_ptr<PLFunction>>(horizon + 1));

  // final inventory at the end of each period
  I = vector<vector<double>>(nbScenario, vector<double>(horizon));

  // replinishment at each period
  quantities = vector<vector<double>>(nbScenario, vector<double>(horizon));
  cost = vector<double>(nbScenario, 0.0);

  breakpoints = vector<vector<std::shared_ptr<Insertion>>>(nbScenario, vector<std::shared_ptr<Insertion>>(horizon));
  
  for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
    for (unsigned int i = 0; i < horizon + 1; i++) {
      C[scenario][i] = std::make_shared<PLFunction>(params); // a piecewise linear function in each period
    }
  }
  const unsigned int GROUP_SIZE_2 = nbScenario / params->nbCores + 1;
	vector<thread> threads;
	for (unsigned int scenario = 0; scenario < nbScenario; scenario += GROUP_SIZE_2) {
		unsigned int end = std::min(scenario + GROUP_SIZE_2, nbScenario);
		threads.emplace_back([this, scenario, end]() {
			for (unsigned int scenario1 = scenario; scenario1 < end; scenario1++) {
				solveOneScenario(scenario1);
        day1(scenario1);
			}
		});
	}
  for (auto& t : threads) {
		t.join();
	}
  
  unsigned int idxInsert = 0;
  double valInsertMin = __DBL_MAX__;
  for (unsigned int k = 0; k < Ck[0].size(); k++) { // loop across all the possible insertions for day 1
    double valInsert = 0;
    for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
      valInsert += Ck[scenario][k]; // best value for insertion k on day1 on this scenario
    }
    valInsert /= (double) nbScenario;
    if (valInsert < valInsertMin) { // we take the min value
      valInsertMin = valInsert;
      idxInsert = k;
    }
  }
  for (unsigned int scenario = 0; scenario < nbScenario; scenario++) {
    backtrack(scenario, idxInsert); // backtrack quantities and inventories for every scenario, for the best insertion index
  }
  return true;
}

void LotSizingSolver::solveOneScenario(unsigned int scenario) {
    double maxInventory;
    double maxDeliverable;
    lastDay(C[scenario][horizon], scenario);

    for (unsigned int t = horizon; t >= 2; t--) {
      if (neq(C[scenario][t]->pieces[0]->p1->x, 0)) {
        std::cout << "Ct(0) doesn't exist!!!: " << C[scenario][t]->pieces[0]->p1->x << std::endl;
        throw string("Ct(0) doesn't exist!!!");
      }
      //piecewise linear functions
      std::shared_ptr<PLFunction> f1;
      std::shared_ptr<PLFunction> f2;
      std::shared_ptr<PLFunction> f3;
      std::shared_ptr<PLFunction> f4;
      std::shared_ptr<PLFunction> fromF1;
      std::shared_ptr<PLFunction> fromF2;
      std::shared_ptr<PLFunction> fromC1;
      std::shared_ptr<PLFunction> fromC2;
      maxInventory = (params->endDayInventories) ? params->cli[client].maxInventory : std::max(0.0, params->cli[client].maxInventory - params->cli[client].dailyDemand[scenario][t - 1]);
      maxDeliverable = (params->endDayInventories) ? params->cli[client].maxInventory + params->cli[client].theoricalMinDemand : params->cli[client].maxInventory;
      
      //(1). First case: no delivery and no stockout
      //q(t) = 0 and  I(t-1) > daily[t] ==> I(t) > 0 
      //dp <--- C(t)(I(t-1) - daily(t)) + inventoryCost * (I(t-1) - daily(t))
      f1 = C[scenario][t]->clone(); //f1(x) = Ct(x)
      f1->addHoldingf(params->cli[client].inventoryCost); //f1(x) = Ct(x) + hi * x
      f1->shiftRight(params->cli[client].dailyDemand[scenario][t]); //f1(x) = Ct(x - d_i^t) + hi * (x - d_i^t)

      for (unsigned int i = 0; i < f1->nbPieces; i++) {
        f1->pieces[i]->fromC_pre = C[scenario][t]->pieces[i]->clone();
        f1->pieces[i]->stockout = false; 
        f1->pieces[i]->fromF = nullptr;
        f1->pieces[i]->fromC = nullptr;
        f1->pieces[i]->fromInst = nullptr;
      }
      f1 = f1->getInBound(params->cli[client].dailyDemand[scenario][t], maxInventory); //x>=d_i^t

      //(2) q(t) == 0,  I(t-1) < daily[t] ==> I(t) = 0 
      //dp <--- C(t)(0) + stockoutCost * (daily-I)
      f2 = std::make_shared<PLFunction>(params);
      shared_ptr<LinearPiece> tmpNoDelivery(make_shared<LinearPiece>(0, params->cli[client].dailyDemand[scenario][t] * params->cli[client].stockoutCost, params->cli[client].dailyDemand[scenario][t], 0));
      f2->append(tmpNoDelivery); //f2(x) = stockoutCost * (daily - x), x < d_i^t
      for (unsigned int i = 0; i < f2->nbPieces; i++){
        f2->pieces[i]->fromC_pre = C[scenario][t]->pieces[i]->clone();
        f2->pieces[i]->stockout = true; 
      }
      f2->moveUp(C[scenario][t]->pieces[0]->p1->y); //f2(x) = C_t(0) + stockoutCost * (daily - x), x < d_i^t
      f2 = f2->getInBound(0, maxInventory);
      f2 = f2->getInBound(0, params->cli[client].dailyDemand[scenario][t]); //x < d_i^t

      // q != 0, I(t-1)+q-daily > 0     ==> I(t)>0
      // C(t-1)(I(t-1)) <----- C(t)(I(t-1) + q - daily(t)) + holdingCost * (I(t-1) + q - daily(t)) + F_t^2(q)
      f3 = std::make_shared<PLFunction>(params);
      fromF2 = F[scenario][t-1]->clone();  //fromF2(q) = F2(q)
      fromF2 = fromF2->getInBound(-maxDeliverable, -1); //q < maxInventory
      
      fromC1 = C[scenario][t]->clone(); //fromC1(x) = C_t(x)
      fromC1->addHoldingf(params->cli[client].inventoryCost); //fromC1(x) = C_t(x) + hi*x
      fromC1->shiftRight(params->cli[client].dailyDemand[scenario][t]); //fromC1(x) = C_t(x-d_i^t) + hi*(x-d_i^t)

      for (unsigned int i = 0; i < fromC1->nbPieces; i++){
        fromC1->pieces[i]->fromC_pre = C[scenario][t]->pieces[i]->clone();
      }
      fromC1 = fromC1->getInBound(params->cli[client].dailyDemand[scenario][t], maxDeliverable);

      f3 = supperposition(fromC1, fromF2); //f3(x) = min_q (fromC1(x-q) + fromF2(q)) (q<0)
      for (unsigned int i = 0; i < f3->nbPieces; i++){
        f3->pieces[i]->fromC_pre = f3->pieces[i]->fromC->fromC_pre;
        f3->pieces[i]->stockout = false; 
      } 
      f3 = f3->getInBound(0, maxInventory);

      // q != 0, I(t-1)+q-daily < 0     ==> I(t)=0
      // C(t-1)(I(t-1)) <----- C(t)(0) + stockoutCost * (d_i^t-q-I(t-1)) + F_t^1(q)

      f4 = std::make_shared<PLFunction>(params);
      fromC2 = std::make_shared<PLFunction>(params);
      shared_ptr<LinearPiece> tmpDelivery(make_shared<LinearPiece>(0, params->cli[client].stockoutCost * params->cli[client].dailyDemand[scenario][t], params->cli[client].dailyDemand[scenario][t], 0));
      tmpDelivery->fromC_pre = C[scenario][t]->pieces[0]->clone();

      fromC2->append(tmpDelivery); //fromC2(x) = stockoutCost * (d_i^t-x)
      fromC2->moveUp(C[scenario][t]->pieces[0]->p1->y); //add C_t(0)
      fromC2 = fromC2->getInBound(0, maxDeliverable);
      
      fromF1 = F[scenario][t-1]->clone(); //fromF1(q) = F_1^t(q)
      fromF1 = fromF1->getInBound(- std::min(params->cli[client].dailyDemand[scenario][t], maxDeliverable), -1); //q < d_i^t
      
      f4 = supperposition(fromC2, fromF1); //f4(x) = min_q (fromC2(x-q) + fromF1(q)) (q<0)

      f4 = f4->getInBound(0, std::min(params->cli[client].dailyDemand[scenario][t], maxInventory));

      for (unsigned int i = 0; i < f4->nbPieces; i++){
        f4->pieces[i]->fromC_pre = f4->pieces[i]->fromC->fromC_pre; 
        f4->pieces[i]->stockout = true; 
      }
      C[scenario][t-1] = min_final(min_final(f1,f3), min_final(f2,f4))->clone(); 

    }
}

void LotSizingSolver::day1(unsigned int scenario) {
  unsigned int nbInsertionsFirstDay = (unsigned int) insertions[0][0].size();
  double maxInventory = params->cli[client].startingInventory;
  double maxDeliverable = (params->endDayInventories) ? params->cli[client].maxInventory + params->cli[client].theoricalMinDemand : params->cli[client].maxInventory;

  for (unsigned int k = 0; k < nbInsertionsFirstDay; k++) {
      //picewise linear functions
      std::shared_ptr<PLFunction> f3;
      std::shared_ptr<PLFunction> f4;
      std::shared_ptr<PLFunction>  fromF1;
      std::shared_ptr<PLFunction>  fromF2;
      std::shared_ptr<PLFunction>  fromC1;
      std::shared_ptr<PLFunction>  fromC2;

      // q != 0, I(t-1)+q-daily > 0     ==> I(t)>0
      // C(t-1)(I(t-1)) <----- C(t)(I(t-1) + q - daily(t)) + holdingCost * (I(t-1) + q - daily(t)) + F_t^2(q)
      f3 = std::make_shared<PLFunction>(params);
      fromF2 = Fk[scenario][k]->clone();  //fromF2(q) = F2(q)
      fromF2 = fromF2->getInBound(-maxDeliverable, -1); //q < maxInventory
      
      fromC1 = C[scenario][1]->clone(); //fromC1(x) = C_t(x)
      fromC1->addHoldingf(params->cli[client].inventoryCost); //fromC1(x) = C_t(x) + hi*x
      fromC1->shiftRight(params->cli[client].dailyDemand[scenario][1]); //fromC1(x) = C_t(x-d_i^t) + hi*(x-d_i^t)

      for (unsigned int i = 0; i < fromC1->nbPieces; i++){
        fromC1->pieces[i]->fromC_pre = C[scenario][1]->pieces[i]->clone();
      }
      fromC1 = fromC1->getInBound(params->cli[client].dailyDemand[scenario][1], maxDeliverable);

      f3 = supperposition(fromC1, fromF2); //f3(x) = min_q (fromC1(x-q) + fromF2(q)) (q<0)
      for (unsigned int i = 0; i < f3->nbPieces; i++){
        f3->pieces[i]->stockout = false; 
        f3->pieces[i]->fromC_pre = f3->pieces[i]->fromC->fromC_pre;
      } 
      f3 = f3->getInBound(0, maxInventory);

      // q != 0, I(t-1)+q-daily < 0     ==> I(t)=0
      // C(t-1)(I(t-1)) <----- C(t)(0) + stockoutCost * (d_i^t-q-I(t-1)) + F_t^1(q)

      f4 = std::make_shared<PLFunction>(params);
      fromC2 = std::make_shared<PLFunction>(params);
      shared_ptr<LinearPiece> tmpDelivery(make_shared<LinearPiece>(0, params->cli[client].stockoutCost * params->cli[client].dailyDemand[scenario][1], params->cli[client].dailyDemand[scenario][1], 0));
      tmpDelivery->fromC_pre = C[scenario][1]->pieces[0]->clone();

      fromC2->append(tmpDelivery); //fromC2(x) = stockoutCost * (d_i^t-x)
      fromC2->moveUp(C[scenario][1]->pieces[0]->p1->y); //add C_t(0)
      fromC2 = fromC2->getInBound(0, maxDeliverable);
      
      fromF1 = Fk[scenario][k]->clone(); //fromF1(q) = F_1^t(q)
      fromF1 = fromF1->getInBound(-std::min(params->cli[client].dailyDemand[scenario][1], maxDeliverable),-1); //q < d_i^t
      
      f4 = supperposition(fromC2, fromF1); //f4(x) = min_q (fromC2(x-q) + fromF1(q)) (q<0)

      f4 = f4->getInBound(0, std::min(params->cli[client].dailyDemand[scenario][1], maxInventory));

      for (unsigned int i = 0; i < f4->nbPieces; i++){
        f4->pieces[i]->stockout = true; 
        f4->pieces[i]->fromC_pre = f4->pieces[i]->fromC->fromC_pre; 
      }
      std::shared_ptr<PLFunction> temp = min_final(f3, f4)->clone(); 
      temp = temp->getInBound(params->cli[client].startingInventory, params->cli[client].startingInventory);
      double useless;
      std::shared_ptr<LinearPiece> tmp = temp->getMinimalPiece(useless, Ck[scenario][k]);
      savedCk[scenario][k] = temp;
    }
      //picewise linear functions
      std::shared_ptr<PLFunction> f1;
      std::shared_ptr<PLFunction> f2;
      
      //(1). First case: no delivery and no stockout
      //q(t) = 0 and  I(t-1) > daily[t] ==> I(t) > 0 
      //dp <--- C(t)(I(t-1) - daily(t)) + inventoryCost * (I(t-1) - daily(t))
      f1 = C[scenario][1]->clone(); //f1(x) = Ct(x)
      f1->addHoldingf(params->cli[client].inventoryCost); //f1(x) = Ct(x) + hi*x
      f1->shiftRight(params->cli[client].dailyDemand[scenario][1]); //f1(x) = Ct(x-d_i^t) + hi*(x-d_i^t)

      for (unsigned int i = 0; i < f1->nbPieces; i++) {
        f1->pieces[i]->fromC_pre = C[scenario][1]->pieces[i]->clone();
        f1->pieces[i]->stockout = false; 
        f1->pieces[i]->fromF = nullptr;
        f1->pieces[i]->fromC = nullptr;
        f1->pieces[i]->fromInst = nullptr;
      }
      f1 = f1->getInBound(params->cli[client].dailyDemand[scenario][1], maxInventory); //x>=d_i^t

      //(2) q(t) == 0,  I(t-1) < daily[t] ==> I(t) = 0 
      //dp <--- C(t)(0) + stockoutCost * (daily-I)
      f2 = std::make_shared<PLFunction>(params);
      shared_ptr<LinearPiece> tmpNoDelivery(make_shared<LinearPiece>(0, params->cli[client].dailyDemand[scenario][1] * params->cli[client].stockoutCost, params->cli[client].dailyDemand[scenario][1], 0));
      f2->append(tmpNoDelivery); //f2(x) = stockoutCost * (daily - x), x < d_i^t
      for (unsigned int i = 0; i < f2->nbPieces; i++){
        f2->pieces[i]->fromC_pre = C[scenario][1]->pieces[i]->clone();
        f2->pieces[i]->stockout = true; 
      }
      f2->moveUp(C[scenario][1]->pieces[0]->p1->y); //f2(x) = C_t(0) + stockoutCost * (daily - x), x < d_i^t
      f2 = f2->getInBound(0, maxInventory);
      f2 = f2->getInBound(0, params->cli[client].dailyDemand[scenario][1]); //x < d_i^t
      std::shared_ptr<PLFunction> temp = min_final(f1,f2)->clone(); 
      temp = temp->getInBound(params->cli[client].startingInventory, params->cli[client].startingInventory);

      double useless;
      std::shared_ptr<LinearPiece> tmp = temp->getMinimalPiece(useless, Ck[scenario][nbInsertionsFirstDay]);
      savedCk[scenario][nbInsertionsFirstDay] = temp;
}

LotSizingSolver::~LotSizingSolver() {}
