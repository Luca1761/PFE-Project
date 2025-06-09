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
LotSizingSolver::LotSizingSolver(Params *params, vector<vector<Insertion>> inst, int client)
    : params(params), insertions(inst), client(client)
{
  horizon = (int)insertions.size();
  
  F1 = vector<std::shared_ptr<PLFunction>>(horizon);
  F2 = vector<std::shared_ptr<PLFunction>>(horizon); 

  for (int t = 0; t < horizon; t++)
  {
    vector<Insertion> tmp = insertions[t]; // all possible place in day t
  
    std::shared_ptr<PLFunction> a = std::make_shared<PLFunction>(params, tmp, t, client, params->cli[client].dailyDemand[t]); 

    std::shared_ptr<PLFunction> b = std::make_shared<PLFunction>(params, tmp, t, client); 
    F1[t] = std::make_shared<PLFunction>(params);
    F2[t] = std::make_shared<PLFunction>(params);
    
    for (int i = a->nbPieces-1; i >=0 ; i--) {
      std::shared_ptr<LinearPiece> tmp = a->pieces[i]->clone();
      tmp->updateLinearPiece(-tmp->p2->x, tmp->p2->y, -tmp->p1->x, tmp->p1->y);
      F1[t]->append(tmp);
    }
    
    for (int i = b->nbPieces-1; i >=0 ; i--) {
      std::shared_ptr<LinearPiece> tmp = b->pieces[i]->clone();
      tmp->updateLinearPiece(-tmp->p2->x, tmp->p2->y, -tmp->p1->x, tmp->p1->y);
      F2[t]->append(tmp);
    }

    for (int i = 0; i < F1[t]->nbPieces; i++){
        F1[t]->pieces[i]->fromF = F1[t]->pieces[i]->clone();
        F1[t]->pieces[i]->replenishment_loss = 0;
        F1[t]->pieces[i]->fromC = nullptr;
        F1[t]->pieces[i]->fromC_pre = nullptr;
        F1[t]->pieces[i]->fromInst = make_shared<Insertion>(
            F1[t]->pieces[i]->fromInst->detour, 
            F1[t]->pieces[i]->fromInst->load,
            F1[t]->pieces[i]->fromInst->place);
    }
    for (int i = 0; i < F2[t]->nbPieces; i++){
        F2[t]->pieces[i]->fromF = F2[t]->pieces[i]->clone(); 
        F2[t]->pieces[i]->replenishment_loss = -1;
        F2[t]->pieces[i]->fromC = nullptr;
        F2[t]->pieces[i]->fromC_pre = nullptr;
        F2[t]->pieces[i]->fromInst = make_shared<Insertion>(
            F2[t]->pieces[i]->fromInst->detour, 
            F2[t]->pieces[i]->fromInst->load,
            F2[t]->pieces[i]->fromInst->place);
    }
      
  }
  
}

std::shared_ptr<PLFunction> LotSizingSolver::copyPLFunction(
    std::shared_ptr<PLFunction> source)
{
  std::shared_ptr<PLFunction> destination(std::make_shared<PLFunction>(params));
  destination->nbPieces = 0;
  destination->pieces = vector<shared_ptr<LinearPiece>>();
  for (int i = 0; i < source->nbPieces; i++)
  {
    destination->append(source->pieces[i]);
  }

  if (source->pieceAt0 != nullptr)
  {
    destination->pieceAt0 = source->pieceAt0->clone();
  }
  destination->valueAt0 = source->valueAt0;

  return destination;
}

void LotSizingSolver::extractBreakpoints(const std::shared_ptr<PLFunction>& f, vector<double>& breakpoints) {
    for (int i = 0; i < f->nbPieces; i++) {
        if( neq(f->pieces[i]->p1->x,f->pieces[i]->p2->x) &&
              ( i== 0 ||  neq(f->pieces[i]->p1->x,f->pieces[i-1]->p2->x) ) )
            breakpoints.push_back(f->pieces[i]->p1->x);
        breakpoints.push_back(f->pieces[i]->p2->x);
    }
}
void LotSizingSolver::extractRepeat(const std::shared_ptr<PLFunction>& f, vector<double>& breakpoints) {
    for (int i = 0; i < f->nbPieces; i++) {
        if(fabs(f->pieces[i]->p1->x-f->pieces[i]->p2->x)<0.00001)
          breakpoints.push_back(f->pieces[i]->p1->x);
    }
}

vector<double> LotSizingSolver::getBreakpoints_final(std::shared_ptr<PLFunction> f1, std::shared_ptr<PLFunction> f2){
    vector<double> breakpoints;
    vector<double> repeat;
    extractBreakpoints(f1, breakpoints);
    extractBreakpoints(f2, breakpoints);
    extractRepeat(f1, repeat);
    extractRepeat(f2, repeat);
    std::sort(breakpoints.begin(), breakpoints.end());
    std::sort(repeat.begin(), repeat.end());
  auto it = std::unique(breakpoints.begin(), breakpoints.end(), 
                      [](double a, double b) { return fabs(a - b) < 0.00001; });
  breakpoints.erase(it, breakpoints.end());

  auto it1 = std::unique(repeat.begin(), repeat.end(), 
                        [](double a, double b) { return fabs(a - b) < 0.00001; });
  repeat.erase(it1, repeat.end());

    breakpoints.insert(breakpoints.end(), repeat.begin(), repeat.end());
    
    std::sort(breakpoints.begin(), breakpoints.end());
    return breakpoints;
}

std::shared_ptr<LinearPiece> LotSizingSolver::createPieceFromLowerY(
    std::shared_ptr<LinearPiece> &chosenPiece,
    double x1, double y1, double x2, double y2) 
{
    std::shared_ptr<LinearPiece> tmpPiece = std::make_shared<LinearPiece>(x1, y1, x2, y2);
    tmpPiece->fromC_pre = chosenPiece->fromC_pre;
    tmpPiece->fromC = chosenPiece->fromC;
    tmpPiece->fromF = chosenPiece->fromF;
    tmpPiece->fromInst = chosenPiece->fromInst;
    tmpPiece->replenishment_loss = chosenPiece->replenishment_loss;

    return tmpPiece;
}


std::shared_ptr<PLFunction> LotSizingSolver::min_final(std::shared_ptr<PLFunction> f1, std::shared_ptr<PLFunction> f2) {
  std::shared_ptr<PLFunction> f = std::make_shared<PLFunction>(params);
    if (f1->nbPieces == 0 && f2->nbPieces == 0) return f;
    if (f1->nbPieces == 0 && f2->nbPieces > 0) return copyPLFunction(f2);
    if (f1->nbPieces > 0 && f2->nbPieces == 0) return copyPLFunction(f1);

    vector<double> breakpoints = getBreakpoints_final(f1, f2);
    int nbBreakPoints = breakpoints.size();
    std::shared_ptr<LinearPiece> lp1 = f1->pieces[0];
    std::shared_ptr<LinearPiece> lp2 = f2->pieces[0];
    int n1=0,n2=0;
    for (int i = 1; i < nbBreakPoints; i++) {
      if(lp1){lp1->updateLinearPiece(lp1->p1->x,lp1->p1->y,lp1->p2->x,lp1->p2->y);}
      if(lp2)lp2->updateLinearPiece(lp2->p1->x,lp2->p1->y,lp2->p2->x,lp2->p2->y);
      std::shared_ptr<LinearPiece> piece1 = nullptr;
      
      if (lp1) piece1 = lp1->getInpiece(breakpoints[i - 1], breakpoints[i]); 
      
      std::shared_ptr<LinearPiece> piece2 = nullptr;
      if (lp2)  piece2 = lp2->getInpiece(breakpoints[i - 1], breakpoints[i]); 
      
      if (!piece1 && !piece2);
      else if (!piece1)  f->append(piece2);
      else if (!piece2) f->append(piece1);
      else {
        double x, y;
        bool intersects = f->intersect(piece1, piece2, x, y);
            if (eq(piece2->p1->x, piece2->p2->x) && eq(piece1->p2->x, piece2->p1->x)) {
              f->append(piece1);
              if (gt(piece1->p2->y, piece2->p1->y)) f->append(piece2);
            } else if (eq(piece1->p1->x, piece1->p2->x) && eq(piece2->p2->x, piece1->p1->x)) {
              f->append(piece2);
              if (gt(piece2->p2->y, piece1->p1->y)) f->append(piece1);
            } 
            else if (intersects) {
                std::shared_ptr<LinearPiece> firstPiece, secondPiece, chosenPiece;
                chosenPiece = ( piece1->p1->y <= piece2->p1->y) ? piece1 : piece2;
                if (round(x) == chosenPiece->p2->x) {
                  firstPiece = createPieceFromLowerY(chosenPiece, chosenPiece->p1->x, chosenPiece->p1->y, chosenPiece->p2->x, chosenPiece->p2->y);
                  f->append(firstPiece);

                  chosenPiece =( (piece1->p2->y < piece2->p2->y) )? piece1 : piece2;
                  secondPiece = createPieceFromLowerY(chosenPiece, chosenPiece->p2->x, chosenPiece->p2->y, chosenPiece->p2->x, chosenPiece->p2->y);
                  f->append(secondPiece);
                } else {
                  firstPiece = createPieceFromLowerY(chosenPiece, chosenPiece->p1->x, chosenPiece->p1->y, x, y);
                  f->append(firstPiece);

                  chosenPiece =( (piece1->p2->y < piece2->p2->y) )? piece1 : piece2;
                  secondPiece = createPieceFromLowerY(chosenPiece, x, y, chosenPiece->p2->x, chosenPiece->p2->y);
                  f->append(secondPiece);
                }
            } 
            else if ( eq(piece1->p1->x ,x) ){
                if(piece1->p2->y < piece2->p2->y) 
                    f->append(piece1);
                 else 
                    f->append(piece2);
            }
            else if ( eq(piece1->p2->x ,x) ){
                if(piece1->p1->y <piece2->p1->y) 
                    f->append(piece1);
                 else 
                    f->append(piece2);
            }
            else {
              if ( (piece1->p1->y -piece2->p1->y) + (piece1->p2->y -piece2->p2->y) < 0) 
                    f->append(piece1);
                else 
                    f->append(piece2);
                
            }
           
        }
        if (lp1&& fabs(breakpoints[i]- lp1->p2->x)<0.00001 ){
            if(n1 + 1 < f1->nbPieces) lp1 = f1->pieces[++n1];
            else lp1 = nullptr;
        }
        if (lp2&& fabs(breakpoints[i]- lp2->p2->x)<0.00001 ){
          if(n2 + 1 < f2->nbPieces) lp2 = f2->pieces[++n2];
          else lp2 = nullptr;
        }
        if (!lp1 && !lp2) break;
    }
    return f;
}





std::shared_ptr<PLFunction> LotSizingSolver::supperpositionPieces(
    std::shared_ptr<LinearPiece> fromPieceC,
    std::shared_ptr<LinearPiece> fromPieceF)
{
  std::shared_ptr<PLFunction> f(std::make_shared<PLFunction>(params));
  vector<Point> points;
  shared_ptr<LinearPiece> tmpPiece,newpiece;
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
    if(gt(minYPoint1.x,minYPoint2.x) )std::swap(minYPoint1,minYPoint2);

    tmpPiece = make_shared<LinearPiece>(minYPoint1.x, minYPoint1.y,minYPoint2.x, minYPoint2.y);
    f->append(tmpPiece); 
    f->pieces[0]->fromC = fromPieceC->clone();
    f->pieces[0]->fromF = fromPieceF->clone();
    f->pieces[0]->fromInst = fromPieceF->fromInst; 
    return f;
  }

  int minYIndex1 = std::find(points.begin(), points.end(), minYPoint1) - points.begin();
  int minYIndex2 = std::find(points.begin(), points.end(), minYPoint2) - points.begin();
  int next1 = (minYIndex1 + 1) % 4,next2 = (minYIndex1+3)%4;

  double slope1 = (eq(points[next1].x, points[minYIndex1].x )||   eq(points[next1].y, points[minYIndex1].y ) )
                    ? 0 : (points[next1].y- points[minYIndex1].y )/(points[next1].x- points[minYIndex1].x );
  double slope2 = ( eq(points[next2].x, points[minYIndex1].x )||  eq(points[next2].y, points[minYIndex1].y ))
                    ? 0 : (points[next2].y- points[minYIndex1].y )/(points[next2].x- points[minYIndex1].x );
  if(eq(slope1, 0) || eq(slope2 ,0)){//one 0
    if(minYPoint1.x>minYPoint2.x )std::swap(minYIndex1,minYIndex2);
    next1 = (minYIndex1 + 1) % 4; next2 = (minYIndex1+3)%4;
    int rightup = 6-minYIndex1-next1-next2;
    int leftup = 6-minYIndex1-minYIndex2-rightup;

    if(lt(points[rightup].x,points[minYIndex2].x)){
        tmpPiece = make_shared<LinearPiece>(points[leftup].x, points[leftup].y,points[minYIndex1].x, points[minYIndex1].y); 
        newpiece = make_shared<LinearPiece>(points[minYIndex1].x, points[minYIndex1].y,  points[minYIndex2].x, points[minYIndex2].y);
    }
    else{
        tmpPiece = make_shared<LinearPiece>(points[minYIndex1].x, points[minYIndex1].y,  points[minYIndex2].x, points[minYIndex2].y); 
        newpiece = make_shared<LinearPiece>(points[minYIndex2].x, points[minYIndex2].y,points[rightup].x, points[rightup].y);
    }
  }  
  else if (lt(slope1,0) && lt(0,slope2) || lt(slope2,0) && lt(0,slope1)){//one negative one postive
    if(points[next1].x>points[next2].x)std::swap(next1,next2);
    tmpPiece = make_shared<LinearPiece>(points[next1].x, points[next1].y,points[minYIndex1].x, points[minYIndex1].y); 
    newpiece = make_shared<LinearPiece>(points[minYIndex1].x, points[minYIndex1].y,  points[next2].x, points[next2].y);
    
  }
  else if(lt(slope1,0)){//both negative
    if(eq(minYPoint2.x,minYPoint3.x) && eq(minYPoint2.y,minYPoint3.y) ){
      tmpPiece = make_shared<LinearPiece>(minYPoint4.x, minYPoint4.y,minYPoint2.x, minYPoint2.y); 
      newpiece = make_shared<LinearPiece>(minYPoint2.x, minYPoint2.y,minYPoint1.x, minYPoint1.y);  
    }
    else{
      int tmp1 = gt(slope1,slope2)?next1:next2;
      int tmp2 = 6-minYIndex1-next1-next2;
      tmpPiece = make_shared<LinearPiece>(points[tmp2].x, points[tmp2].y,points[tmp1].x, points[tmp1].y); 
      newpiece = make_shared<LinearPiece>(points[tmp1].x, points[tmp1].y,points[minYIndex1].x, points[minYIndex1].y);  
    }
  }
  else{//both positive
    if(eq(minYPoint2.x,minYPoint3.x) && eq(minYPoint2.y,minYPoint3.y) ){//line
      tmpPiece = make_shared<LinearPiece>(minYPoint1.x, minYPoint1.y,minYPoint2.x, minYPoint2.y); 
      newpiece = make_shared<LinearPiece>(minYPoint2.x, minYPoint2.y,minYPoint4.x, minYPoint4.y);  
    }
    else{
      int tmp1 = gt(slope1,slope2)?next2:next1;
      int tmp2 = 6-minYIndex1-next1-next2;
     
      tmpPiece = make_shared<LinearPiece>(points[minYIndex1].x, points[minYIndex1].y,points[tmp1].x, points[tmp1].y); 
      newpiece = make_shared<LinearPiece>(points[tmp1].x, points[tmp1].y,points[tmp2].x, points[tmp2].y);  
    }
  }
  f->append(tmpPiece);
  f->append(newpiece);
  for(int i = 0 ; i < f->nbPieces ; i++){
    f->pieces[i]->fromC = fromPieceC->clone();
    f->pieces[i]->fromF = fromPieceF->clone();
    f->pieces[i]->fromInst = fromPieceF->fromInst; 
  }
  return f;
}

std::shared_ptr<PLFunction> LotSizingSolver:: supperposition(
    std::shared_ptr<PLFunction> fromC, std::shared_ptr<PLFunction> fromF)
{
 std::shared_ptr<PLFunction> f(std::make_shared<PLFunction>(params));

  for (int i = 0; i < fromC->nbPieces; i++){
    for (int j = 0; j < fromF->nbPieces; j++){
      
      std::shared_ptr<PLFunction> tmpF = std::make_shared<PLFunction>(params);
      std::shared_ptr<LinearPiece> fromPieceC = fromC->pieces[i]->clone();
      
      std::shared_ptr<LinearPiece> fromPieceF = fromF->pieces[j]->clone();

      tmpF = supperpositionPieces(fromPieceC, fromPieceF);
      
      if (f->nbPieces > 0 || tmpF->nbPieces > 0){
        std::shared_ptr<PLFunction> fmin;
        fmin = min_final(f, tmpF);
        f = copyPLFunction(fmin);
      }

      tmpF.reset();
      fromPieceC.reset();
      fromPieceF.reset();
    }
   
  }
  
  return f;
}

void LotSizingSolver::solveEquationSystemHoldingBackward(std::shared_ptr<LinearPiece> C,
                                          std::shared_ptr<LinearPiece> fromC,
                                          std::shared_ptr<LinearPiece> fromF,
                                          double I, double demand,
                                          double &fromI, double &quantity,std::shared_ptr<LinearPiece> Cpre,double hodling){
  if(eq(fromC->p1->x,fromC->p2->x )){
    fromI = round(fromC->p2->x - demand);
    quantity = round(fromI + demand - I);
    return;
  }
  if (eq(fromF->p1->x, fromF->p2->x)) {
    quantity = round(-fromF->p1->x);
    fromI = round(I - demand + quantity);
    return;
  }

  double slopeC = (fromC->p2->y - fromC->p1->y) / (fromC->p2->x - fromC->p1->x);
  double slopeF = -(fromF->p2->y - fromF->p1->y) / (fromF->p2->x - fromF->p1->x);
  
  if(eq(slopeC,-slopeF)){
    double upperbound = std::min<double>(fromC->p2->x-I,-fromF->p1->x);
    double lowerbound = std::max<double>(fromC->p1->x-I,-fromF->p2->x);
    
    if (ge(upperbound ,lowerbound) ){
        quantity = std::floor(upperbound);
        fromI = round(I - demand + quantity);
    }
    else{
       quantity = std::round(upperbound);
        fromI = round(I - demand + quantity);
    }
    return;
  } 
  slopeC*=10000;slopeF*=10000;
  double x1 = fromC->p2->x, y1 = fromC->p2->y,x2 = -fromF->p1->x,y2 = fromF->p1->y;
  double numerator = C->cost(std::max<double>(0,I)) *10000- y1*10000 - y2*10000;
    numerator -= slopeC * (I-x1);
    numerator += slopeF * x2;
    quantity = round(numerator / (slopeF + slopeC));
    double left = round(std::max<double>(fromC->p1->x-I,-fromF->p2->x));
    double right = round(std::min<double>(fromC->p2->x-I,-fromF->p1->x));
    if(gt( quantity ,right )||lt(quantity,left)  ){
      if (gt(quantity ,right) )  quantity = right;
      if (lt(quantity,left) )  quantity = left;
    }
    
  fromI = round(I - demand + quantity);
}

bool LotSizingSolver::backtrackingStockoutBackward() {
  // initialization
  for (int i = 0; i < horizon; i++){
    quantities[i] = 0.0;
    breakpoints[i] = nullptr;
    I[i] = 0;
  }
  I.push_back(0.0);
  int day = 0;
  
  if (C[day]->nbPieces == 0)
    return false;
  std::shared_ptr<LinearPiece> tmp = 
      C[day]->getMinimalPiece_stockout(client, I[day], objective);

  ///////////////////////////////////////////
  while (tmp != nullptr){

    // if do not delivery any thing, then inventory at the end of next day
    // equals this day demand
    if (eq(tmp->replenishment_loss,-1) && !tmp->fromF ) //f1 q = 0
    {
      I[day + 1] = I[day] - params->cli[client].dailyDemand[day];
    }
    
    else if (neq(tmp->replenishment_loss,-1) && !tmp->fromF) //f2  q= 0, not enough
    { 
      if(I[day] - params->cli[client].dailyDemand[day] >= 0){
        // cout <<"lotsizing :: backward error on f2!!!!Iday - d >=0: "<<I[day] <<endl;
      }    
      I[day + 1] = 0.0;
    }
    
    else if (eq(tmp->replenishment_loss,-1) && tmp->fromF ) //f3
    {
      std::shared_ptr<LinearPiece> fromC = tmp->fromC->clone();
      std::shared_ptr<LinearPiece> fromF = tmp->fromF->clone();
      std::shared_ptr<LinearPiece> Cpre = tmp->fromC_pre->clone();
      solveEquationSystemHoldingBackward(tmp, fromC, fromF, I[day],
                            params->cli[client].dailyDemand[day], I[day + 1],
                            quantities[day],Cpre,params->cli[client].inventoryCost);
      
      shared_ptr<LinearPiece> tmpF = tmp->fromF;
      breakpoints[day] = tmpF->fromInst;
      if (std::isnan(quantities[day])) {
        std::cout << "from no stockout" << std::endl;
        fromC->print();
        fromF->print();
      }
    }
    else if (neq(tmp->replenishment_loss,-1) && tmp->fromF )//f4
    {
      std::shared_ptr<LinearPiece> fromC = tmp->fromC->clone();
      std::shared_ptr<LinearPiece> fromF = tmp->fromF->clone();
      std::shared_ptr<LinearPiece> Cpre = tmp->fromC_pre->clone();
      solveEquationSystemHoldingBackward(tmp, fromC, fromF, I[day],
                      params->cli[client].dailyDemand[day], I[day + 1],
                      quantities[day],Cpre,params->cli[client].inventoryCost);
      
      I[day + 1] = 0.0;

      shared_ptr<LinearPiece> tmpF = tmp->fromF;
      breakpoints[day] = tmpF->fromInst;
      if (std::isnan(quantities[day])) {
        std::cout << "from stockout" << std::endl;
        fromC->print();
        fromF->print();
      }
    }
    tmp = tmp->fromC_pre;

    day = day + 1;
    if(day == horizon) break;
  }
  return true;
}

void  LotSizingSolver::Lastday(vector<std::shared_ptr<PLFunction>> &C){//for the last day
  double maxInventory = params->cli[client].maxInventory - params->cli[client].dailyDemand[horizon - 1];
  C[horizon]->nbPieces = 0;
  C[horizon]->pieces = vector<shared_ptr<LinearPiece>>();
  double finalI = params->cli[client].startingInventory;
  double eps = 0.000000001;
  double realCost = 1000000000;
  std::shared_ptr<LinearPiece> tmp = std::make_shared<LinearPiece>(0, 0, maxInventory, 0);
  tmp->fromC = nullptr;
  tmp->fromC_pre = nullptr;
  tmp->replenishment_loss = -finalI;
  tmp->fromF = nullptr;
  tmp->fromInst = nullptr;
  C[horizon]->append(tmp); 
  return;

  std::shared_ptr<LinearPiece> tmptmp = std::make_shared<LinearPiece>(finalI-eps, 0, finalI + eps, 0);
  tmptmp->fromC = nullptr;
  tmptmp->fromC_pre = nullptr;
  tmptmp->replenishment_loss = -finalI;
  tmptmp->fromF = nullptr;
  tmptmp->fromInst = nullptr;
  C[horizon]->append(tmptmp); 

  std::shared_ptr<LinearPiece> tmptmptmp = std::make_shared<LinearPiece>(finalI + eps, 0, maxInventory, 0);
  tmptmptmp->fromC = nullptr;
  tmptmptmp->fromC_pre = nullptr;
  tmptmptmp->replenishment_loss = -finalI;
  tmptmptmp->fromF = nullptr;
  tmptmptmp->fromInst = nullptr;
  C[horizon]->append(tmptmptmp);
}

bool LotSizingSolver::solveStockoutBackward()
{
  double minInventory, maxInventory;

  C = vector<std::shared_ptr<PLFunction>>(horizon + 1);

  // final inventory at the end of each period
  I = vector<double>(horizon);

  // replinishment at each period
  quantities = vector<double>(horizon);

  breakpoints = vector<std::shared_ptr<Insertion>>(horizon);
  
  for (int i = 0; i < horizon + 1; i++) {
    C[i] = std::make_shared<PLFunction>(params); // a piecewise linear function in each period
  }

  Lastday(C);

  for (int t = horizon; t >= 1; t--)
  {
    if (neq(C[t]->pieces[0]->p1->x,0)) {
      throw string("Ct(0) doesn't exist!!!");
    }
    //picewise linear functions
    std::shared_ptr<PLFunction> f1;
    std::shared_ptr<PLFunction> f2;
    std::shared_ptr<PLFunction> f3;
    std::shared_ptr<PLFunction> f4;
    std::shared_ptr<PLFunction>  fromF;
    std::shared_ptr<PLFunction>  fromF1;
    std::shared_ptr<PLFunction>  fromF2;
    std::shared_ptr<PLFunction>  fromC1;
    std::shared_ptr<PLFunction>  fromC2;
    minInventory = params->cli[client].minInventory;
    maxInventory = (t==1) ? params->cli[client].startingInventory : params->cli[client].maxInventory - params->cli[client].dailyDemand[t - 2];
    
    //(1). First case: no delivery and no stockout
    //q(t) = 0 and  I(t-1) > daily[t] ==> I(t) > 0 
    //dp <--- C(t)(I(t-1) - daily(t)) + inventoryCost * (I(t-1) - daily(t))

    f1 = copyPLFunction(C[t]); //f1(x) = Ct(x)
    f1->addHoldingf(params->cli[client].inventoryCost); //f1(x) = Ct(x) + hi*x
    f1->shiftRight(params->cli[client].dailyDemand[t-1]); //f1(x) = Ct(x-d_i^t) + hi*(x-d_i^t)

    for (int i = 0; i < f1->nbPieces; i++)
    {
      f1->pieces[i]->fromC_pre = C[t]->pieces[i]->clone();
      f1->pieces[i]->replenishment_loss = -1; 
      f1->pieces[i]->fromF = nullptr;
      f1->pieces[i]->fromC = nullptr;
      f1->pieces[i]->fromInst = nullptr;
    }
    f1 = f1->getInBound(params->cli[client].dailyDemand[t-1], maxInventory); //x>=d_i^t


    //(2) q(t) == 0,  I(t-1) < daily[t] ==> I(t) = 0 
    //dp <--- C(t)(0) + stockoutCost * (daily-I)
    f2 = std::make_shared<PLFunction>(params);
    shared_ptr<LinearPiece> tmpNoDelivery(make_shared<LinearPiece>(0, params->cli[client].dailyDemand[t-1] * params->cli[client].stockoutCost, params->cli[client].dailyDemand[t-1], 0));
    f2->append(tmpNoDelivery); //f2(x) = stockoutCost * (daily - x), x < d_i^t
    for (int i = 0; i < f2->nbPieces; i++){
      f2->pieces[i]->fromC_pre = C[t]->pieces[i]->clone();
      f2->pieces[i]->replenishment_loss = 0; 
      f2->pieces[i]->fromC = nullptr;
      f2->pieces[i]->fromF = nullptr;
      f2->pieces[i]->fromInst = nullptr;
    }
    f2->moveUp(C[t]->pieces[0]->p1->y); //f2(x) = C_t(0) + stockoutCost * (daily - x), x < d_i^t
    
    f2 = f2->getInBound(0, params->cli[client].dailyDemand[t-1]); //x < d_i^t

    // q != 0, I(t-1)+q-daily > 0     ==> I(t)>0
    // C(t-1)(I(t-1)) <----- C(t)(I(t-1) + q - daily(t)) + holdingCost * (I(t-1) + q - daily(t)) + F_t^2(q)
    f3 = std::make_shared<PLFunction>(params);
    fromF2 = copyPLFunction(F2[t-1]);  //fromF2(q) = F2(q)
    fromF2 = fromF2->getInBound(-params->cli[client].maxInventory,-1); //q < maxInventory
    
    fromC1 = copyPLFunction(C[t]); //fromC1(x) = C_t(x)
    fromC1->addHoldingf(params->cli[client].inventoryCost); //fromC1(x) = C_t(x) + hi*x
    fromC1->shiftRight(params->cli[client].dailyDemand[t-1]); //fromC1(x) = C_t(x-d_i^t) + hi*(x-d_i^t)

    for (int i = 0; i < fromC1->nbPieces; i++){
      fromC1->pieces[i]->fromC_pre = C[t]->pieces[i]->clone();
    }
    fromC1 = fromC1->getInBound(params->cli[client].dailyDemand[t-1], params->cli[client].maxInventory);

    f3 = supperposition(fromC1, fromF2); //f3(x) = min_q (fromC1(x-q) + fromF2(q)) (q<0)
    for (int i = 0; i < f3->nbPieces; i++){
      f3->pieces[i]->replenishment_loss = -1; 
      f3->pieces[i]->fromC_pre = f3->pieces[i]->fromC->fromC_pre;
    } 
    f3 = f3->getInBound(0,maxInventory);

    // q != 0, I(t-1)+q-daily < 0     ==> I(t)=0
    // C(t-1)(I(t-1)) <----- C(t)(0) + stockoutCost * (d_i^t-q-I(t-1)) + F_t^1(q)

    f4 = std::make_shared<PLFunction>(params);
    fromC2 = std::make_shared<PLFunction>(params);
    shared_ptr<LinearPiece> tmpDelivery(make_shared<LinearPiece>(0, params->cli[client].stockoutCost * params->cli[client].dailyDemand[t-1], params->cli[client].dailyDemand[t-1], 0));
    fromC2->append(tmpDelivery); //fromC2(x) = stockoutCost * (d_i^t-x)
    fromC2 = fromC2->getInBound(0,params->cli[client].dailyDemand[t-1]); // x < d_i^t
    fromC2->moveUp(C[t]->pieces[0]->p1->y); //add C_t(0)
    
    fromF1 = copyPLFunction(F1[t-1]); //fromF1(q) = F_1^t(q)
    fromF1 = fromF1->getInBound(-(params->cli[client].dailyDemand[t-1]),-1); //q < d_i^t
    

    for (int i = 0; i < fromC2->nbPieces; i++){
      fromC2->pieces[i]->fromC_pre = C[t]->pieces[i]->clone();
    }
    
    f4 = supperposition(fromC2, fromF1); //f4(x) = min_q (fromC2(x-q) + fromF1(q)) (q<0)

    f4 = f4->getInBound(0,params->cli[client].dailyDemand[t-1]);

    for (int i = 0; i < f4->nbPieces; i++){
      f4->pieces[i]->replenishment_loss = 0; 
      f4->pieces[i]->fromC_pre = f4->pieces[i]->fromC->fromC_pre; 
    }

    C[t-1] = min_final(min_final(f1,f3),min_final(f2,f4)); 

  }
  C[0] = C[0]->getInBound(params->cli[client].startingInventory, params->cli[client].startingInventory);

  bool ok = backtrackingStockoutBackward();

  if (!ok)
  {
    return false;
  }
  return true;

}

LotSizingSolver::~LotSizingSolver()
{
  
}
