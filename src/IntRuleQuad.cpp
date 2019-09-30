/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "tpanic.h"

IntRuleQuad::IntRuleQuad(): IntRule(){
}

IntRuleQuad::IntRuleQuad(int order): IntRule() {
    SetOrder(order);
}

void IntRuleQuad::SetOrder(int order) {
    fOrder = order;

    IntRule1d dummy(order);
    int nPoints = dummy.NPoints();

    fPoints.Resize(nPoints*nPoints,2);
    fWeights.resize(nPoints*nPoints);

    VecDouble x_coord(1),y_coord(1);
    double x_weight,y_weight;
    for(int i=0;i<nPoints;i++){  dummy.Point(i,x_coord,x_weight); std::cout <<"x_coord  " <<x_coord[0] << "weight  " << x_weight << std::endl;
        for(int j=0;j<nPoints;j++){
            dummy.Point(j,y_coord,y_weight);
            fPoints(nPoints*i+j,0) = x_coord[0];
            fPoints(nPoints*i+j,1) = y_coord[0];
            fWeights[nPoints*i+j] = x_weight*y_weight;
        }
    }
}

void IntRuleQuad::gaulegQuad(const double x1, const double x2, VecDouble &co, VecDouble &w) {
    //fPoints and fWeights are specified according to the order using the 1d Quadrature (see SetOrder above, IntRule1d::SetOrder and IntRule1d::gauleg). Therefore, as the code was constructed, a method gaulegQuad is not necessary and even meaningless.
    IntRule1d dummy;
    dummy.gauleg(x1,x2,co,w);
}
