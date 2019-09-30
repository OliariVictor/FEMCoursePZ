/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "IntRule.h"
#include "tpanic.h"

using namespace std;

IntRule::IntRule() : fOrder(0), fPoints( 0,0), fWeights(0){

}

IntRule::IntRule(int order) : fOrder(order), fPoints(0,0), fWeights(0) {
}

IntRule::~IntRule() {

}

IntRule::IntRule(const IntRule& copy) : fOrder(copy.fOrder),fPoints(copy.fPoints),fWeights(copy.fWeights) {
}

IntRule &IntRule::operator=(const IntRule &cp) {
    fOrder = cp.fOrder;
    fPoints = cp.fPoints;
    fWeights = cp.fWeights;
    return *this;
}

int IntRule::NPoints() const {
    int nPoints;

    //If fOrder is even...Else...
    if (fOrder%2 == 0){
        nPoints = (int)(fOrder/2)+1;
    }
    else{
        nPoints = (int)((fOrder-1)/2+1);
    }

    return nPoints;
}

void IntRule::Print(std::ostream &out) {
    out << "\nGaussian Legendre Quadrature:\n";
    out << "Polynomial Order = " << fOrder << "\n";
    int colnNum = fPoints.Cols();
    int nPoints = fPoints.Rows();
    double sum = 0;
    switch(colnNum) {
        case 1:
            int i;
            for (i = 0; i < nPoints; i++) {
                out << "x_coord = " << fPoints(i, 0) << "  ||  weight =" << fWeights[i] << "\n";
                sum +=fWeights[i];
            } out << "Sum of all weights = " << sum;
            break;
        case 2:
            for (int i = 0; i < nPoints; i++) {
                out << "x_coord = " << fPoints(i, 0) << "  &&  y_coord = " << fPoints(i, 1) << "  ||  weight ="
                    << fWeights[i] << "\n";
                sum +=fWeights[i];
            } out << "Sum of all weights = " << sum;
            break;
        case 3:
            for (int i = 0; i < nPoints; i++) {
                out << "x_coord = " << fPoints(i, 0) << "  &&  y_coord = " << fPoints(i, 1) << "  &&  z_coord = "
                    << fPoints(i, 2) << "\n";
                sum +=fWeights[i];
            } out << "Sum of all weights = " << sum;
            break;
    }
}


void IntRule::Point(int p, VecDouble& co, double& w) const {
    int coSize = fPoints.Cols();
    for (int i = 0; i < coSize; i++) co[i] = fPoints.GetVal(p, i);
    w = fWeights[p];
}

