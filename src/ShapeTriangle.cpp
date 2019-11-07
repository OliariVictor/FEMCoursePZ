/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "ShapeTriangle.h"
#include "tpanic.h"

void ShapeTriangle::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi) {
    int no = 3, counter = 0, size = NShapeFunctions(orders);
    phi.resize(size); dphi.Resize(size,2);

    phi[0] = 1 - xi[0] - xi[1]; phi[1] = xi[0]; phi[2] = xi[1];
    dphi(0,0) = dphi(0,1) = -1;
    dphi(1,0) = dphi(2,1) = 1;
    dphi(1,1) = dphi(2,0) = 0;

    auto edge = [&](int i, int j, int z){
        phi[z] = phi[i]*phi[j];
        dphi(z,0) = dphi(i,0)*phi[j]+phi[i]*dphi(j,0);
        dphi(z,1) = dphi(i,1)*phi[j]+phi[i]*dphi(j,1);
    };

    while (no < nSides) {
        if (orders[no] != 1 && orders[no] != 2) { std::cout << "Requested order not defined in this scope";DebugStop();}
        if(no < 6){
            if (orders[no] == 1) counter++;
            else if(orders[no] == 2) {
                switch(no){
                    case 3: edge(0,1,3-counter);break;
                    case 4: edge(1,2,4-counter);break;
                    case 5: edge(2,0,5-counter);break;
                }
            } else {std::cout << "Please choose a valid shape order"; DebugStop;}
        }else if (no == 6){
            if (orders[no] == 2){
                phi[6-counter] = phi[0]*phi[1]*phi[2];
                dphi(6-counter,0) = dphi(0,0)*phi[1]*phi[2]+phi[0]*dphi(1,0)*phi[2]+phi[0]*phi[1]*dphi(2,0);
                dphi(6-counter,1) = dphi(0,1)*phi[1]*phi[2]+phi[0]*dphi(1,1)*phi[2]+phi[0]*phi[1]*dphi(2,1);
            }
        }
        no++;
    }
}

int ShapeTriangle::NShapeFunctions(int side, int order) {
    switch(side) {
        case 0:
        case 1:
        case 2:
            if(order == 1) return (1);
            else {std::cout << "Corner nodes shape are strictly defined at order 1"; DebugStop();}
        case 3:
        case 4:
        case 5:
        case 6:
            if(order == 2) return(1);
            else return (0);
        default:
            std::cout << "side out of bounds"; DebugStop();
    }
}

int ShapeTriangle::NShapeFunctions(VecInt &orders) {
    int fsides = 0, count = 0;
    for (int i: orders) {
        count += NShapeFunctions(fsides,i);
        fsides++;
    } return(count);
}
