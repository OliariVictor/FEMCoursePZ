/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Shape1d.h"
#include "tpanic.h"

void Shape1d::Shape(const VecDouble& xi, VecInt& orders, VecDouble& phi, Matrix& dphi) {
    int no = 0;
    for(int i :orders){
        if (orders[no] != 1 && orders[no] != 2) { std::cout << "Requested order not defined in this scope";DebugStop();}
        if (no < 2){
            if (orders[no] == 1 || orders[no] == 2) {
                phi[0] = 0.5*(1 - xi[0]); dphi(0,0) = -0.5;
                phi[1] = 0.5*(1 + xi[0]); dphi(1,0) = +0.5;
            } else {std::cout << "Shape functions on corner nodes are strictly linear"; DebugStop();}
        }
        else if (no == 2) {
            if (orders[no] == 2 ) {phi[2] = phi[0]*phi[1]; dphi(2,0) = dphi(0,0)*phi[1] + dphi(1,0)*phi[0]; }
            else { int size = NShapeFunctions(orders); phi.resize(size);dphi.Resize(size,1);}
        } no++;
      }
}

int Shape1d::NShapeFunctions(int side, int order) {
    switch(side) {
        case 0:
        case 1:
            if(order == 1 || order == 2) return (1);
            else {std::cout << "Corner nodes are strictly defined at order 1"; DebugStop();}
        case 2:
            if(order == 2) return(1);
            else return (0);
    }
}

int Shape1d::NShapeFunctions(VecInt &orders) {
    int fsides = 0, count = 0;
    for (int i: orders) {
        count += NShapeFunctions(fsides,i);
        fsides++;
    } return(count);
}
