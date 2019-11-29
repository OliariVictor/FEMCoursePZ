/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "ShapeQuad.h"
#include "tpanic.h"

void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi) {
    int no = 0, counter = 0, size = NShapeFunctions(orders);
    phi.resize(size); dphi.Resize(size,2);

    auto f = [&](int i, int j, int k, int z) -> void {       //look for nodes 4-7
        phi[z] = (phi[i] + phi[j]) * phi[k];
        dphi(z, 0) = (phi[i] + phi[j]) * dphi(k, 0) + (dphi(i, 0) + dphi(j, 0)) * phi[k];
        dphi(z, 1) = (phi[i] + phi[j]) * dphi(k, 1) + (dphi(i, 1) + dphi(j, 1)) * phi[k];
    };

    Matrix phiBase(2,2,0), dphiBase(2,2,0);
    for(int i=0; i <2 ; i++) {
        phiBase(i,0) =  0.5*(1 - xi[i]) ; dphiBase(i,0) = -0.5;
        phiBase(i,1) =  0.5*(1 + xi[i]) ; dphiBase(i,1) = +0.5;
    }

    while (no < nSides) {
        if (orders[no] != 1 && orders[no] != 2) { std::cout << "Requested order not defined in this scope";DebugStop();}
        if (no < 4) {
            if (orders[no] == 1 || orders[no] == 2) {
                switch (no) {
                    case 0:
                        phi[0] = phiBase(0, 0) * phiBase(1, 0);
                        dphi(0, 0) = phiBase(1, 0) * dphiBase(0, 0);
                        dphi(0, 1) = phiBase(0, 0) * dphiBase(1, 0);
                        break;
                    case 1:
                        phi[1] = phiBase(0, 1) * phiBase(1, 0);
                        dphi(1, 0) = dphiBase(0, 1) * phiBase(1, 0);
                        dphi(1, 1) = phiBase(0, 1) * dphiBase(1, 0);
                        break;
                    case 2:
                        phi[2] = phiBase(0, 1) * phiBase(1, 1);
                        dphi(2, 0) = dphiBase(0, 1) * phiBase(1, 1);
                        dphi(2, 1) = phiBase(0, 1) * dphiBase(1, 1);
                        break;
                    case 3:
                        phi[3] = phiBase(0, 0) * phiBase(1, 1);
                        dphi(3, 0) = dphiBase(0, 0) * phiBase(1, 1);
                        dphi(3, 1) = phiBase(0, 0) * dphiBase(1, 1);
                        break;
                }
            } else std::cout << "Shape of corner nodes are defined only on first order";
        }else if(no < 8){
            if (orders[no] == 1) counter++;
            else if(orders[no] == 2) {
                switch(no){
                    case 4: f(0,3,1,4-counter);break;
                    case 5: f(1,0,2,5-counter);break;
                    case 6: f(2,1,3,6-counter);break;
                    case 7: f(3,2,0,7-counter);break;
                }
            } else {std::cout << "Please choose a valid shape order"; DebugStop;}
        }else if (no == 8){
            if (orders[no] == 2){
                phi[8-counter] = phi[0]*phi[2];
                dphi(8-counter,0) = dphi(0,0)*phi[2]+phi[0]*dphi(2,0);
                dphi(8-counter,1) = dphi(0,1)*phi[2]+phi[0]*dphi(2,1);
            }
        }
        no++;
    }
}

int ShapeQuad::NShapeFunctions(int side, int order) {
    switch(side) {
        case 0:
        case 1:
        case 2:
        case 3:
            if(order == 1 || order == 2) return (1);
            else {std::cout << "Corner nodes shape are strictly defined at order 1"; DebugStop();}
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
            if(order == 2) return(1);
            else return (0);
        default:
            std::cout << "side out of bounds"; DebugStop();
    }
}

int ShapeQuad::NShapeFunctions(VecInt &orders) {
    int fsides = 0, count = 0;
    for (int i: orders) {
        count += NShapeFunctions(fsides,i);
        fsides++;
    } return(count);
}
