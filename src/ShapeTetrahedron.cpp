/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "ShapeTetrahedron.h"
#include "tpanic.h"

void ShapeTetrahedron::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi) {
    int no = 4, counter = 0, size = NShapeFunctions(orders);
    phi.resize(size); dphi.Resize(size,3);

    phi[0] = 1 - xi[0] - xi[1] - xi[2];
    for(int i = 0; i<3 ; i++) phi[i+1] = xi[i];

    dphi(0,0) = dphi(0,1) = dphi(0,2) =  -1;
    for (int i = 1; i < 4; i++) for(int j = 0; j<3;j++) {
        if(i == j+1) dphi(i,j) = 1;
        else dphi(i,j) = 0;
    }

    auto edge = [&](int i, int j, int z) {
        phi[z] = phi[i]*phi[j];
        for (int k = 0; k < 3; k++) dphi(z,k) = dphi(i,k)*phi[j] + phi[i]*dphi(j,k);
    };

    auto face = [&](int i, int j, int k, int z){
        phi[z] = phi[i]*phi[j]*phi[k];
        for (int l = 0; l <3; l++) dphi(z,l) = dphi(i,l)*phi[j]*phi[k]+ phi[i]*dphi(j,l)*phi[k] + phi[i]*phi[j]*dphi(k,l);
    };

    auto volume = [&](int i, int j, int k, int l, int z){
        phi[z] = phi[i]*phi[j]*phi[k]*phi[l];
        for (int h = 0; h < 3 ; h++) dphi(z,h) = dphi(i,h)*phi[j]*phi[k]*phi[l]+phi[i]*dphi(j,h)*phi[k]*phi[l]+phi[i]*phi[j]*dphi(k,h)*phi[l]+phi[i]*phi[j]*phi[k]*dphi(l,h);
    };

    while (no < nSides) {
        if (orders[no] != 1 && orders[no] != 2) { std::cout << "Requested order not defined in this scope";DebugStop();}
        if(no < 10){
            if (orders[no] == 1) counter++;
            else if(orders[no] == 2) {
                switch(no){
                    case 4: edge(0,1,4-counter);break;
                    case 5: edge(1,2,5-counter);break;
                    case 6: edge(2,0,6-counter);break;
                    case 7: edge(3,0,7-counter);break;
                    case 8: edge(1,3,8-counter);break;
                    case 9: edge(2,3,9-counter);break;
                }
            } else {std::cout << "Please choose a valid shape order"; DebugStop;}
        }else if (no < 14){
            if (orders[no] == 1) counter++;
            else if(orders[no] == 2) {
                switch(no) {
                    case 10: face(0,1,2,10-counter);break;
                    case 11: face(0,1,3,11-counter);break;
                    case 12: face(1,2,3,12-counter);break;
                    case 13: face(0,2,3,13-counter);break;
                }
            }
        }else if(no == 14){
            volume(0,1,2,3,14-counter);break;
        }
        no++;
    }
}

int ShapeTetrahedron::NShapeFunctions(int side, int order) {
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
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
            if(order == 2) return(1);
            else return (0);
        default:
            std::cout << "side out of bounds"; DebugStop();
    }
}

int ShapeTetrahedron::NShapeFunctions(VecInt &orders) {
    int fsides = 0, count = 0;
    for (int i: orders) {
        count += NShapeFunctions(fsides,i);
        fsides++;
    } return(count);
}
