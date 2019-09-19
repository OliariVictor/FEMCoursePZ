/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Topology1d.h"
#include "tpanic.h"

int Topology1d::NSideNodes(int side) {
    switch(side) {
        case 0:
        case 1:
            return 1;
        case 2:
            return 2;
        default: {
            std::cout << " Topology1d::NSideNodes: 1D elements have only 3 nodes.\n";
            DebugStop();
        }
    }
}

int Topology1d::SideNodeIndex(int side, int node) {
    if(node >= NSideNodes(side)) DebugStop();
    switch(side) {
        case 0:
        case 1:
        return side;
    case 2:
        return node;
    default:
        std::cout << "Topology1d::SideNodeIndex: Invalid Side Number.\n";
        DebugStop;
    }
}

ElementType Topology1d::Type() {
    return EOned;
}
