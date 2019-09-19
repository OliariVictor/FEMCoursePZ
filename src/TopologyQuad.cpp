/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "TopologyQuad.h"
#include "tpanic.h"

static int nsidenodes[TopologyQuad::nSides] = {1,1,1,1,2,2,2,2,4};

int TopologyQuad::NSideNodes(int side) {
    if(side >= TopologyQuad::nSides) DebugStop();
    else return nsidenodes[side];
}

int TopologyQuad::SideNodeIndex(int side, int node) {
    if(side >= TopologyQuad::nSides || side <0) DebugStop;
    if(node <0) DebugStop;
    if(side <4) {
        if (node == 0) return side;
        else DebugStop();
    }
    else if(side < 8) {
        if (node == 0) return (side - 4);
        else if (node == 1 && side <7 ) return (side - 3);
        else if (node == 1 && side == 7 ) return (0);
        else DebugStop();
    }

    else if (side ==8) {
        if (node < 4) return node;
        else DebugStop();
    }
    else DebugStop();

}

ElementType TopologyQuad::Type() {
    return EQuadrilateral;
}
