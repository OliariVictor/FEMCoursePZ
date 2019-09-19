/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "TopologyTriangle.h"
#include "tpanic.h"

static int nsidenodes[7]= {1,1,1,
                          2,2,2,
                          3};

int TopologyTriangle::NSideNodes(int side) {
    return(nsidenodes[side]);
}

int TopologyTriangle::SideNodeIndex(int side, int node) {
    if (side <0 || node <0) DebugStop();
    else if (side <3 && node ==0) return(side);
    else if (side <3) DebugStop();
    else if (side < 6 && node <2) return((side+node)%3);
    else if (side <6) DebugStop();
    else if (side ==6 && node <4) return(node);
    DebugStop();
}

ElementType TopologyTriangle::Type() {
    return(ETriangle);
}
