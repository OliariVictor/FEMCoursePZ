/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "TopologyTetrahedron.h"
#include "tpanic.h"

static int nsidenodes[15] = {1,1,1,1,
                             2,2,2,2,2,2,
                             3,3,3,3,
                             4};

int TopologyTetrahedron::NSideNodes(int side) {
    if (side < 0 || side > TopologyTetrahedron::nSides) DebugStop();
    return(nsidenodes[side]);
}

int TopologyTetrahedron::SideNodeIndex(int side, int node) {
    if (side < 0 ) DebugStop();
    else if(side <4 && node == 0) return(side);
    else if(side <4) DebugStop();
    else if(side <6 && node <2)return((side+node)%4);
    else if(side ==6 && node <2) return((side+2*node)%4);
    else if(side <10 && node ==0)return(side-7);
    else if(side <10 && node ==1)return(3);
    else if(side <10) DebugStop();
    else if(side ==10 && node <3)return(node);
    else if(side <14 && node ==0 ||node ==1)return((side-11+node)%3);
    else if(side <14 && node ==2)return(3);
    else if(side <14) DebugStop();
    else if(side ==14 && node <4)return(node);
    DebugStop();
}

ElementType TopologyTetrahedron::Type() {
    return(ETetraedro);
}
