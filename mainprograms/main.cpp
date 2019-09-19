

#include <iostream>
#include <cmath.h>
#include <TMatrix.h>
#include <Topology1d.h>
#include <TopologyQuad.h>
#include <TopologyTetrahedron.h>
#include <TopologyTriangle.h>

using std::cout;
using std::endl;
using std::cin;

int main ()
{
    //testing Topology 1D
    /*int side = 2, node = 1;
    int Nsides = Topology1d::NSideNodes(side);
    int localNode = Topology1d::SideNodeIndex(side,node);
    int type = Topology1d::Type();
    cout <<"Local node: "<< type;*/

    //testing TopologyQuad.h
    /*std::cout << "{,";
    for(int i =0 ; i < TopologyQuad::nSides; ++i ) {
        std::cout <<" "<< TopologyQuad::NSideNodes(i) << ",";
    };
    std::cout << "}\n";
    int localNode = TopologyQuad::SideNodeIndex(2,0);
    std::cout << localNode << " " << EQuadrilateral;*/


    //testing TopologyTetrahedron
    /*int Nsides = TopologyTetrahedron::NSideNodes(10);
    int node = TopologyTetrahedron::SideNodeIndex(11,2);
    std::cout << Nsides << "\n" << node;*/

    int Nsides = TopologyTriangle::NSideNodes(6);
    int node = TopologyTriangle::SideNodeIndex(3,0);
    std::cout << Nsides << "\n" << node;

    return 0;
}
