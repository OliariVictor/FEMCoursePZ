

#include <iostream>
#include <TMatrix.h>
#include "Topology1d.h"
#include <TopologyQuad.h>
#include <TopologyTetrahedron.h>
#include <TopologyTriangle.h>
#include "DataTypes.h"
//#include "IntRule1d.h"

using std::cout;
using std::endl;
using std::cin;
using namespace std;

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
    int localNode = TopologyQuad::SideNodeIndex(2,https://dev.to/sandordargo/lambda-expressions-in-c-4pj40);
    std::cout << localNode << " " << EQuadrilateral;*/


    //testing TopologyTetrahedron
    /*int Nsides = TopologyTetrahedron::NSideNodes(10);
    int node = TopologyTetrahedron::SideNodeIndex(11,2);
    std::cout << Nsides << "\n" << node;*/

    //Testing TopologyTriangle
    /*int Nsides = TopologyTriangle::NSideNodes(6);
    int node = TopologyTriangle::SideNodeIndex(3,0);
    std::cout << Nsides << "\n" < 0.57735 || weight =1

Process finished with exit code 0< node;*/

    //Testing IntRule::Constructors
    /*IntRule emptyConstructor;
    IntRule singularConstructor(30);
    IntRule copiedConstructor(singularConstructor);
    emptyConstructor.SetOrder(15);
    IntRule assinedConstructor = emptyConstructor;

    //Testing IntRule::NPoints
    std::cout << emptyConstructor.NPoints()<< std::endl << singularConstructor.NPoints() << std::endl << copiedConstructor.NPoints() << std::endl << assinedConstructor.NPoints() << std::endl;

    //Testing IntRule::GetOrder()
    std::cout << emptyConstructor.GetOrder()<< std::endl << singularConstructor.GetOrder() << std::endl << copiedConstructor.GetOrder() << std::endl << assinedConstructor.GetOrder() << std::endl;

    //singularConstructor.Print(std::cout);*/

    //Testing IntRule1d
    /*IntRule1d emptyConstructor1d;
    IntRule1d singularConstructor1d(4);
    emptyConstructor1d.SetOrder(5);

    std::cout << emptyConstructor1d.NPoints();
    singularConstructor1d.Print(std::cout);

    int IntNumber = singularConstructor1d.NPoints();
    VecDouble coord(IntNumber),weight(IntNumber);
    singularConstructor1d.gauleg(-1.,1.,coord,weight);

    int i;
    for (i =0; i< coord.size(); i++){
        std::cout <<"Coordinate: "<< coord[i] << " || Weight: " << weight[i]<< endl;
    }*/

    //Testing InRuleQuadr
    /*IntRuleQuad emptyConstructorQuad;
    IntRuleQuad singularConstructorQuad(4);
    emptyConstructorQuad.SetOrder(15);

    emptyConstructorQuad.Print(std::cout);*/

    //Testing InRuleTriang
    /*IntRuleTriangle emptyConstructorTri;
    IntRuleTriangle singularConstructorTri(4);
    emptyConstructorTri.SetOrder(12);

    emptyConstructorTri.Print(std::cout);*/

    //Testing InRuleTetra
    IntRuleTetrahedron emptyConstructorTet;
    IntRuleTetrahedron singularConstructorTet(4);
    emptyConstructorTet.SetOrder(12);

    emptyConstructorTet.Print(std::cout);
}
