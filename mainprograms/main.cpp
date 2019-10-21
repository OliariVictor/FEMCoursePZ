

#include <iostream>
#include <TMatrix.h>
#include "Topology1d.h"
#include <TopologyQuad.h>
#include <TopologyTetrahedron.h>
#include <TopologyTriangle.h>
#include <Geom1d.h>
#include <GeomQuad.h>
#include <GeomTriangle.h>
#include <GeomTetrahedron.h>
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
        std::cout <<" "<< TopologyQuad::NSideut << x[0] << std::endl;Nodes(i) << ",";
    };
    std::cout << "}\n";
    int localNode = TopologyQuad::SideNodeIndex(2,https://dev.to/sandordargo/lambda-expressions-in-c-4pj40);
    std::cout << localNode << " " << EQuadrilateral;*/


    //testing TopologyTetrahedron
    /*int NsidesTriangle = TopologyTetrahedron::NSideNodes(10);
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
IntRule emptyConstructor;
    IntRule singularConstructor(30);
    IntRule copiedConstructor(singularConstructor);
    emptyConstructor.SetOrder(15);
    IntRule assinedConstructor = emptyConstructor;

    //Testing IntRule::NPoints
    std::cout << emptyConstructor.NPoints()<< std::endl << singularConstructor.NPoints() << std::endl << copiedConstructor.NPoints() << std::endl << assinedConstructor.NPoints() << std::endl;

    //Testing IntRule::GetOrder()
    std::cout << emptyConstructor.GetOrder()<< std::endl << singularConstructor.GetOrder() << std::endl << copiedConstructor.GetOrder() << std::endl << assinedConstructor.GetOrder() << std::endl;

    //singularConstructor.Print(std::cout);
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
/*//tetrahedral element (pztopology::TPZTetrahedra) has 15 sides (4 points 6 lines 4 triangles and one volume)
//point element (pztopology::TPZPoint) has one side: the point itself
    //Testing InRuleQuadr
    /*IntRuleQuad emptyConstructorQuad;
    IntRuleQuad singularConstructorQuad(4);
    emptyConstructorQuad.SetOrder(15);

    emptyConstructorQuad.Print(std::cout);*/

    //Testing InRuleTriangut << x[0] << std::endl;
    /*IntRuleTriangle emptyConstructorTri;
    IntRuleTriangle singularConstructorTri(4);
    emptyConstructorTri.SetOrder(12);

    emptyConstructorTri.Print(std::cout);

    //Testing InRuleTetra
    IntRuleTetrahedron emptyConstructorTet;
    IntRuleTetraheTriangledron singularConstructorTet(4);
    emptyConstructorTet.SetOrder(12);

    emptyConstructorTet.Print(std::cout);*/

    //Testing Geom1d
    /*Geom1d line;
    VecDouble coord = {0}, pos = {0.0};
    Matrix nodes(1,2,0),gradx(1,1); nodes(0,0) = 15; nodes(0,1) = 22;
    line.X(coord, nodes, pos);
    line.GradX(coord, nodes, pos, gradx);
    //cout << gradx(0,0);

    VecInt kkk = {34,35};
    line.SetNodes(kkk);
    int num = line.NodeIndex(1);
    //cout << line.NumNodes();*/

    /*GeomQuad quad;
    VecDouble coord = {-1,0.8}, pos = {0,0};
    Matrix nodes(4,2,0), grad(2,2,0); nodes(0,0) = nodes(3,0) = 0;
    nodes(1,0) = nodes(2,0) = 10;
    nodes(0,1) = nodes(1,1) = 0;
    nodes(2,1) = nodes(3,1) = 10;
    quad.X(coord, nodes, pos);
    quad.GradX(coord, nodes, pos,grad);
    cout << pos[0] << endl;  cout << pos[1] << endl;
    cout << grad(0,0) << "\t" << grad(0,1) << endl;
    cout << grad(1,0) << "\t" << grad(1,1) << endl;*/

    /*GeomTriangle tri;
    VecDouble coord = {0.5,0.5}, pos = {0,0};
    Matrix nodes(3,2,0), grad(2,2,0);
    nodes(0,0) = nodes(0,1) = nodes(1,1) = nodes(2,0)=0;
    nodes(1,0) = nodes(2,1) = 10;
    tri.X(coord, nodes, pos);
    tri.GradX(coord, nodes, pos,grad);
    cout << pos[0] << "\t";  cout << pos[1] << endl;
    cout << grad(0,0) << "\t" << grad(0,1) << endl;
    cout << grad(1,0) << "\t" << grad(1,1) << endl;
    tri.SetNodes({4,5,6});VecInt no;tri.GetNodes(no); cout <<no[2];*/

    GeomTetrahedron tetra;
    VecDouble coord = {0.2,0.5,0.2}, pos = {0,0,0};
    Matrix nodes(4,3,0), grad(3,3,0), fdphi(4,3,0); VecDouble fphi(4,0);
    nodes(1,0) = nodes(2,1) = nodes(3,2) =  15; nodes(0,0) = nodes(0,1) = nodes(0,2) = -1; nodes(0,2) = -3;
    //tetra.Shape(coord, fphi, fdphi) ; for(int i = 0; i < 4; i++) cout << fphi[i] << "\t"; cout <<endl; for(int i = 0; i < 4; i++) cout << nodes(i,0) << "\t"; cout <<endl;
    tetra.X(coord, nodes, pos);
    tetra.GradX(coord, nodes, pos,grad);
    cout << pos[0] << "\t";  cout << pos[1] << "\t";   cout << pos[2] << endl;
    cout << grad(0,0) << "\t" << grad(0,1) << "\t" << grad(0,2) << endl;
    cout << grad(1,0) << "\t" << grad(1,1) << "\t" << grad(1,2) << endl;
    cout << grad(2,0) << "\t" << grad(2,1) << "\t" << grad(2,2) << endl;
    tetra.SetNodes({4,5,6,7});VecInt no;tetra.GetNodes(no); cout <<no[3];
}
