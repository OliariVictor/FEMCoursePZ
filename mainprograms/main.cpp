

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
#include "Shape1d.h"
#include "IntRule1d.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "ShapeTetrahedron.h"
#include "GeoNode.h"
#include "GeoMesh.h"
#include "GeoElementTemplate.h"
#include "CompElementTemplate.h"
#include "CompElement.h"
#include "CompMesh.h"
#include "VTKGeoMesh.h"
#include "MathStatement.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "PostProcessTemplate.h"
#include "PostProcess.h"
#include "DOF.h"

using std::cout;
using std::endl;
using std::cin;
using namespace std;

void force2D(const VecDouble &co, VecDouble &result){
    result.resize(0); // nstate = 1
    if (co.size() != 2) {std::cout << "force2D: coordinate must have two dimensions"; DebugStop()};
    result[0] = sin(co[0])+cos(co[1]);
}

void lapace2D(const VecDouble &co, VecDouble &result, Matrix &deriv){
    result.resize(0); // nstate = 1
    Matrix.Resize(2,0);
    if (co.size() != 2) {std::cout << "force2D: coordinate must have two dimensions"; DebugStop()};
    result[0] = sin(co[0])+cos(co[1]);
    deriv(0,0) = cos(co[0]);
    deriv(1,0) = sin(co[1]);
}

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
    for(int i =0 ; i < TopologyQuad::nSides Shape1d teste;
    VecDouble coord = {0.0}, fphi(3,0); VecInt forders = {1,1,2};
    Matrix fdphi(1,3,0);
    teste.Shape(coord,forders,fphi,fdphi);
    cout << fphi[0] << "\t" << fphi[1] << "\t" << fphi.size() << endl;
    fdphi.Print();
    cout << teste.NShapeFunctions(forders);; ++i ) {
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

    //Testing IntRule::Constructors Shape1d teste;
    /*VecDouble coord = {0.0}, fphi(3,0); VecInt forders = {1,1,2};
    Matrix fdphi(1,3,0);
    teste.Shape(coord,forders,fphi,fdphi);
    cout << fphi[0] << "\t" << fphi[1] << "\t" << fphi.size() << endl;
    fdphi.Print();
    cout << teste.NShapeFunctions(forders);*/
    /*IntRule emptyConstructor;
    IntRule singularConstructor(30);
    IntRule copiedConstructor(singularConstructor);
    emptyConstructor.SetOrder(15);
    IntRule assinedConstructor = emptyConstructor;

    //Testing IntRule::NPoints
    std::cout << emptyConstructor.NPoints()<< std::endl << singularConstructor.NPoints() << std::endl << copiedConstructor.NPoints() << std::endl << assinedConstructor.NPoints() << std::endl;

    //Testing IntRuleinclude "tpanic.h"::GetOrder()
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
    emptyConstructorQuad.SetOrder(15)=;

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

    /*GeomTetrahedron tetra;
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
    tetra.SetNodes({4,5,6,7});VecInt no;tetra.GetNodes(no); cout <<no[3];*/

    /*Shape1d teste;
    VecDouble coord = {0.0}, fphi(3,0); VecInt forders = {1,1,2};
    Matrix fdphi(1,3,0);
    teste.Shape(coord,forders,fphi,fdphi);
    cout << fphi[0] << "\t" << fphi[1] << "\t" << fphi.size() << endl;
    fdphi.Print();
    cout << teste.NShapeFunctions(forders);*/

    /*ShapeQuad teste;
    VecDouble coord = {0.4,0.6}, fphi(teste.nSides,0); VecInt forders = {1,1,1,1,1,1,1,2,2};
    Matrix fdphi(teste.nSides,2,0);
    teste.Shape(coord,forders,fphi,fdphi);
    for (int i= 0 ; i < fphi.size();i++) cout << fphi[i] << "\t";cout << endl;
    fdphi.Print();
    //cout << teste.NShapeFunctions(forders);*/

    /*ShapeTriangle teste;
    VecDouble coord = {0.4,0.7}, fphi(teste.nSides,0); VecInt forders = {1,1,1,2,2,2,2};
    Matrix fdphi(teste.nSides,2,0);
    teste.Shape(coord,forders,fphi,fdphi);
    for (int i= 0 ; i < fphi.size();i++) cout << fphi[i] << "\t";cout << endl;
    fdphi.Print();*/

    /*ShapeTetrahedron teste;
    VecDouble coord = {0.4,0.7,0.2}, fphi(teste.nSides,0); VecInt forders = {1,1,1,1,2,2,2,2,2,2,2,2,2,2,2};
    Matrix fdphi(teste.nSides,3,0);
    teste.Shape(coord,forders,fphi,fdphi);
    for (int i= 0 ; i < fphi.size();i++) cout << fphi[i] << "\t";cout << endl;
    fdphi.Print();*/


    //Setting Geomesh;
    GeoMesh mesh;
    mesh.SetNumNodes(9);
    mesh.SetNumElements(4);
    mesh.SetDimension(2);

    //Defining node coordinates
    double index = 0;
    for(int i =0; i <mesh.NumNodes(); i +=3) {index = double(i/3); mesh.Node(i).SetCo({0,index*5});}
    for(int i =1; i <mesh.NumNodes(); i +=3) {index = double(i/3); mesh.Node(i).SetCo({index+5,index*5});}
    for(int i =2; i <mesh.NumNodes(); i +=3) {index = double(i/3); mesh.Node(i).SetCo({index+10,index*5});}

    VecInt oneIndices = {0,1,4,3},
           twoIndices = {1,2,5,4},
           threeIndices = {4,5,8,7},
           fourIndices = {3,4,7,6};

    GeoElementTemplate<GeomQuad> one(oneIndices,0,&mesh,0),
                                 two(twoIndices,1,&mesh,1),
                                 three(threeIndices,2,&mesh,2),
                                 four(fourIndices,3,&mesh,3);
    mesh.BuildConnectivity();

    mesh.Print(std::cout);

    Poisson math;
    CompMesh cmesh(&mesh);
    int numElem = cmesh.GetElementVec().size();
    for(int elemInd = 0; elemInd < numElem; elemInd++) {
        GeoElement *gel = mesh.Element(elemInd);
        int nSides = gel->NSides();
        for(int sideInd = 0; sideInd < nSides ; sideInd++){
            GeoElementSide *gelSide = new GeoElementSide(gel,sideInd);
            gel->
    }
    vector<DOF> dofVec1;



    cmesh.AutoBuild(); //Create the necessary CompElements && set firstEquation values.



    };
    //std::function<void(const VecDouble &co, VecDouble &result)> forceFunction;

    //std::function<void(const VecDouble &loc, VecDouble &result, Matrix &deriv)> SolutionExact;


    //VTKGeoMesh::PrintGMeshVTK(&mesh,"FirstVTK_1.vtk");
    //Computational Mesh

    /*CompMesh cmesh(&mesh);

    CompElementTemplate<ShapeQuad> cone(0,&cmesh, &one),
                                   ctwo(0,&cmesh, &two),
                                   cthree(0,&cmesh, &three),
                                   cfour(0,&cmesh, &four);

    vector<DOF> dofs;



    cmesh.SetMathVec()
    cmesh.SetDOFVec()
    //
    //
    //
    //
    //
    // cmesh.AutoBuild();
    //cmesh.Print(std::cout);*/




}




