

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
    result.resize(1); // nstate = 1
    if (co.size() != 2) {std::cout << "force2D: coordinate must have two dimensions"; DebugStop();};
    result[0] = sin(co[0])+cos(co[1]);
}

void laplace2D(const VecDouble &co, VecDouble &result, Matrix &deriv){
    result.resize(1); // nstate = 1
    deriv.Resize(2,1);
    if (co.size() != 2) {std::cout << "force2D: coordinate must have two dimensions"; DebugStop();};
    result[0] = sin(co[0])+cos(co[1]);
    deriv(0,0) = cos(co[0]);
    deriv(1,0) = -sin(co[1]);
}

void force3D(const VecDouble &co, VecDouble &result){
    result.resize(1); // nstate = 1
    if (co.size() != 3) {std::cout << "force3D: coordinates must have three dimensions"; DebugStop();};
    result[0] = sin(co[0])+cos(co[1]) + sin(co[2]);
}

void laplace3D(const VecDouble &co, VecDouble &result, Matrix &deriv){
    result.resize(1); // nstate = 1
    deriv.Resize(3,1);
    if (co.size() != 3) {std::cout << "force3D: coordinate must have two dimensions"; DebugStop();};
    result[0] = sin(co[0])+cos(co[1]) + sin(co[2]);
    deriv(0,0) = cos(co[0]);
    deriv(1,0) = -sin(co[1]);
    deriv(2,0) = cos(co[0]);
}

GeoMesh *CreateGmeshQuad(int nx);
GeoMesh *CreateGmeshTri(int nx);
GeoMesh *CreateGmeshTetra(int nx);

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
    /*GeoMesh *gmesh = new GeoMesh();
    gmesh->SetNumNodes(9);
    gmesh->SetNumElements(6);
    gmesh->SetDimension(2);

    //Defining node coordinates
    double index = 0;
    for(int i =0; i <gmesh->NumNodes(); i +=3) {index = double(i/3); gmesh->Node(i).SetCo({0,index*.5});}
    for(int i =1; i <gmesh->NumNodes(); i +=3) {index = double(i/3); gmesh->Node(i).SetCo({0.5,index*.5});}
    for(int i =2; i <gmesh->NumNodes(); i +=3) {index = double(i/3); gmesh->Node(i).SetCo({1,index*.5});}

    VecInt oneIndices = {0,1,4,3},
           twoIndices = {1,2,5,4},
           threeIndices = {4,5,8,7},
           fourIndices = {3,4,7,6},
           fiveIndices = {0,1},
           sixIndices = {1,2},
           sevenIndices = {2,5},
           eightIndices = {5,8},
           nineIndices = {8,7},
           tenIndices ={7,6},
           elevenIndices ={6,3},
           twelveIndices ={3,0};

    GeoElementTemplate<GeomQuad> one(oneIndices,1,gmesh,0),
                                 two(twoIndices,1,gmesh,1),
                                 three(threeIndices,1,gmesh,2),
                                 four(fourIndices,1,gmesh,3);

    GeoElementTemplate<Geom1d> five(fiveIndices,-1,gmesh,4),
                               six(sixIndices,-1,gmesh,5);
                               //seven(sevenIndices,-1,gmesh,6),
                               //eight(eightIndices,-1,gmesh,7),
                               //nine(nineIndices,-1,gmesh,8),
                               //ten(tenIndices,-1,gmesh,9),
                               //eleven(elevenIndices,-1,gmesh,10),
                               //twelve(twelveIndices,-1,gmesh,11);

    gmesh->BuildConnectivity();*/

    //GeoMesh *gmesh = CreateGmeshQuad(8);
    //GeoMesh *gmesh = CreateGmeshTri(2);

    //GeoElementTemplate<GeomTetrahedron> geo();
    GeoMesh *gmesh = CreateGmeshTetra(2);


    gmesh->Print(std::cout); VTKGeoMesh::PrintGMeshVTK(gmesh,"GMESH.vtk");


    CompMesh *cmesh = new CompMesh(gmesh);
    int numElem = gmesh->NumElements();
    int poissonId = 1;
    Matrix perm(0,0), proj(0,0);
    cmesh->SetNumberMath(numElem);
    cmesh->SetNumberElement(numElem);

    for(int elemInd = 0; elemInd < numElem ; elemInd++){
        //Set cmesh's statement...
        GeoElement *gel = gmesh->Element(elemInd);
        int matId = gel->Material(), matDim = -1;

        switch(gel-> Type()){
            case EOned: matDim =1; break;
            case ETriangle: case EQuadrilateral: matDim = 2; break;
            case ETetraedro: matDim = 3; break;
            default: std::cout << "Invalid Element Type"; DebugStop();
        }

        if(matId == poissonId){    //Interior Element
            Poisson *mat = new Poisson(poissonId,perm);
            mat->SetDimension(matDim);
            if(gmesh->Dimension() == 2) { mat->SetForceFunction(force2D); mat->SetExactSolution(laplace2D);}
            else {mat->SetForceFunction(force3D); mat->SetExactSolution(laplace3D);}
            cmesh->SetMathStatement(elemInd,mat);
        }
        else{    //Outer Element
            //L2Projection(int bctype, int materialid, Matrix &proj, Matrix Val1, Matrix Val2);
            Matrix Val1(1,1),Val2(1,1);
            L2Projection *mat = new L2Projection(0,matId,proj,Val1,Val2);
            mat->SetDimension(matDim);
            if(gmesh->Dimension() == 2) {mat->SetForceFunction(force2D); mat->SetExactSolution(laplace2D);}
            else {mat->SetForceFunction(force3D); mat->SetExactSolution(laplace3D);}
            cmesh->SetMathStatement(elemInd,mat);
        }
    }
    cmesh->SetDefaultOrder(1);
    cmesh->AutoBuild();
    cmesh->Print(std::cout);

    Analysis analysis(cmesh);
    analysis.RunSimulation();
    VecDouble Sol = cmesh->Solution();

    PostProcessTemplate<Poisson> process;
    if(gmesh->Dimension() == 2) process.SetExact(laplace2D); else process.SetExact(laplace3D);
    analysis.PostProcessError(std::cout,process);
    /*Matrix nodes(9,9,0);
    nodes(0,0) = 0; nodes(0,1)= 0;
    nodes(1,0) = 0.5; nodes(1,1)= 0.;
    nodes(2,0) = 0.5; nodes(2,1)= 0.5;
    nodes(3,0) = 0.; nodes(3,1)= 0.5;
    nodes(4,0) = 1.; nodes(4,1)= 0.;
    nodes(5,0) = 1.; nodes(5,1)= 0.5;
    nodes(6,0) = 0.5; nodes(6,1)= 1.;
    nodes(7,0) = 0.; nodes(7,1)= 1.;
    nodes(8,0) = 1.; nodes(8,1)= 1.;
    vector<VecDouble> exact(9); Matrix vide(0,0);
    for(int ni = 0; ni < nodes.Rows(); ni++){
        VecDouble coord = {nodes(ni,0),nodes(ni,1)};
        laplace2D(coord, reinterpret_cast<VecDouble &>(exact[ni]), vide); std::cout << "\n exact:   " << exact[ni][0];
    }*/


    process.AppendVariable("Sol");
    process.AppendVariable("Sol_Exact");
    process.AppendVariable("Force");

    VTKGeoMesh::PrintGMeshVTK(gmesh,"GMESH.vtk");
    VTKGeoMesh::PrintCMeshVTK(cmesh,2,"CMESH.vtk"); //
    //analysis.PostProcessSolution("SOLUTION.vtk",process);
}

GeoMesh *CreateGmeshQuad(int nx ){
    GeoMesh *gmesh = new GeoMesh();
    int numNodes = (nx+1)*(nx+1);
    gmesh->SetNumNodes(numNodes);
    int numElem = nx*nx+4*nx;
    gmesh->SetNumElements(numElem);
    gmesh->SetDimension(2);

    //Defining node coordinates
    double h = 1./nx;
    for(int i =0; i < nx+1;i++) for(int j = 0; j < nx+1; j++) gmesh->Node(i*(nx+1)+j).SetCo({h*j,h*i});

    VecInt indices(4,0);int index = 0;
    for(int xInd = 0; xInd < nx; xInd++) for(int yInd = 0; yInd < nx ; yInd++){
        indices[0] = (nx+1)*(xInd) + yInd;
        indices[1] = (nx+1)*(xInd) + yInd+1;
        indices[3] = (nx+1)*(xInd+1) + yInd;
        indices[2] = (nx+1)*(xInd+1) + yInd+1;
        GeoElementTemplate<GeomQuad> *quad = new GeoElementTemplate<GeomQuad>(indices,1,gmesh,index);
        gmesh->SetElement(index,quad);
        index++;
    }
    //Bottom
    VecInt indicesBC(nx,0);
    for(int BCind = 0; BCind < nx; BCind++){
        indicesBC[0] = BCind;
        indicesBC[1] = BCind+1;
        GeoElementTemplate<Geom1d> *BC = new GeoElementTemplate<Geom1d>(indicesBC,-1,gmesh,index);
        gmesh->SetElement(index,BC);
        index++;
    }
    //Right
    for(int BCind = 0; BCind < nx; BCind++){
        indicesBC[0] = nx+(nx+1)*BCind;
        indicesBC[1] = indicesBC[0]+(nx+1);
        GeoElementTemplate<Geom1d> *BC = new GeoElementTemplate<Geom1d>(indicesBC,-1,gmesh,index);
        gmesh->SetElement(index,BC);
        index++;
    }
    //TOP
    for(int BCind = 0; BCind < nx; BCind++){
        indicesBC[0] = (numNodes-1)-BCind;
        indicesBC[1] = indicesBC[0]-1;
        GeoElementTemplate<Geom1d> *BC = new GeoElementTemplate<Geom1d>(indicesBC,-1,gmesh,index);
        gmesh->SetElement(index,BC);
        index++;
    }
    //Left
    for(int BCind = 0; BCind < nx; BCind++){
        indicesBC[0] = numNodes-(nx+1)*(1+BCind);
        indicesBC[1] = indicesBC[0]-(nx+1);
        GeoElementTemplate<Geom1d> *BC = new GeoElementTemplate<Geom1d>(indicesBC,-1,gmesh,index);
        gmesh->SetElement(index,BC);
        index++;
    }

    gmesh->BuildConnectivity();
    return(gmesh);
}

GeoMesh *CreateGmeshTri(int nx ){
    GeoMesh *gmesh = new GeoMesh();
    int numNodes = (nx+1)*(nx+1);
    gmesh->SetNumNodes(numNodes);
    int numElem = 2*nx*nx+4*nx;
    gmesh->SetNumElements(numElem);
    gmesh->SetDimension(2);

    //Defining node coordinates
    double h = 1./nx;
    for(int i =0; i < nx+1;i++) for(int j = 0; j < nx+1; j++) gmesh->Node(i*(nx+1)+j).SetCo({h*j,h*i});

    //Volume Elements
    VecInt L(3,0),U(3,0);
    int indL = 0, indU = 1;
    for(int j = 0; j < nx; j++) for(int i = 0; i < nx ; i++){
            L[0] =  i + (nx+1)*j;
            L[1] = L[0] + 1;
            L[2] = L[0]+(nx+1);
            GeoElementTemplate<GeomTriangle> *lower = new GeoElementTemplate<GeomTriangle>(L,1,gmesh,indL);
            gmesh->SetElement(indL,lower);
            indL+=2;

            U[0] = L[2]+1;
            U[1] = L[2];
            U[2] = L[1];
            GeoElementTemplate<GeomTriangle> *upper = new GeoElementTemplate<GeomTriangle>(U,1,gmesh,indU);
            gmesh->SetElement(indU,upper);
            indU+=2;
        }
    //Bottom
    int index = indL;
    VecInt indicesBC(nx,0);
    for(int BCind = 0; BCind < nx; BCind++){
        indicesBC[0] = BCind;
        indicesBC[1] = BCind+1;
        GeoElementTemplate<Geom1d> *BC = new GeoElementTemplate<Geom1d>(indicesBC,-1,gmesh,index);
        gmesh->SetElement(index,BC);
        index++;
    }
    //Right
    for(int BCind = 0; BCind < nx; BCind++){
        indicesBC[0] = nx+(nx+1)*BCind;
        indicesBC[1] = indicesBC[0]+(nx+1);
        GeoElementTemplate<Geom1d> *BC = new GeoElementTemplate<Geom1d>(indicesBC,-1,gmesh,index);
        gmesh->SetElement(index,BC);
        index++;
    }
    //TOP
    for(int BCind = 0; BCind < nx; BCind++){
        indicesBC[0] = (numNodes-1)-BCind;
        indicesBC[1] = indicesBC[0]-1;
        GeoElementTemplate<Geom1d> *BC = new GeoElementTemplate<Geom1d>(indicesBC,-1,gmesh,index);
        gmesh->SetElement(index,BC);
        index++;
    }
    //Left
    for(int BCind = 0; BCind < nx; BCind++){
        indicesBC[0] = numNodes-(nx+1)*(1+BCind);
        indicesBC[1] = indicesBC[0]-(nx+1);
        GeoElementTemplate<Geom1d> *BC = new GeoElementTemplate<Geom1d>(indicesBC,-1,gmesh,index);
        gmesh->SetElement(index,BC);
        index++;
    }

    gmesh->BuildConnectivity();
    return(gmesh);
}

GeoMesh *CreateGmeshTetra(int nx){
    GeoMesh *gmesh = new GeoMesh();
    int numNodes = (nx+1)*(nx+1)*(nx+1);
    gmesh->SetNumNodes(numNodes);
    int numElem = 6*nx*nx*nx +6*nx*nx;
    gmesh->SetNumElements(numElem);
    gmesh->SetDimension(3);

    //Defining node coordinates
    double h = 1./nx;
    for(int k =0; k < nx+1;k++) for(int j = 0; j < nx+1; j++) for(int i = 0; i < nx+1; i++) gmesh->Node(k*(nx+1)*(nx+1)+j*(nx+1)+i).SetCo({h*i,h*j,h*k});

    //Volume Elements
    VecInt T(4,0);

    int index = 0;
    int kp = (nx+1)*(nx+1), jp = (nx+1), ip = 1;
    for(int k =0; k < nx;k++) for(int j = 0; j < nx; j++) for(int i = 0; i < nx ; i++){
            int comum = k*(nx+1)*(nx+1)+j*(nx+1)+i;

            T[0] = comum+ kp;
            T[1] = comum;
            T[2] = comum+ip+kp;
            T[3] = comum+ip + jp+kp;
            GeoElementTemplate<GeomTetrahedron> *tet1 = new GeoElementTemplate<GeomTetrahedron>(T,1,gmesh,index);
            gmesh->SetElement(index,tet1);
            index++;

            T[0] = comum+kp;
            T[1] = comum;
            T[2] = comum+ jp+kp;
            T[3] = comum+ip+jp+kp;
            GeoElementTemplate<GeomTetrahedron> *tet2 = new GeoElementTemplate<GeomTetrahedron>(T,1,gmesh,index);
            gmesh->SetElement(index,tet2);
            index++;

            T[0] = comum;
            T[1] = comum+jp;
            T[2] = comum + jp + kp;
            T[3] = comum + jp + kp + ip;
            GeoElementTemplate<GeomTetrahedron> *tet3 = new GeoElementTemplate<GeomTetrahedron>(T,1,gmesh,index);
            gmesh->SetElement(index,tet3);
            index++;

            T[0] = comum + ip;
            T[1] = comum;
            T[2] = comum +ip + kp;
            T[3] = comum + ip + kp + jp;
            GeoElementTemplate<GeomTetrahedron> *tet4 = new GeoElementTemplate<GeomTetrahedron>(T,1,gmesh,index);
            gmesh->SetElement(index,tet4);
            index++;

            T[0] = comum + ip + jp;
            T[1] = comum;
            T[2] = comum + ip;
            T[3] = comum + ip + jp + kp;
            GeoElementTemplate<GeomTetrahedron> *tet5 = new GeoElementTemplate<GeomTetrahedron>(T,1,gmesh,index);
            gmesh->SetElement(index,tet5);
            index++;

            T[0] = comum;
            T[1] = comum + jp;
            T[2] = comum + jp + ip;
            T[3] = comum + jp+ip+kp;
            GeoElementTemplate<GeomTetrahedron> *tet6 = new GeoElementTemplate<GeomTetrahedron>(T,1,gmesh,index);
            gmesh->SetElement(index,tet6);
            index++;
        }
    //BOUNDARY CONDITION
    //XY BOTTOM
    VecInt LBC(3,0),UBC(3,0);
    for(int j = 0; j < nx; j++) for(int i = 0; i < nx ; i++){
            int Lcomum = i + (nx+1)*j;
            int ip = 1, jp = nx+1;

            LBC[0] = Lcomum;
            LBC[1] = Lcomum + ip;
            LBC[2] = Lcomum+jp;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + ip;

            UBC[0] = Ucomum+jp;
            UBC[1] = Ucomum;
            UBC[2] = Ucomum+jp-ip;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index;
        }
    //XY TOP
    for(int j = 0; j < nx; j++) for(int i = 0; i < nx ; i++){
            int Lcomum = i + (nx+1)*j+nx*(nx+1)*(nx+1);
            int ip = 1, jp = nx+1;

            LBC[0] = Lcomum;
            LBC[1] = Lcomum + ip;
            LBC[2] = Lcomum+jp;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + ip;

            UBC[0] = Ucomum+jp;
            UBC[1] = Ucomum;
            UBC[2] = Ucomum+jp-ip;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index;
        }
    //ZX FRONT
    for(int k = 0; k < nx; k++) for(int i = 0; i < nx ; i++){
            int Lcomum = i + (nx+1)*(nx+1)*k;
            int ip = 1, kp = (nx+1)*(nx+1);

            LBC[0] = Lcomum;
            LBC[1] = Lcomum + ip;
            LBC[2] = Lcomum+kp;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + ip;

            UBC[0] = Ucomum+kp;
            UBC[1] = Ucomum;
            UBC[2] = Ucomum+kp-ip;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index;
        }
    //ZX BACK
    for(int k = 0; k < nx; k++) for(int i = 0; i < nx ; i++){
            int Lcomum = i + (nx+1)*(nx+1)*k+nx*(nx+1);
            int ip = 1, kp = (nx+1)*(nx+1);

            LBC[0] = Lcomum;
            LBC[1] = Lcomum + ip;
            LBC[2] = Lcomum+kp;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + ip;

            UBC[0] = Ucomum+kp;
            UBC[1] = Ucomum;
            UBC[2] = Ucomum+kp-ip;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index;
        }
    //YZ LEFT
    for(int k = 0; k < nx; k++) for(int j = 0; j < nx ; j++){
            int Lcomum =  j*(nx+1)+k*(nx+1)*(nx+1);
            int jp = (nx+1), kp = (nx+1)*(nx+1);

            LBC[0] = Lcomum;
            LBC[1] = Lcomum + jp;
            LBC[2] = Lcomum+kp;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + jp;

            UBC[0] = Ucomum+kp;
            UBC[1] = Ucomum;
            UBC[2] = Ucomum+kp-jp;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index;
        }
    //YZ RIGHT
    for(int k = 0; k < nx; k++) for(int j = 0; j < nx ; j++){
            int Lcomum =  j*(nx+1)+k*(nx+1)*(nx+1) +nx*1;
            int jp = (nx+1), kp = (nx+1)*(nx+1);

            LBC[0] = Lcomum;
            LBC[1] = Lcomum + jp;
            LBC[2] = Lcomum+kp;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + jp;

            UBC[0] = Ucomum+kp;
            UBC[1] = Ucomum;
            UBC[2] = Ucomum+kp-jp;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index;
        }
    gmesh->BuildConnectivity();
    return(gmesh);
}