

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
    deriv(2,0) = cos(co[2]);
}

void laplace3D2P(const VecDouble &co, VecDouble &result, Matrix &deriv){
    result.resize(1);
    deriv.Resize(3,1);
    if (co.size() != 3) {std::cout << "force3D: coordinate must have two dimensions"; DebugStop();};
    result[0] = 0.5*(co[0]*co[0]+co[1]*co[1]+co[2]*co[2]);//co[0]+co[1]+co[2];//
    deriv(0,0) = co[0];//1;
    deriv(1,0) = co[1];//1;
    deriv(2,0) = co[2];//1;
}

void force3D2P(const VecDouble &co, VecDouble &result){
    result.resize(1);
    if (co.size() != 3) {std::cout << "force3D: coordinates must have three dimensions"; DebugStop();};
    result[0] = -3;//-3;//0;
}

void laplace3D1P(const VecDouble &co, VecDouble &result, Matrix &deriv){
    result.resize(1);
    deriv.Resize(3,1);
    if (co.size() != 3) {std::cout << "force3D: coordinate must have two dimensions"; DebugStop();};
    result[0] = (co[0]+co[1]+co[2]);
    deriv(0,0) = 1;
    deriv(1,0) = 1;
    deriv(2,0) = 1;
}

void force3D1P(const VecDouble &co, VecDouble &result){
    result.resize(1);
    if (co.size() != 3) {std::cout << "force3D: coordinates must have three dimensions"; DebugStop();};
    result[0] = 0;
}

void force2D2P(const VecDouble &co, VecDouble &result){
    result.resize(1);
    if (co.size() != 2) {std::cout << "force2D: coordinate must have two dimensions"; DebugStop();};
    result[0] = -2;
}

void laplace2D2P(const VecDouble &co, VecDouble &result, Matrix &deriv){
    result.resize(1);
    deriv.Resize(2,1);
    if (co.size() != 2) {std::cout << "force2D: coordinate must have two dimensions"; DebugStop();};
    result[0] = 0.5*(co[0]*co[0]+co[1]*co[1]);
    deriv(0,0) = co[0];
    deriv(1,0) = co[1];
}

void force2D1P(const VecDouble &co, VecDouble &result){
    result.resize(1);
    if (co.size() != 2) {std::cout << "force2D: coordinate must have two dimensions"; DebugStop();};
    result[0] = 0;
}

void laplace2D1P(const VecDouble &co, VecDouble &result, Matrix &deriv){
    result.resize(1);
    deriv.Resize(2,1);
    if (co.size() != 2) {std::cout << "force2D: coordinate must have two dimensions"; DebugStop();};
    result[0] = (co[0]+co[1]);
    deriv(0,0) = 1;
    deriv(1,0) = 1;
}

void setFuncPoisson(const int dim, const int fType, Poisson *pot);
void setFuncL2(const int dim, const int fType, L2Projection *l2);
void setFuncPostProcess(const int dim, const int fType, PostProcessTemplate<Poisson> *pp);

GeoMesh *CreateGmeshQuad(int nx);
GeoMesh *CreateGmeshTri(int nx);
GeoMesh *CreateGmeshTetra(int nx);

int main ()
{
    int pOrder = 2; // 1 or 2
    int solType = 0; // 0: Laplace Problem ; 1: Sol exact for pOrder == 1; 2: Sol Exact for pOrder == 2;

    //GeoMesh *gmesh = CreateGmeshQuad(1);
    //GeoMesh *gmesh = CreateGmeshTri(8);
    GeoMesh *gmesh = CreateGmeshTetra(8);

    //gmesh->Print(std::cout);
    std::cout << "\ngeometric mesh executed\n";

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
            setFuncPoisson(gmesh->Dimension(),solType,mat);
            cmesh->SetMathStatement(elemInd,mat);
        }
        else{    //Outer Element
            Matrix Val1(1,1),Val2(1,1);
            L2Projection *mat = new L2Projection(0,matId,proj,Val1,Val2);
            mat->SetDimension(matDim);
            setFuncL2(gmesh->Dimension(),solType,mat);
            cmesh->SetMathStatement(elemInd,mat);
        }
    }
    cmesh->SetDefaultOrder(pOrder);
    cmesh->AutoBuild();
    std::cout << "\nComputational mesh executed\n";
    cmesh->Print(std::cout);

    Analysis analysis(cmesh);
    analysis.RunSimulation();
    VecDouble Sol = cmesh->Solution();
    std::cout << "\nAnalysis Executed\n";

    PostProcessTemplate<Poisson> process;
    setFuncPostProcess(gmesh->Dimension(),solType,&process);
    analysis.PostProcessError(std::cout,process);
    std::cout << "\nPostProcessError Executed\n";

    process.AppendVariable("Sol");
    process.AppendVariable("SolExact");
    //process.AppendVariable("Force");

    VTKGeoMesh::PrintGMeshVTK(gmesh,"GMESH.vtk");
    VTKGeoMesh::PrintCMeshVTK(cmesh,2,"CMESH.vtk");
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
    int numElem = 6*nx*nx*nx +12*nx*nx;
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
            T[2] = comum;
            T[1] = comum+ip+kp;
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

            T[2] = comum;
            T[0] = comum+jp;
            T[1] = comum + jp + kp;
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

            T[1] = comum + ip + jp;
            T[2] = comum;
            T[0] = comum + ip;
            T[3] = comum + ip + jp + kp;
            GeoElementTemplate<GeomTetrahedron> *tet5 = new GeoElementTemplate<GeomTetrahedron>(T,1,gmesh,index);
            gmesh->SetElement(index,tet5);
            index++;

            T[1] = comum;
            T[0] = comum + jp;
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
            int Lcomum = i + (nx+1)*j+(nx+1);
            int ip = 1, jp = nx+1;

            LBC[0] = Lcomum;
            LBC[1] = Lcomum - jp;
            LBC[2] = Lcomum + ip;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + ip-jp;

            UBC[0] = Ucomum;
            UBC[1] = Ucomum+jp;
            UBC[2] = Ucomum-ip;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index++;
        }
    //XY TOP
    for(int j = 0; j < nx; j++) for(int i = 0; i < nx ; i++){
            int Lcomum =(nx+1)+ i + (nx+1)*j+nx*(nx+1)*(nx+1);
            int ip = 1, jp = nx+1;

            LBC[0] = Lcomum;
            LBC[1] = Lcomum - jp;
            LBC[2] = Lcomum + ip;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + ip-jp;

            UBC[0] = Ucomum;
            UBC[1] = Ucomum+jp;
            UBC[2] = Ucomum-ip;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index++;
        }
    //ZX FRONT
    for(int k = 0; k < nx; k++) for(int i = 0; i < nx ; i++){
            int Lcomum = (nx+1)*(nx+1)+i + (nx+1)*(nx+1)*k;
            int ip = 1, kp = (nx+1)*(nx+1);

            LBC[0] = Lcomum;
            LBC[1] = Lcomum + ip;
            LBC[2] = Lcomum-kp;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + ip-kp;

            UBC[0] = Ucomum;
            UBC[1] = Ucomum+kp;
            UBC[2] = Ucomum-ip;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index++;
        }
    //ZX BACK
    for(int k = 0; k < nx; k++) for(int i = 0; i < nx ; i++){
            int Lcomum = (nx+1)*(nx+1)+ i + (nx+1)*(nx+1)*k+nx*(nx+1);
            int ip = 1, kp = (nx+1)*(nx+1);

            LBC[0] = Lcomum;
            LBC[1] = Lcomum + ip;
            LBC[2] = Lcomum-kp;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + ip-kp;

            UBC[0] = Ucomum;
            UBC[1] = Ucomum+kp;
            UBC[2] = Ucomum-ip;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index++;
        }
    //YZ LEFT
    for(int k = 0; k < nx; k++) for(int j = 0; j < nx ; j++){
            int Lcomum = (nx+1)*(nx+1)+ j*(nx+1)+k*(nx+1)*(nx+1);
            int jp = (nx+1), kp = (nx+1)*(nx+1);

            LBC[0] = Lcomum;
            LBC[1] = Lcomum -kp;
            LBC[2] = Lcomum+jp;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + jp-kp;

            UBC[0] = Ucomum;
            UBC[1] = Ucomum+kp;
            UBC[2] = Ucomum-jp;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index++;
        }
    //YZ RIGHT
    for(int k = 0; k < nx; k++) for(int j = 0; j < nx ; j++){
            int Lcomum = (nx+1)*(nx+1) + j*(nx+1)+k*(nx+1)*(nx+1) +nx*1;
            int jp = (nx+1), kp = (nx+1)*(nx+1);

            LBC[0] = Lcomum;
            LBC[1] = Lcomum -kp;
            LBC[2] = Lcomum+jp;
            GeoElementTemplate<GeomTriangle> *lBC = new GeoElementTemplate<GeomTriangle>(LBC,-1,gmesh,index);
            gmesh->SetElement(index,lBC);
            index++;

            int Ucomum = Lcomum + jp - kp;

            UBC[0] = Ucomum;
            UBC[1] = Ucomum+kp;
            UBC[2] = Ucomum-jp;
            GeoElementTemplate<GeomTriangle> *uBC = new GeoElementTemplate<GeomTriangle>(UBC,-1,gmesh,index);
            gmesh->SetElement(index,uBC);
            index++;
        }



    //XY BOTTOM
    /*VecInt LBC(3,0),UBC(3,0);
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
            index++;
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
            index++;
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
            index++;
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
            index++;
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
            index++;
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
            index++;
        }*/
    gmesh->BuildConnectivity();
    return(gmesh);
}

void setFuncPoisson(const int dim, const int fType, Poisson *pot) {
    if (dim == 2) {
        switch (fType) {
            case 0:
                pot->SetForceFunction(force2D);
                pot->SetExactSolution(laplace2D);
                break;
            case 1:
                pot->SetForceFunction(force2D1P);
                pot->SetExactSolution(laplace2D1P);
                break;
            case 2:
                pot->SetForceFunction(force2D2P);
                pot->SetExactSolution(laplace2D2P);
                break;
            default:
                std::cout << "Undefined solution type";
        }
    } else {
        switch (fType) {
            case 0:
                pot->SetForceFunction(force3D);
                pot->SetExactSolution(laplace3D);
                break;
            case 1:
                pot->SetForceFunction(force3D1P);
                pot->SetExactSolution(laplace3D1P);
                break;
            case 2:
                pot->SetForceFunction(force3D2P);
                pot->SetExactSolution(laplace3D2P);
                break;
            default:
                std::cout << "Undefined solution type";
        }
    }
}

void setFuncL2(const int dim, const int fType, L2Projection *l2) {
    if (dim == 2) {
        switch (fType) {
            case 0:
                l2->SetForceFunction(force2D);
                l2->SetExactSolution(laplace2D);
                break;
            case 1:
                l2->SetForceFunction(force2D1P);
                l2->SetExactSolution(laplace2D1P);
                break;
            case 2:
                l2->SetForceFunction(force2D2P);
                l2->SetExactSolution(laplace2D2P);
                break;
            default:
                std::cout << "Undefined solution type";
        }
    } else {
        switch (fType) {
            case 0:
                l2->SetForceFunction(force3D);
                l2->SetExactSolution(laplace3D);
                break;
            case 1:
                l2->SetForceFunction(force3D1P);
                l2->SetExactSolution(laplace3D1P);
                break;
            case 2:
                l2->SetForceFunction(force3D2P);
                l2->SetExactSolution(laplace3D2P);
                break;
            default:
                std::cout << "Undefined solution type";
        }
    }
}

void setFuncPostProcess(const int dim, const int fType, PostProcessTemplate<Poisson> *pp) {
    if (dim == 2) {
        switch (fType) {
            case 0:
                pp->SetExact(laplace2D);
                break;
            case 1:
                pp->SetExact(laplace2D1P);
                break;
            case 2:
                pp->SetExact(laplace2D2P);
                break;
            default:
                std::cout << "Undefined solution type";
        }
    } else {
        switch (fType) {
            case 0:
                pp->SetExact(laplace3D);
                break;
            case 1:
                pp->SetExact(laplace3D1P);
                break;
            case 2:
                pp->SetExact(laplace3D2P);
                break;
            default:
                std::cout << "Undefined solution type";
        }
    }
}