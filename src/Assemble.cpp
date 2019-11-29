/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Assemble.h"
#include "tpanic.h"
#include "CompMesh.h"
#include "CompElement.h"
#include "GeoMesh.h"
#include "GeoElementSide.h"
#include "DOF.h"

Assemble::Assemble() : cmesh() {
}

Assemble::Assemble(CompMesh *mesh) {
    cmesh = mesh;
}

Assemble::Assemble(const Assemble &copy) {
    cmesh = copy.cmesh;
}

Assemble &Assemble::operator=(const Assemble &copy) {
    cmesh = copy.cmesh;
    return(*this);
}

void Assemble::SetMesh(CompMesh *mesh) {
    cmesh = mesh;
}

int64_t Assemble::NEquations() {
    //WRONG IMPLEMENTATION!!!!! MUST TAKE INTO ACCOUNT NEIGHBOURING SIDES
    /*int ncel = cmesh->GetElementVec().size();
    int sumNodes = 0;
    int numNodes = cmesh->GetGeoMesh()->NumNodes();
    int numElem = cmesh->GetGeoMesh()->NumElements();

    //Geometric Solution
    for (int i = 0; i < ncel; i++)
        sumNodes += cmesh->GetElement(i)->NDOF();

    int cOrder = cmesh->GetDefaultOrder();
    if (cOrder != 1 && cOrder != 2) {
        std::cout << "Assemble::NEquations(): Computation Order not implemented";
        DebugStop();
    }

    int nodeEq= 0, repeated = 0;
    //Iterating through each element and computing how many times each node repeat itself
    VecInt Nodes(numNodes, -1);
    for (int elemIndex = 0; elemIndex < numElem; elemIndex++) {
        GeoElement *gel = cmesh->GetGeoMesh()->Element(elemIndex);
        int numCorners = gel->NCornerNodes();
        for (int locCorner = 0; locCorner < numCorners; locCorner++) {
            int globalIndex = gel->NodeIndex(locCorner);
            Nodes[globalIndex] += 1;
        }
    }
    for (int Node : Nodes) repeated += Node;
    nodeEq = sumNodes - repeated;

    Matrix M(0, 2);
    int sideCount = 0, mCountI = 0, mCountF = 0, neiSize = 0;
    bool repElemSide;
    std::vector<GeoElementSide> allNei(0);

    for (int elemIndex = 0; elemIndex < numElem; elemIndex++) {
        GeoElement *gel = cmesh->GetGeoMesh()->Element(elemIndex);
        for (int locSide = gel->NCornerNodes(); locSide < gel->NSides(); locSide++) {
            GeoElementSide gelSide(gel, locSide);
            repElemSide = false;
            for (int mIndex = 0; mIndex < M.Rows(); mIndex++)
                if (M(mIndex, 0) == elemIndex && M(mIndex, 1) == locSide) repElemSide = true;
            if (!repElemSide) {
                sideCount += 1;
                allNei.resize(0);
                gelSide.AllNeighbours(allNei);
                neiSize = allNei.size();
                mCountF += neiSize;
                M.Resize(mCountF, 2);
                for (int neiIndex = mCountI; neiIndex < mCountF; neiIndex++) {
                    M(neiIndex, 0) = allNei[neiIndex].Element()->GetIndex();
                    M(neiIndex, 1) = allNei[neiIndex].Side();
                }
                mCountI = mCountF;
            }
        }
    }
    int total1 =0;
    if (cOrder == 1) total1 = nodeEq;
    else total1 = sideCount + nodeEq;
    */
    //Computational solution
    int64_t nDof = cmesh->GetNumberDOF();
    int total2 = 0;
    for (int64_t dofIndex = 0; dofIndex < nDof; dofIndex++) {
        DOF dof = cmesh->GetDOF(dofIndex);
        total2 += dof.GetNState()*dof.GetNShape();
    }
    ///if(total1 != total2) {std::cout<<"Assemble::NEquations(): Warning: Geometric solution is different from Computational Solution\n"; DebugStop;}
    return(total2);
}

void Assemble::OptimizeBandwidth() { // This method is not necessary to solve the given problem. Hence It won`t be filled up at this moment.
}

void Assemble::Compute(Matrix &globmat, Matrix &rhs) {
    int size = NEquations();
    globmat.Resize(size,size); globmat.Zero();
    rhs.Resize(size,1); rhs.Zero();
    Matrix force(size,1,0);

    int cOrder = cmesh->GetDefaultOrder();
    if (cOrder != 1 && cOrder !=2) {std::cout << "Computation Order not implemented"; DebugStop();}

    int ncel = cmesh->GetElementVec().size();
    for(int64_t celIndex = 0; celIndex < ncel ; celIndex ++){
        CompElement *cel = cmesh->GetElement(celIndex);
        Matrix EK(0,0),EF(0,0);
        cel->CalcStiff(EK,EF); std::cout << "\nElem :" << celIndex << "\nEK: \n"; EK.Print(std::cout);std::cout << "\n\nEF: \n"; EF.Print(std::cout); std::cout <<std::endl;


        int ndof = cel->NDOF();
        VecInt dofEquations(EK.Rows());
        for(int dofIndex =0; dofIndex < ndof ; dofIndex++){
            DOF dof = cmesh->GetDOF(cel->GetDOFIndex(dofIndex));
            int size = dof.GetNState()*dof.GetNShape();
            for(int index = 0; index < size; index++)
                dofEquations[dofIndex+index] = dof.GetFirstEquation()+index;
        }

        int i,j;
        for(int row =0; row <EK.Rows(); row++) {
            i = dofEquations[row];
            rhs(i,0) += EF(row,0);
            for(int col =0; col < EK.Cols() ; col++){
                j = dofEquations[col];
                globmat(i,j) += EK(row,col);
            }
        }
    }
}
