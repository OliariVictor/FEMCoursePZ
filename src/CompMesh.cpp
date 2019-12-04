/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoElement.h"
#include "CompElementTemplate.h"
#include "CompMesh.h"
#include "GeoMesh.h"
#include "tpanic.h"
#include "MathStatement.h"
#include "ShapeTetrahedron.h"
#include "ShapeTriangle.h"
#include "ShapeQuad.h"
#include "Shape1d.h"
#include "DOF.h"

CompMesh::CompMesh() {
    geomesh = 0;

}

CompMesh::CompMesh(const CompMesh &copy): geomesh(), compelements(0), dofs(0), mathstatements(0), solution(0,0) {
    DefaultOrder = 1;
}

CompMesh::CompMesh(GeoMesh *gmesh): geomesh(gmesh), compelements(0), dofs(0), mathstatements(0), solution(0,0) {
    DefaultOrder = 1;
}

CompMesh::~CompMesh() {
}

GeoMesh *CompMesh::GetGeoMesh() const {
    return(geomesh);
}

void CompMesh::SetGeoMesh(GeoMesh *gmesh) {
    geomesh = gmesh;
}

void CompMesh::SetNumberElement(int64_t nelem) {
    compelements.resize(nelem);
}

void CompMesh::SetNumberDOF(int64_t ndof) {
    dofs.resize(ndof);
}

void CompMesh::SetNumberMath(int nmath) {
    mathstatements.resize(nmath);
}

void CompMesh::SetElement(int64_t elindex, CompElement *cel) {
    compelements[elindex] = cel;
}

void CompMesh::SetDOF(int64_t index, const DOF &dof) {
    dofs[index] = dof;
}

void CompMesh::SetMathStatement(int index, MathStatement *math) {
    mathstatements[index] = math;
}

DOF &CompMesh::GetDOF(int64_t dofindex) {
    return(dofs[dofindex]);
}

CompElement *CompMesh::GetElement(int64_t elindex) const {
    return(compelements[elindex]);
}

MathStatement *CompMesh::GetMath(int matindex) const {
    return(mathstatements[matindex]);
}

std::vector<CompElement *> CompMesh::GetElementVec() const {
    return(compelements);
}

std::vector<DOF> CompMesh::GetDOFVec() const {
    return(dofs);
}

std::vector<MathStatement *> CompMesh::GetMathVec() const {
    return(mathstatements);
}

void CompMesh::SetElementVec(const std::vector<CompElement *> &vec) {
    compelements = vec;
}

void CompMesh::SetDOFVec(const std::vector<DOF> &dofvec) {
    dofs = dofvec;
}

void CompMesh::SetMathVec(const std::vector<MathStatement *> &mathvec) {
    mathstatements = mathvec;
}

void CompMesh::AutoBuild() {
    int numElem = GetGeoMesh()->NumElements();
    SetNumberElement(numElem);
    for(int i =0; i<numElem; i++){
        GeoElement *geo = GetGeoMesh()->Element(i);
        CompElement *cel = geo->CreateCompEl(this,i);
        SetElement(i,cel);
    }
    Resequence();
}
//Initialize first equation of DOFs;
void CompMesh::Resequence() {
    int64_t nDof = GetNumberDOF();
    int64_t firstEq = 0;
    int dofSize;
    for(int iDof = 0; iDof < nDof; iDof++){
        GetDOF(iDof).SetFirstEquation(firstEq);
        dofSize = GetDOF(iDof).GetNShape()*GetDOF(iDof).GetNState();
        firstEq += dofSize;
    }
}
// Initialize the datastructure FirstEquation of the DOF objects in the order specified by the vector
void CompMesh::Resequence(VecInt &DOFindices) {
    int64_t nDof = GetNumberDOF();
    for (int iDof = 0; iDof < nDof; iDof++) {
        GetDOF(iDof).SetFirstEquation(DOFindices[iDof]);
    }
}

std::vector<double> &CompMesh::Solution() {
    return(solution);
}

void CompMesh::LoadSolution(std::vector<double> &Sol) {
    solution = Sol;
}

void CompMesh::Print(std::ostream & out) {
    out << "\n\t\tCOMPUTABLE GRID INFORMATIONS:\n\n";
    int numConnect;
    int numDOFtotal = GetNumberDOF();
    int numDOFp2 = 0;
    for (int64_t dofInd =0; dofInd < numDOFtotal; dofInd++){
        DOF dof = this->GetDOF(dofInd);
        numDOFp2 += dof.GetNShape()*dof.GetNState();
    }
    if(DefaultOrder ==1){
        numConnect = GetGeoMesh()->NumNodes(); out << "number of connects            = " << numConnect  << std::endl;
    }
    else {
        numConnect = numDOFp2; out << "number of connects            = " << numConnect << std::endl;
    }

    out << "number of elements            = " << this->GetElementVec().size() << std::endl;
    out << "number of materials           = " << this->GetMathVec().size() << std::endl;
    out << "dimension of the mesh         = " << this->GetMathVec()[0]->Dimension() << std::endl;

    out << "\n\t Connect Information:\n\n";
    int64_t i, nelem = GetNumberDOF();
    for (i = 0; i < nelem; i++) {
        out << "\n---------DOF---------\n";
        out << "Index: " << i << std::endl;
        GetDOFVec()[i].Print(*this, out);
    }
    out << "\n\t Computable Element Information:\n\n";
    nelem = GetElementVec().size();
    for (i = 0; i < nelem; i++) {
        if (!GetElement(i)) continue;
        CompElement *el = GetElement(i);
        out << "\n------Element " << i << " -------\n";
        el->Print(out);
    }
    out << "\n\t Material Information:\n\n";
    nelem = GetMathVec().size();
    for (i = 0; i < nelem; i++) {
        MathStatement *mat = GetMathVec()[i];
        if (!mat) DebugStop();
        out << "Element: " << i<< "\n";
        mat->Print(out);
    }
}