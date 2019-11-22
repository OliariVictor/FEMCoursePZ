/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "CompElementTemplate.h"
#include "CompElement.h"
#include "CompMesh.h"
#include "DOF.h"
#include "MathStatement.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "GeoNode.h"
#include "GeoMesh.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "ShapeTetrahedron.h"
#include "tpanic.h"

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(): CompElement() {
    int shapeNumber = NShapeFunctions();
    dofindexes.resize(shapeNumber);
}

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, CompMesh *cmesh, GeoElement *geo): CompElement(ind,cmesh,geo) {
    int shapeNumber = NShapeFunctions();
    dofindexes.resize(shapeNumber);
}

template<class Shape >
CompElementTemplate<Shape>::CompElementTemplate(const CompElementTemplate & copy): CompElement(copy.GetIndex(),copy.GetCompMesh(),copy.GetGeoElement()) {
    int shapeNumber = NShapeFunctions();
    dofindexes.resize(shapeNumber);
}

template<class Shape >
CompElementTemplate<Shape> &CompElementTemplate<Shape>::operator=(const CompElementTemplate & copy) {
    SetGeoElement(copy.GetGeoElement());
    SetCompMesh(copy.GetCompMesh());
    SetIndex(copy.GetIndex());
    SetIntRule(copy.GetIntRule());
    SetStatement(copy.GetStatement());
    dofindexes = copy.dofindexes;
    return(*this);
}

template<class Shape >
CompElementTemplate<Shape>::~CompElementTemplate() {
}

template<class Shape >
CompElement * CompElementTemplate<Shape>::Clone() const {
    CompElement *result = new CompElementTemplate<Shape>(*this);
    return (result);
}

template<class Shape>
void CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix & dphi) const {
    int N = NDOF();
    VecInt Order(N);
    for(int i = 0; i< N; i++) { Order[i] = GetCompMesh()->GetDOF(dofindexes[i]).GetOrder();}
    Shape::Shape(intpoint,Order,phi,dphi);
}

template<class Shape>
void CompElementTemplate<Shape>::GetMultiplyingCoeficients(VecDouble & coefs) const {
    VecDouble T;
    T = GetCompMesh()->Solution();

    int length = dofindexes.size();
    coefs.resize(length);

    for(int i; i< length;i++) coefs[i] = T[dofindexes[i]];
}

template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions() const{
    //int shapeOrder = GetCompMesh()->GetDefaultOrder();
    //if(shapeOrder != 1 && shapeOrder != 2) {std::cout << "CompElementTemplate<Shape>::NShapeFunctions(): Non-defined shape order";DebugStop();}
    int total = 0,N,cap = 0;
    N = NDOF();
    for(int i = 0; i< N; i++) {total += NShapeFunctions(i);}
    return(total);
}

template<class Shape>
void CompElementTemplate<Shape>::SetNDOF(int64_t ndof) {
    dofindexes.resize(ndof);
}

template<class Shape>
void CompElementTemplate<Shape>::SetDOFIndex(int i, int64_t dofindex) {
    dofindexes[i] = dofindex;
}

template<class Shape>
int64_t CompElementTemplate<Shape>::GetDOFIndex(int i) {
    return(dofindexes[i]);
}

template<class Shape>
int CompElementTemplate<Shape>::NDOF() const {
    return dofindexes.size();
}

template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex) const {
   int globalInd =  dofindexes[doflocindex];
   GetCompMesh()->GetDOF(globalInd).GetNShape();
}

template<class Shape>
int CompElementTemplate<Shape>::ComputeNShapeFunctions(int doflocindex, int order) {
    return(Shape::NShapeFunctions(doflocindex,order));
}

template<class Shape>
void CompElementTemplate<Shape>::Print(std::ostream &out) {
    //Copied from PZ`s CompEl Print

    out << "\nOutput for a computable element index: " << GetIndex();
    //out << "\nfReferenceIndex " << fReferenceIndex;
    if(this->GetGeoElement())
    {
        out << "\nGeometry: Center Coordinates: ";
        VecDouble coord(3,0);
        int dim = Shape::Dimension;
        for(int i =0; i < Shape::nCorners;i++) for(int j =0; j< dim;j++) coord[j] += GetGeoElement()->GetMesh()->Node(GetGeoElement()->NodeIndex(i)).Coord(j);
        for(int i =0; i <dim; i++) std::cout << "[" << i << "] ,";
    }
    if(this->GetStatement())
    {
        out << "\nMaterial id " << GetStatement()->GetMatID() << "\n";
    }
    else {
        out << "\nNo material\n";
    }

    out << "Number of connects = " << NDOF();
    out<< "\nConnect indexes : ";
    int nod;
    for(nod=0; nod< NDOF(); nod++)
    {
        out << GetDOFIndex(nod) <<  ' ' ;
    }
    GetIntRule()->Print(std::cout);
    out << std::endl;

}



template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;
