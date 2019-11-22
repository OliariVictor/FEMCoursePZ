/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoElementSide.h"
#include "GeoElement.h"
#include "tpanic.h"
#include "GeoElementTemplate.h"
#include "GeoMesh.h"
#include "CompElementTemplate.h"
#include "CompElement.h"
#include "CompMesh.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "ShapeTetrahedron.h"

GeoElement::GeoElement() : MaterialId(-1), Index(-1) { GMesh = 0; Reference = 0;
}

GeoElement::GeoElement(int materialid, GeoMesh *mesh, int index) : Reference() {
    MaterialId = materialid;
    GMesh = mesh;
    Index = index;

    GMesh->SetElement(index,this);
}

GeoElement::GeoElement(const GeoElement &copy) {
    MaterialId = copy.MaterialId;
    GMesh = copy.GMesh;
    Reference = copy.Reference;
    Index = copy.Index;
}

GeoElement::~GeoElement() {
}

CompElement *GeoElement::CreateCompEl(CompMesh *mesh, int64_t index) {
    CompElement *cel;
    switch(Type()) {
        case EOned:
            cel = new CompElementTemplate<Shape1d>(index,mesh,this);
            Reference = cel;
            return cel;
        case ETriangle:
            cel = new CompElementTemplate<ShapeTriangle>(index,mesh,this);
            Reference = cel;
            return cel;
        case ETetraedro:
            cel = new CompElementTemplate<ShapeTetrahedron>(index,mesh,this);
            Reference = cel;
            return cel;
        case EQuadrilateral:
            cel = new CompElementTemplate<ShapeQuad>(index,mesh,this);
            Reference = cel;
            return cel;
        default:
            std::cout<< "CompMesh::AutoBuild(): Invalid Type()\n";DebugStop();
    }
}

void GeoElement::Print(std::ostream &out) {
    std::cout << std::endl <<"Index: "<< Index << "\tMaterialId: "<< MaterialId << "\n";

    for (int i = 0;i < NSides();i++) {
        out << "Neighbours for side   " << i << " : ";
        GeoElementSide neighbour = Neighbour(i);
        GeoElementSide thisside(this,i);
        if (!neighbour.Exists())
        {
            out << "No neighbour\n";
        }
        else {
            while (neighbour != thisside ) {
                out << "Element neighbour: "<< neighbour.Element()->Index << "\t Neighbour`s side: " << neighbour.Side() << " || ";
                neighbour = neighbour.Neighbour();
            }
            out << std::endl;
        }
    }
}
