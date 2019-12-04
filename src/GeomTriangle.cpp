/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTriangle.h"
#include "tpanic.h"

GeomTriangle::GeomTriangle(): fNodeIndices({0,1,2}),fNeighbours()  {
}

GeomTriangle::~GeomTriangle() {
}

GeomTriangle::GeomTriangle(const GeomTriangle &copy): fNodeIndices(copy.fNodeIndices) {
    for(int i = 0; i <nSides; i++) fNeighbours[i] = copy.fNeighbours[i];
}

GeomTriangle& GeomTriangle::operator=(const GeomTriangle& copy) {
    fNodeIndices = copy.fNodeIndices;
    for(int i = 0; i <nSides; i++) fNeighbours[i] = copy.fNeighbours[i];
    return *this;
}

void GeomTriangle::Shape(const VecDouble& xi, VecDouble& phi, Matrix& dphi) {
    phi[0] = 1 - xi[0] - xi[1]; dphi(0,0) = -1; dphi(0,1) = -1;
    phi[1] = xi[0]; dphi(1,0) = 1; dphi(1,1) = 0;
    phi[2] = xi[1]; dphi(2,0) = 0; dphi(2,1) = 1;
}

void GeomTriangle::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    VecDouble fphi(3,0); Matrix fdphi(3,2,0);
    Shape(xi,fphi,fdphi);

    int dim = NodeCo.Cols();
    x.resize(dim); for(int i =0 ; i < dim ; i++) x[i] = 0; //std::cout << "NodeCoord\n"; NodeCo.Print();

    for(int i = 0 ; i<3 ; i++){
        for(int j =0; j < dim; j++)
            x[j] += NodeCo(i,j)*fphi[i];
    }
}

void GeomTriangle::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    VecDouble fphi(3,0); Matrix fdphi(3,2,0);
    Shape(xi,fphi,fdphi);

    int dim = NodeCo.Cols(); //std::cout << "\n\nTriangle Node coordinate:\n";NodeCo.Print(); std::cout << "\n\nfdphi\n";fdphi.Print();
    gradx.Resize(dim,2);
    for(int i =0 ; i < dim ; i++) for(int j =0; j<2 ; j++) gradx(i,j) = 0;

    for(int i = 0; i < 2; i++) for(int j =0 ; j < 3;j++) {
        for(int k = 0; k < dim ; k++)
            gradx(k,i) += NodeCo(j,k)*fdphi(j,i);
    }
}

void GeomTriangle::SetNodes(const VecInt &nodes) {
    fNodeIndices = nodes;
}

void GeomTriangle::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
}

int GeomTriangle::NodeIndex(int node) {
    return(fNodeIndices[node]);
}

int GeomTriangle::NumNodes() {
    return(3);
}

GeoElementSide GeomTriangle::Neighbour(int side) {
    return(fNeighbours[side]);
}

void GeomTriangle::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
