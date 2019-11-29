/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "GeomQuad.h"
#include "tpanic.h"

GeomQuad::GeomQuad(): fNodeIndices({0,1,2,3}),fNeighbours() {
}

GeomQuad::~GeomQuad() {
}

GeomQuad::GeomQuad(const GeomQuad &copy) : fNodeIndices(copy.fNodeIndices) {
    for(int i = 0; i <nSides; i++) fNeighbours[i] = copy.fNeighbours[i];
}

GeomQuad& GeomQuad::operator=(const GeomQuad& copy) {
    fNodeIndices = copy.fNodeIndices;
    for(int i = 0; i <nSides; i++) fNeighbours[i] = copy.fNeighbours[i];
    return *this;
}

void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi) {
    phi[0] = 0.25*(1-xi[0])*(1-xi[1]); dphi(0,0) = -.25*(1-xi[1]); dphi(0,1) = -.25*(1-xi[0]);
    phi[1] = 0.25*(1+xi[0])*(1-xi[1]); dphi(1,0) = +.25*(1-xi[1]); dphi(1,1) = -.25*(1+xi[0]);
    phi[2] = 0.25*(1+xi[0])*(1+xi[1]); dphi(2,0) = +.25*(1+xi[1]); dphi(2,1) = +.25*(1+xi[0]);
    phi[3] = 0.25*(1-xi[0])*(1+xi[1]); dphi(3,0) = -.25*(1+xi[1]); dphi(3,1) = +.25*(1-xi[0]);
}

void GeomQuad::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    VecDouble fphi(4,0); Matrix fdphi(4,2,0);
    Shape(xi,fphi,fdphi);
    int dim = NodeCo.Cols();
    x.resize(dim); for(int i =0 ; i < dim ; i++) x[i] = 0;

    for(int i = 0 ; i<4 ; i++){
        for(int j =0; j<dim ; j++)
            x[j] += NodeCo(i,j)*fphi[i];
    }
}

void GeomQuad::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    VecDouble fphi(4,0); Matrix fdphi(4,2,0);
    Shape(xi,fphi,fdphi);

    int dim = NodeCo.Cols();
    gradx.Resize(dim,2);
    for(int i =0 ; i < dim ; i++) for(int j =0; j<2 ; j++) gradx(i,j) = 0;
    //std::cout<<"\n\n\nNodeCo\n\n\n";NodeCo.Print(std::cout);
    //std::cout<<"\n\n\ndPhi\n\n\n"; fdphi.Print(std::cout);
    for(int i = 0; i < 2; i++) for(int j =0 ; j < 4;j++) {
        for(int k = 0; k < dim ; k++)
            gradx(k,i) += NodeCo(j,k)*fdphi(j,i);
    }
}

void GeomQuad::SetNodes(const VecInt &nodes) {
    fNodeIndices = nodes;
}

void GeomQuad::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
}

int GeomQuad::NodeIndex(int node) {
    return(fNodeIndices[node]);
}

int GeomQuad::NumNodes() {
    return(4);
}

GeoElementSide GeomQuad::Neighbour(int side) {
    return(fNeighbours[side]);
}

void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
