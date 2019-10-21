/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTetrahedron.h"
#include "tpanic.h"

GeomTetrahedron::GeomTetrahedron(): fNodeIndices({0,1,2,3,4}),fNeighbours()  {
}

GeomTetrahedron::~GeomTetrahedron() {
}

GeomTetrahedron::GeomTetrahedron(const GeomTetrahedron &copy): fNodeIndices(copy.fNodeIndices) {
    for(int i = 0; i <nSides; i++) fNeighbours[i] = copy.fNeighbours[i];
}

GeomTetrahedron& GeomTetrahedron::operator=(const GeomTetrahedron& copy){
    fNodeIndices = copy.fNodeIndices;
    for(int i = 0; i <nSides; i++) fNeighbours[i] = copy.fNeighbours[i];
    return *this;
}

void GeomTetrahedron::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi) {
    phi[0] = -xi[0] - xi[1] - xi[2] +1; dphi(0,0) = dphi(0,1) = dphi(0,2) = -1;
    phi[1] = xi[0]; dphi(1,0) = 1; dphi(1,1) = dphi(1,2) = 0;
    phi[2] = xi[1]; dphi (2,1) = 1; dphi(2,0) = dphi(2,2) = 0;
    phi[3] = xi[2]; dphi (3,2) = 1; dphi(3,0) = dphi(3,1) = 0;

}

void GeomTetrahedron::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    VecDouble fphi(4,0); Matrix fdphi(4,3,0);
    Shape(xi,fphi,fdphi); x[0] = x[1] = x[2] = 0;
    for(int i = 0 ; i<4 ; i++){
        x[0] += NodeCo(i,0)*fphi[i];
        x[1] += NodeCo(i,1)*fphi[i];
        x[2] += NodeCo(i,2)*fphi[i];
    }
}


void GeomTetrahedron::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    VecDouble fphi(4,0); Matrix fdphi(4,3,0);
    Shape(xi,fphi,fdphi);
    for (int i =0; i<3; i++) for(int j =0; j<3 ; j++) gradx(i,j) = 0;
    for(int i = 0; i < 3; i++) for(int j =0 ; j < 4;j++) {
            gradx(0,i) += NodeCo(j,0)*fdphi(j,i);
            gradx(1,i) += NodeCo(j,1)*fdphi(j,i);
            gradx(2,i) += NodeCo(j,2)*fdphi(j,i);
        }
}


void GeomTetrahedron::SetNodes(const VecInt &nodes) {
    fNodeIndices = nodes;
}

void GeomTetrahedron::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
}

int GeomTetrahedron::NodeIndex(int node) {
    return(fNodeIndices[node]);
}

int GeomTetrahedron::NumNodes() {
    return(5);

}

GeoElementSide GeomTetrahedron::Neighbour(int side) {
    return(fNeighbours[side]);
}

void GeomTetrahedron::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
