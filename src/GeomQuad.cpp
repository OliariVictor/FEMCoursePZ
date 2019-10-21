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
    Shape(xi,fphi,fdphi); x[0] = 0;x[1] = 0;
    for(int i = 0 ; i<4 ; i++){
        x[0] += NodeCo(i,0)*fphi[i];
        x[1] += NodeCo(i,1)*fphi[i];
    }
}

void GeomQuad::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    VecDouble fphi(4,0); Matrix fdphi(4,2,0);
    Shape(xi,fphi,fdphi);
    gradx(0,0) = gradx(0,1) = gradx(1,0) =gradx(1,1) = 0;
    for(int i = 0; i < 2; i++) for(int j =0 ; j < 4;j++) {
        gradx(0,i) += NodeCo(j,0)*fdphi(j,i);
        gradx(1,i) += NodeCo(j,1)*fdphi(j,i);
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
