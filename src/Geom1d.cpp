/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Geom1d.h"
#include "tpanic.h"

Geom1d::Geom1d(): fNodeIndices({0,1}),fNeighbours() {
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d &copy) : fNodeIndices(copy.fNodeIndices) {
    for(int i = 0; i <nSides; i++) fNeighbours[i] = copy.fNeighbours[i];
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
    for(int i = 0; i <nSides; i++) fNeighbours[i] = copy.fNeighbours[i];
    return *this;
}

void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi) {
    phi[0] = 0.5*(1-xi[0]); phi[1] = 0.5*(1+xi[0]);
    dphi(0,0) = -0.5; dphi(1,0) = 0.5;
}

void Geom1d::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    //x[0] = 0.5*(xi[0]*(NodeCo[1]-NodeCo[0])+(NodeCo[1]+NodeCo[0]));  Simplified relationship.
    VecDouble fphi(2,0); Matrix fdphi(2,1,0);
    Shape(xi,fphi,fdphi); //std::cout << "\nPrint Node Coordinates\n";NodeCo.Print(std::cout); std::cout<< "\nPrint Finished\n";
    int dim = NodeCo.Cols();
    for(int ind =0; ind < dim; ind++ ) x[ind] = NodeCo(0,ind)*fphi[0]+ NodeCo(1,ind)*fphi[1];
}

void Geom1d::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    VecDouble fphi(2,0); Matrix fdphi(2,1,0);
    Shape(xi,fphi,fdphi);
    int dim = NodeCo.Cols();
    gradx.Resize(1,dim);
    for( int spacialInd = 0; spacialInd <dim ; spacialInd++) {
        gradx(0, spacialInd) = NodeCo(0, spacialInd) * fdphi(0, 0) + NodeCo(1, spacialInd) * fdphi(1, 0);
    }
}

void Geom1d::SetNodes(const VecInt &nodes) {
    fNodeIndices.resize(NumNodes());
    for(auto i:fNodeIndices) fNodeIndices[i] = nodes[i];
}

void Geom1d::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
    //nodes.resize(fNodeIndices.size());
    //for(auto i:fNodeIndices) nodes[i] = fNodeIndices[i];
}

int Geom1d::NodeIndex(int node) {
    return(fNodeIndices[node]);
}

int Geom1d::NumNodes() {
    return (2);
}

GeoElementSide Geom1d::Neighbour(int side) {
    return(fNeighbours[side]);
}

void Geom1d::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
