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
    dphi(0,0) = -0.5; dphi(0,1) = 0.5;
}

void Geom1d::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    //x[0] = 0.5*(xi[0]*(NodeCo[1]-NodeCo[0])+(NodeCo[1]+NodeCo[0]));  Simplified relationship.
    VecDouble fphi(2,0); Matrix fdphi(1,2,0);
    Shape(xi,fphi,fdphi);
    x[0] = NodeCo(0,0)*fphi[0]+ NodeCo(0,1)*fphi[1];
}

void Geom1d::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    VecDouble fphi(2,0); Matrix fdphi(1,2,0);
    Shape(xi,fphi,fdphi);
    gradx(0,0) = NodeCo(0,0)*fdphi(0,0)+NodeCo(0,1)*fdphi(0,1);
    //gradx(1,1) = 0.5*(NodeCo[1]-NodeCo[0]);
}

void Geom1d::SetNodes(const VecInt &nodes) {
    for(auto i:fNodeIndices) fNodeIndices[i] = nodes[i];
}

void Geom1d::GetNodes(VecInt &nodes) {
    for(auto i:fNodeIndices) nodes[i] = fNodeIndices[i];
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
