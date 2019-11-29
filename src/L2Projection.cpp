/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "L2Projection.h"
#include "tpanic.h"
#include "MathStatement.h"
#include <string.h>

L2Projection::L2Projection() {
    BCVal1.Resize(0,0);
    BCVal2.Resize(0,0);
    projection.Resize(0,0);
    SetMatID(-1);
}

L2Projection::L2Projection(int bctype, int materialid, Matrix &proj, Matrix Val1, Matrix Val2) {
    BCType = bctype;
    SetMatID(materialid);
    projection = proj;
    BCVal1 = Val1;
    BCVal2 = Val2;
}

L2Projection::L2Projection(const L2Projection &copy) {
    BCType = copy.GetBCType();
    projection = copy.GetProjectionMatrix();
    BCVal1 = copy.Val1();
    BCVal2 = copy.Val2();
    SetMatID(-1);

    forceFunction = copy.GetForceFunction();
}

L2Projection &L2Projection::operator=(const L2Projection &copy) {
    BCType = copy.GetBCType();
    projection = copy.GetProjectionMatrix();
    BCVal1 = copy.Val1();
    BCVal2 = copy.Val2();
    SetMatID(-1);

    forceFunction = copy.GetForceFunction();
}

L2Projection *L2Projection::Clone() const {
}

L2Projection::~L2Projection() {
}

Matrix L2Projection::GetProjectionMatrix() const {
    return projection;
}

void L2Projection::SetProjectionMatrix(const Matrix &proj) {
    projection = proj;
}

void L2Projection::Contribute(IntPointData &integrationpoint, double weight, Matrix &EK, Matrix &EF) const {
    //Copied from Poisson
    VecDouble ksi = integrationpoint.ksi;
    VecDouble X = integrationpoint.x;
    int size = integrationpoint.phi.size(); // Implementation assumes single state variable

    int type = GetBCType(); // 0 : Dirichlet ; 1 : Neumann ;
    double alfa;
    VecDouble res(0); Matrix dRes(0,0);
    SolutionExact(X,res,dRes);

    switch(type){
        case 0: // Dirichlet;
            for(int i =0; i < size ; i++) {
                alfa = res[0]; //std::cout << "\n\n\nVal1:" << alfa << "\n\n\n";    //Val1 i
                EF(i,0) += weight * gBigNumber * alfa * integrationpoint.phi[i];
                for(int j =0 ; j< size; j++){
                    EK(i,j) += weight*gBigNumber*integrationpoint.phi[i]*integrationpoint.phi[j];
                }
            }
            break;
        case 1: //Neumann;
            for(int i =0; i< size;i++){
                alfa = Val1()(0,0);
                EF(i,0) +=weight*alfa*integrationpoint.phi[i];
            }
            break;
        default: std::cout << "Error: Unknown boundary condition"; DebugStop();
    }

    // From Poisson
    VecDouble f(NState(),0);
    for(int i =0; i< size;i++){
        forceFunction(X,f);
        for(int k = 0; k <NState();k++) EF(i,0) += weight*(f[k]*integrationpoint.phi[i]);    //Loop is necessary to taken into account multiple state variables.
        for(int j =0; j < size; j++){
            EK(i,j) += weight*(integrationpoint.dphidx(i,0)*integrationpoint.dphidx(j,0)+integrationpoint.dphidx(i,1)*integrationpoint.dphidx(j,1));
        }
    }
}

int L2Projection::NEvalErrors() const {
    return 3;
}

void L2Projection::ContributeError(IntPointData &integrationpoint, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
    VecDouble uh; PostProcessSolution(integrationpoint,ESol,uh);
    Matrix duhDx = integrationpoint.dsoldx;// PostProcessSolution(data,EDSol,duhDx);
    errors.resize(3);

    for(int i=0; i <errors.size();i++) errors[i] = 0;

    //Energy norm calculation
    double e;
    errors[2] = 0;
    for (int i =0; i < uh.size(); i++){
        e = uh[i] - u_exact[i];
        errors[0] += e*e;
        if(abs(e) > errors[2]) errors[2] = abs(e);
    }
    errors[1] = errors[0];
    double de;
    for(int i = 0 ; i< duhDx.Rows();i++) for(int j = 0; j <duhDx.Cols(); j++) {
            de = duhDx(i,j) - du_exact(i,j);
            errors[0] += de*de;
        }
}

int L2Projection::VariableIndex(const PostProcVar var) const {
    int nEnum = 3;
    for(int i =0; i< nEnum; i++) if(var == PostProcVar(i)) return i;
}

L2Projection::PostProcVar L2Projection::VariableIndex(const std::string & name) {
    //  enum PostProcVar {ENone, ESol, EDSol};
    if(!strcmp("None",name.c_str())) return ENone;
    if(!strcmp("Sol",name.c_str())) return ESol;
    if(!strcmp("EDSol",name.c_str())) return EDSol;
}

int L2Projection::NSolutionVariables(const PostProcVar var) {
    switch (var) {
        case ESol:
            return NState();
        case EDSol:
            return Dimension();
        case ENone:
            return 0;
    }
}

void L2Projection::PostProcessSolution(const IntPointData &integrationpoint, const int var, VecDouble &sol) const {
    //Modified from Poisson
    VecDouble u_h = integrationpoint.solution;
    Matrix du_hdx = integrationpoint.dsoldx;
    int cols = du_hdx.Cols(), rows = du_hdx.Rows();
    int tot_size = cols*rows;

    switch(var){
        case ESol: //ux, uy and uz
            sol.resize(NState());
            for(int i =0; i< NState();i++) sol[i] = u_h[i]; break;
        case EDSol: //gradU
            sol.resize(tot_size);
            for (int i = 0 ;i< rows ; i++) for(int j = 0 ; j< cols ; j++) sol[j+cols*i] = du_hdx(i,j); break;
        case ENone: break;
        default: std::cout << "Uncovered PostProcVar Value";
    };
}
