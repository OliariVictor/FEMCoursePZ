/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Poisson.h"
#include "tpanic.h"
#include "string.h"

Poisson::Poisson() {
    SetMatID(1);
}

Poisson::Poisson(int materialid, Matrix &perm) {
    SetMatID(materialid);
    permeability = perm;
}

Poisson::Poisson(const Poisson &copy) {
    forceFunction = copy.GetForceFunction();
    permeability = copy.permeability;
    SolutionExact = copy.SolutionExact;
    SetMatID(1);
}

Poisson &Poisson::operator=(const Poisson &copy) {
    forceFunction = copy.GetForceFunction();
    permeability = copy.permeability;
    SolutionExact = copy.SolutionExact;

    SetMatID(1);
    return(*this);
}

Poisson *Poisson::Clone() const {
}

Poisson::~Poisson() {
}

Matrix Poisson::GetPermeability() const {
    return permeability;
}

void Poisson::SetPermeability(const Matrix &perm) {
    permeability = perm;
}

int Poisson::NEvalErrors() const {
    return Dimension();
}

int Poisson::VariableIndex(const PostProcVar var) const {
    int nEnum = 7;
    for(int i =0; i< nEnum; i++) if(var == PostProcVar(i)) return i;
}

Poisson::PostProcVar Poisson::VariableIndex(const std::string &name) {
    //  enum PostProcVar {ENone, ESol, EDSol, EFlux, EForce, ESolExact, EDSolExact};
    if(!strcmp("None",name.c_str())) return ENone;
    if(!strcmp("Sol",name.c_str())) return ESol;
    if(!strcmp("EDSol",name.c_str())) return EDSol;
    if(!strcmp("Flux",name.c_str())) return EFlux;
    if(!strcmp("Force",name.c_str())) return EForce;
    if(!strcmp("SolExact",name.c_str())) return ESolExact;
    if(!strcmp("DSolExact",name.c_str())) return EDSolExact;
}

int Poisson::NSolutionVariables(const PostProcVar var) {
    switch(var) {
        case ESol: return NState();
        case EDSol: return Dimension();
        case EFlux: return Dimension();
        case EForce: return NState();
        case ESolExact: return NState();
        case EDSolExact: return Dimension();
        case ENone: return 0;
    }
}

void Poisson::ContributeError(IntPointData &data, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
    VecDouble uh; PostProcessSolution(data,ESol,uh);
    Matrix duhDx = data.dsoldx;// PostProcessSolution(data,EDSol,duhDx);
    errors.resize(3);

    for(int i=0; i <errors.size();i++) errors[i] = 0;

    //Energy norm calculation
    double e;
    errors[2] = 0;
    for (int i =0; i < uh.size(); i++){
        e = uh[i] - u_exact[i];
        errors[0] += e*e;
        if (e < 0) {if (-1*e>errors[2]) errors[2] = -1*e;}
        else if(e > errors[2]) errors[2] = e;
    }
    errors[1] = errors[0]; double sum;
    double de; //std::cout << "\n\nDuDx\n"; duhDx.Print(std::cout); std::cout << "\n\nDu_Exact\n"; du_exact.Print(std::cout);
    for(int i = 0 ; i< duhDx.Rows();i++) for(int j = 0; j <duhDx.Cols(); j++) {
        de = duhDx(i,j) - du_exact(i,j);
        errors[0] += de*de;
    }
}

void Poisson::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const {
    //Poisson's problem requires no K coefficent and no Permeability Matrix
    VecDouble ksi = data.ksi;
    VecDouble X = data.x;
    int size = data.phi.size();
    //std::cout << "\ndphidx:\n";
    //data.dphidx.Print(std::cout);
    VecDouble f(NState(),0);
    for(int i =0; i< size;i++){
        forceFunction(X,f);
        for(int k = 0; k <NState();k++) {EF(i,0) += weight*(f[k]*data.phi[i]);}   //std::cout<<"\nf[k] = "<<f[k]<<"\nphi[i] =" << data.phi[i]<< "\n\nDetJac: "<< data.detjac<<  std::endl;}//Loop is necessary to taken into account multiple state variables.
        for(int j =0; j < size; j++){
            EK(i,j) += weight*(data.dphidx(i,0)*data.dphidx(j,0)+data.dphidx(i,1)*data.dphidx(j,1));
        }
    }
}
//Contribute is only implemented for ESol, nonetheless, all other PostProcVar have no application and were filled only for didactic purposes
void Poisson::PostProcessSolution(const IntPointData &data, const int var, VecDouble &Solout) const {
    PostProcVar v = PostProcVar(var);

    VecDouble u_h = data.solution;
    Matrix du_hdx = data.dsoldx;
    int cols = du_hdx.Cols(), rows = du_hdx.Rows();
    int tot_size = cols*rows;
    //Declarations for EFlux
    Matrix Flux(permeability.Rows(),cols,0);
    VecDouble f(NState(),0);
    //Declarations for ESolExact and EDSolExact.
    VecDouble Sol(NState(),0.);
    Matrix dSol(NState(),1,0.);

    switch(var){
        case ESol: //ux, uy and uz
            Solout.resize(NState());
            for(int i =0; i< NState();i++) Solout[i] = u_h[i]; break;
        case EDSol: //gradU
            Solout.resize(tot_size);
            for (int i = 0 ;i< rows ; i++) for(int j = 0 ; j< cols ; j++) Solout[j+cols*i] = du_hdx(i,j); break;
        case EFlux: //Flux = Perm x Grad u
            Solout.resize(permeability.Rows()*cols);

            //Computing the Flux Matrix
            for (int i = 0 ; i< permeability.Rows(); i++) for(int j =0; j<cols ; j++) for(int k =0 ; k<permeability.Cols();k++) Flux(i,j) += permeability.GetVal(i,k)*du_hdx(k,j);
            //Reshaping Flux Matrix into SoloutVector
            for (int i = 0 ;i< rows ; i++) for(int j = 0 ; j< cols ; j++) Solout[j+cols*i] = Flux(i,j); break;
        case EForce: //f : Solout returns forceVector
            Solout.resize(NState());
             forceFunction(data.x,f);

            for(int i =0; i< NState(); i++) Solout[i] = f[i]; break;
        case ESolExact: //u_exact
        case EDSolExact: //du_exact
            Solout.resize(NState());

            if(SolutionExact) SolutionExact(data.x,Sol,dSol);

            if(var == ESolExact) {for(int i =0; i<NState(); i++) Solout[i] = Sol[i]; break;}
            else for(int i =0; i<dSol.Rows(); i++) for(int j =0; j<dSol.Cols(); j++)  Solout[j+dSol.Cols()*i] = dSol(i,j); break;
        case ENone: break;
        default: std::cout << "Uncovered PostProcVar Value";
    }

}

double Poisson::Inner(Matrix &S, Matrix & T) const {
    double sum = 0;
    if(S.Cols() != T.Cols() || S.Rows() != T.Rows()) { std::cout << "Poisson::Inner: Matrix must have the same size"; DebugStop(); }
    for(int i =0 ; i< S.Cols(); i++) for(int j =0; j< S.Rows(); j++) sum += S(i,j)*T(i,j);
}
