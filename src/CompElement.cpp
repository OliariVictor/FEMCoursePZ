/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "CompElement.h"
#include "tpanic.h"
#include "CompMesh.h"
#include "CompElementTemplate.h"
#include "GeoMesh.h"
#include "MathStatement.h"
#include "PostProcess.h"
#include "PostProcessTemplate.h"
#include "IntPointData.h"
#include "GeoElement.h"
#include "DataTypes.h"

using namespace std;

CompElement::CompElement() { //Variables initialized at header
    geoel->SetReference(this);
}

CompElement::CompElement(int64_t ind, CompMesh *cmesh, GeoElement *geo) {
    index = ind;
    compmesh = cmesh;
    geoel = geo;
    geoel->SetReference(this);
}

CompElement::CompElement(const CompElement &copy) {
    index = copy.index;
    compmesh = copy.compmesh;
    geoel = copy.geoel;
    intrule = copy.intrule;
    mat = copy.mat;
    geoel->SetReference(this);
}

CompElement &CompElement::operator=(const CompElement &copy) {
    index = copy.index;
    compmesh = copy.compmesh;
    geoel = copy.geoel;
    intrule = copy.intrule;
    mat = copy.mat;
    geoel->SetReference(this);
}

CompElement::~CompElement() {
}

CompElement *CompElement::Clone() const {
    //CompElement *result = new CompElement(*this);
    //return (*this);
}

MathStatement *CompElement::GetStatement() const {
    return(mat);
}

void CompElement::SetStatement(MathStatement *statement) {
    mat = statement;
}

IntRule *CompElement::GetIntRule() const {
    return(intrule);
}

void CompElement::SetIntRule(IntRule *irule) {
    intrule = irule;
}

void CompElement::SetIndex(int64_t ind) {
    index = ind;
}

GeoElement *CompElement::GetGeoElement() const {
    return(geoel);
}

void CompElement::SetGeoElement(GeoElement *element) {
    geoel = element;
}

CompMesh *CompElement::GetCompMesh() const {
   return(compmesh);
}

void CompElement::SetCompMesh(CompMesh *mesh) {
    compmesh = mesh;
}

void CompElement::InitializeIntPointData(IntPointData &data) const {
    int dim = geoel->GetMesh()->Dimension();int eleDim = geoel->GetReference()->Dimension();
    int shapeNum = this->NShapeFunctions();
    int nState = this->mat->NState();
    data.ksi.resize(dim,-10);
    data.weight = -10;
    data.phi.resize(shapeNum);
    data.dphidksi.Resize(shapeNum,eleDim);
    data.x.resize(dim);
    data.gradx.Resize(eleDim,dim);
    data.axes.Resize(eleDim,eleDim);
    data.detjac = 0;
    data.dphidx.Resize(shapeNum,dim);
    data.solution.resize(nState);
    data.dsoldksi.Resize(dim,nState);   // "dim" is in rows to honor the pattern defined;
    data.dsoldx.Resize(dim,nState);   // "dim" is in rows to honor the pattern defined;
    data.coefs.resize(shapeNum*nState);
}

void CompElement::ComputeRequiredData(IntPointData &data, VecDouble &intpoint) const {
    data.ksi = intpoint;
    //filling weight
    int n = intrule->NPoints();
    int dim = intpoint.size();
    VecDouble coor(dim,0);

    //Weight is always computed elsewhere
    /*double wei = 0; int p = 0;
    for(p=0; p<n; p++){
        intrule->Point(p,coor,wei);
        double epsilon = 0.001; int i =0;
        for(i=0; i<dim;i++) {
            if(coor[i] - intpoint[i] <  epsilon) continue;
            else break;
        } if (i == dim) {data.weight = wei;break;}
    } if (p == n) std::cout << "CompElement::ComputeRequiredData: Weight Unknown";*/

    ShapeFunctions(intpoint,data.phi,data.dphidksi);

    geoel->GradX(intpoint,data.x,data.gradx); std::cout << "\n\nGradX\n\n"; data.gradx.Print();
    geoel->X(intpoint,data.x);
    Matrix jacinv(dim,dim,0), jac(dim,dim,0);
    geoel->Jacobian(data.gradx,jac,data.axes,data.detjac,jacinv); std::cout << "\n\nJacInv\n\n";jacinv.Print();
    //std::cout << "\n\nCompElement::Convert2Axes: Print dphi\n"; data.dphidksi.Print(std::cout); std::cout << "\n\n";
    Convert2Axes(data.dphidksi, jacinv,data.dphidx); //std::cout << "\n\ndphidksi: \n"; data.dphidksi.Print(std::cout);std::cout << "\n\ndphidx: \n"; data.dphidx.Print(std::cout);

    //GetMultiplyingCoeficients(data.coefs);
    //data.ComputeSolution(); //fill solution, dsoldksi and dsoldx
}

void CompElement::Convert2Axes(const Matrix &dphi, const Matrix &jacinv, Matrix &dphidx) const {
    //Convert parameter shape derivative into x,y,z coordinates shape derivative;
    int dim = dphi.Cols(); dphidx.Zero();
    for (int i = 0; i < dphi.Rows(); i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                dphidx(i, j) += dphi.GetVal(i, k) * jacinv.GetVal(j, k);
            }
        }
    }
}

void CompElement::CalcStiff(Matrix &ek, Matrix &ef) const {
    int nshape = NShapeFunctions();
    int dim = GetStatement()->Dimension();
    int nPoints = intrule->NPoints();
    int ndof = NDOF();

    ek.Resize(ndof,ndof); ef.Resize(ndof,1);
    for(int i =0; i< ndof; i++) { ef(i,0) = 0; for(int j = 0; j < ndof; j++) ek(i,j) = 0;}

    VecDouble intPoint(dim,0);
    double weights(0);
    IntPointData data[nPoints];
    for(int i = 0; i< nPoints ; i++){
        InitializeIntPointData(data[i]);
        intrule->Point(i,intPoint,weights);
        ComputeRequiredData(data[i],intPoint);
        if (data[i].detjac <= 0) std::cout <<"Warning: Non-Positive Jacobian determinant\n";
        GetStatement()->Contribute(data[i], weights*data[i].detjac, ek,ef);
    }
}

void CompElement::EvaluateError(std::function<void(const VecDouble &loc, VecDouble &val, Matrix &deriv) > fp, VecDouble &errors) const {

    MathStatement *material = GetStatement();
    IntRule *intRuleError = GetIntRule();
    intRuleError->SetOrder(15); //15 is the maximum order

    int numErrors =  2; //errors.size();
    // 0: Energy Norm;
    // 1: Zero Norm;
    // 2: Infinite Norm
    errors.resize(3);
    for(int i = 0;i <numErrors; i++) errors[i] = 0.0;

    int nstate = material->NState();
    int dim = Dimension();
    VecDouble u(nstate), val(numErrors);
    Matrix du(nstate,dim);

    IntPointData data;
    InitializeIntPointData(data);
    double weight = 0;
    int nintrulepoints = intRuleError->NPoints();

    for(int i = 0; i<nintrulepoints ; i++){
        //if(this->GetStatement()->GetMatID() != 1) continue;
        intRuleError->Point(i,data.ksi,data.weight);
        ComputeRequiredData(data,data.ksi);
        GetMultiplyingCoeficients(data.coefs);
        data.ComputeSolution();
        if(data.detjac < 0) std::cout << "CompElement::EvaluateError: Waring: detjac is not strictly positive";
        weight = data.weight*data.detjac;

        if(fp) { fp(data.x,u,du); };

        mat->ContributeError(data,u,du,val);

        for(int ier = 0; ier < numErrors; ier++) errors[ier] += val[ier]*weight;
        if (val[2] < 0) {if (-1*val[2]>errors[2]) errors[2] = -1*val[2];}
        else if(val[2] > errors[2]) errors[2] = val[2];
    }

    for(int ier = 0; ier < numErrors ; ier++) {
        errors[ier] = sqrt(errors[ier]);
    }
}

void CompElement::Solution(VecDouble &intpoint, int var, VecDouble &sol) const {
    IntPointData data;
    InitializeIntPointData(data);
    ComputeRequiredData(data,intpoint);

    GetMultiplyingCoeficients(data.coefs);
    data.ComputeSolution();
    sol.resize(GetStatement()->Dimension());

    GetStatement()->PostProcessSolution(data,var,sol);
}
