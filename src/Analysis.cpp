/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Analysis.h"
#include "tpanic.h"
#include "Assemble.h"
#include "CompMesh.h"
#include "VTKGeoMesh.h"
#include "CompElement.h"
#include "GeoElement.h"
#include "GeoMesh.h"
#include "MathStatement.h"


using namespace std;

Analysis::Analysis(): cmesh(), Solution(0,0), GlobalSystem(0,0), RightHandSide(0,0) {
}

Analysis::Analysis(const Analysis &cp) {
    cmesh = cp.cmesh;
    Solution = cp.Solution;
    GlobalSystem = cp.GlobalSystem;
    RightHandSide = cp.RightHandSide;
}

Analysis &Analysis::operator=(const Analysis &cp) {
    cmesh = cp.cmesh;
    Solution = cp.Solution;
    GlobalSystem = cp.GlobalSystem;
    RightHandSide = cp.RightHandSide;
    return(*this);
}

Analysis::~Analysis() {
}

Analysis::Analysis(CompMesh *mesh): Analysis() {
    cmesh = mesh;
}

void Analysis::SetMesh(CompMesh *mesh) {
    cmesh = mesh;
}

CompMesh *Analysis::Mesh() const {
    return cmesh;
}

void Analysis::RunSimulation() {
    Matrix K,F;
    Assemble ass(cmesh);
    ass.Compute(K,F); std::cout << "\nGlobK: \n"; K.Print(std::cout);std::cout << "\n\nGlobF: \n";F.Print(std::cout);

    GlobalSystem = K;
    RightHandSide = F;

    K.Solve_LU(F);

    Solution = F; std::cout << "\n\nComputedSol\n\n"; F.Print(std::cout);
    VecDouble sol(F.Rows(),0); for(int i =0; i< F.Rows() ; i++) sol[i] = F(i,0);
    cmesh->LoadSolution(sol);
}

void Analysis::PostProcessSolution(const std::string &filename, PostProcess &defPostProc) const {
    VTKGeoMesh::PrintSolVTK(cmesh,defPostProc,filename);
}

VecDouble Analysis::PostProcessError(std::ostream &out, PostProcess &defPostProc) const {
    VecDouble uex;
    Matrix duex;
    VecDouble errors(3,0),val(3,0);
    int numElem = cmesh->GetElementVec().size();

    for(int64_t elemInd = 0; elemInd < numElem; elemInd++){
        CompElement *cel = cmesh->GetElement(elemInd);  if(cel->GetStatement()->GetMatID() != 1) continue;
        cel->EvaluateError(defPostProc.GetExact(),val); std::cout << "\n\ncel : " << cel->GetIndex();
        for(int i =0; i < 2 ; i++){
            errors[i] += val[i]*val[i];
        }
        if(abs(val[2]) > errors[2]) errors[2] = abs(val[2]);
    }
    for (int i =0; i< 2; i++) errors[i] = sqrt(errors[i]);

    out << "--------------ERRORS-------------" << endl;
    out <<"Energy norm -> sqrt[integral[e²+de²]] = "  << errors[0] << endl;
    out <<"L2 Norm -> sqrt[integral[e²]] = "     << errors[1] << endl;
    out << "Infinite Norm -> max[abs[e]] = "    << errors[2]  <<endl;

    return(errors);
}

