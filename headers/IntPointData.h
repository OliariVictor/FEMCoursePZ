//
//  IntPointData.h
//  FemCoursePZ
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef IntPointData_h
#define IntPointData_h

#include "DataTypes.h"
#include "tpanic.h"

class IntPointData
{
public:
    
    VecDouble ksi;
    
    double weight;
    
    VecDouble phi;
    
    Matrix dphidksi;
    
    VecDouble x;
    
    Matrix gradx;
    
    Matrix axes;
    
    double detjac;
    
    Matrix dphidx;
    
    VecDouble solution;
    
    Matrix dsoldksi;
    
    Matrix dsoldx;
    
    VecDouble coefs;
    
    void ComputeSolution()
    {
        if(coefs.size()%phi.size())
        {
        //    DebugStop();
        }
        int nstate = coefs.size()/phi.size();
        solution.resize(nstate);
        for(auto &solval:solution) solval = 0.;
        dsoldksi.Resize(dphidx.Cols(), nstate); // Rows replaced by Cols
        dsoldx.Resize(dphidx.Cols(), nstate);   // Rows replaced by Cols
        dsoldx.Zero();
        dsoldksi.Zero();
        int dim = dphidx.Cols();                      // Rows replaced by Cols
        for (int iphi=0; iphi<phi.size(); iphi++) {
            double phival = phi[iphi];
            for (int istate=0; istate<nstate; istate++) {
                solution[istate] += phival*coefs[iphi*nstate+istate]; if(iphi ==0) {std::cout << "\n\nksi\n\n";for(int j=0;j<ksi.size();j++) std::cout << ksi[j] <<std::endl;std::cout << "\n\nX\n\n";for(int j=0;j<x.size();j++) std::cout << x[j] <<std::endl; std::cout << "\n\ndphidx\n\n";dphidx.Print();std::cout << "\n\ndphidksi\n\n";dphidksi.Print();std::cout << "\n\ngradX\n\n";gradx.Print();std::cout << "\n\nCoefs\n\n";for(int j=0;j<coefs.size();j++) std::cout << coefs[j] <<std::endl;}
                for (int d=0; d < dim; d++) {
                    dsoldksi(d,istate) += coefs[iphi*nstate+istate]*dphidksi(iphi,d); // (d,iphi) replaced by (iphi,d)
                    dsoldx(d,istate) += coefs[iphi*nstate+istate]*dphidx(iphi,d) ;     // (d,iphi) replaced by (iphi,d)
                }
            }
        } std::cout << "\n\ndsoldx\n"; dsoldx.Print(std::cout);
    }
    
};
#endif /* IntPointData_h */