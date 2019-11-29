/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "IntRule1d.h"
#include "tpanic.h"
#include <cmath>

using namespace std;

#define PI 3.141592654

IntRule1d::IntRule1d() : IntRule()
{
}

IntRule1d::IntRule1d(int order) : IntRule()
{
    SetOrder(order);
}

void IntRule1d::gauleg(const double x1, const double x2, VecDouble &co, VecDouble &w){
    {
        const double EPS=1.0e-10;
        double z1,z,xm,xl,pp,p3,p2,p1;
        int n=co.size();
        int m=(n+1)/2;
        xm=0.5*(x2+x1);
        xl=0.5*(x2-x1);
        for (int i=0;i<m;i++) {
            z=cos(3.141592654*(i+0.75)/(n+0.5));
            do {
                p1=1.0;
                p2=0.0;
                for (int j=0;j<n;j++) {
                    p3=p2;
                    p2=p1;
                    p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
                }
                pp=n*(z*p1-p2)/(z*z-1.0);
                z1=z;
                z=z1-p1/pp;
            } while (abs(z-z1) > EPS);
            co[i]=xm-xl*z;
            co[n-1-i]=xm+xl*z;
            w[i]=2.0*xl/((1.0-z*z)*pp*pp);
            w[n-1-i]=w[i];
        }
    }
}

void IntRule1d::SetOrder(int order){
    fOrder = order;
    //int IntOrder = IntRule::NPoints();
    int nPoints;

    //If fOrder is even...Else...
    if (fOrder%2 == 0){
        nPoints = (int)(fOrder/2)+1;
    }
    else{
        nPoints = (int)((fOrder-1)/2+1);
    }
    if (nPoints <1) DebugStop();
    switch(nPoints)
    {
        case 1:
            fPoints.Resize(1,1);
            fWeights.resize(1);
            fPoints(0,0) = 0.;
            fWeights[0] = 2.;
        break;
        case 2:
            fPoints.Resize(2,1);
            fWeights.resize(2);

            fPoints(0,0) = -1/sqrt(3);
            fPoints(1,0) = -fPoints(0,0);
            fWeights[0] = 1.0;
            fWeights[1] = fWeights[0];
            break;
        case 3:
            fPoints.Resize(3,1);
            fWeights.resize(3);

            fPoints(0,0) = -sqrt(3./5.);
            fPoints(1,0) = 0;
            fPoints(2,0) = -fPoints(0,0);

            fWeights[0] = 5./9.;
            fWeights[1] = 8./9.;
            fWeights[2] = fWeights[0];
            break;

        case 4:
            fPoints.Resize(4,1);
            fWeights.resize(4);

            fPoints(0,0) = -.8611363116;
            fPoints(1,0) = -.3399810436;
            fPoints(2,0) = -fPoints(1,0);
            fPoints(3,0) = -fPoints(0,0);

            fWeights[0] = 0.3478548451;
            fWeights[1] = -.3399810436;
            fWeights[2] = fWeights[1];
            fWeights[3] = fWeights[0];
            break;
        default:
            fPoints.Resize(nPoints,1);
            fWeights.resize(nPoints);

            VecDouble coordinates(nPoints), weights(nPoints);
            gauleg(-1,1,coordinates,weights);
            for (int i =0 ; i< nPoints; i++) {
                fPoints(i,0) = coordinates[i];
                fWeights[i] = weights[i];
            }
            break;
    }
}
