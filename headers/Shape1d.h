//
//  Shape1d.h
//  FEMCoursePZ
//
//  Created by Philippe Devloo on 03/04/18.
//

#ifndef Shape1d_h
#define Shape1d_h

#include "DataTypes.h"
#include "Topology1d.h"

class Shape1d : public Topology1d
{
public:
    // Computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
    static void Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi);

    
    // Returns the number of shape functions associated with a side
    static int NShapeFunctions(int side, int order);
    
    // Returns the total number of shape functions
    static int NShapeFunctions(VecInt &orders);

};

#endif /* Shape1d_h */