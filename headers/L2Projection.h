//
//  L2Projection.h
//  FemCoursePZ
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef L2Projection_h
#define L2Projection_h

#include "MathStatement.h"
#include "DataTypes.h"
#include  "IntPointData.h"
#include <functional>

class L2Projection : public MathStatement
{
    // Boundary condition ID
    int BCType = 0;
    
    // L2 projection matrix
    Matrix projection;
    
    // First value of boundary condition
    Matrix BCVal1;
    
    // Second value of boundary condition
    Matrix BCVal2;
    
    // Force funtion related to L2 projection math statement
    std::function<void(const VecDouble &co, VecDouble &result)> forceFunction;
    
    std::function<void(const VecDouble &loc, VecDouble &result, Matrix &deriv)> SolutionExact;
    
public:
    
    enum PostProcVar {ENone, ESol, EDSol};
    
    // Default constructor of L2Projection
    L2Projection();
    
    // Constructor of L2Projection
    L2Projection(int bctype, int materialid, Matrix &proj, Matrix Val1, Matrix Val2);
    
    // Copy constructor of L2Projection
    L2Projection(const L2Projection &copy);
    
    // Operator of copy
    L2Projection &operator=(const L2Projection &copy);
    
    int GetBCType() const { return BCType; }
    
    // Method for creating a copy of the element
    virtual L2Projection *Clone() const;
    
    // Default contructor of L2Projection
    virtual ~L2Projection();
    
    // Return the L2 projection matrix
    Matrix GetProjectionMatrix() const;
    
    // Set the L2 projection matrix
    void SetProjectionMatrix(const Matrix &proj);

    // Return Val1
    Matrix Val1() const    //I modified this method: Inconsistent arguments;
    {
        return BCVal1;
    }
    
    // Return Val2
    Matrix Val2() const
    {
        return BCVal2;
    }
    
    // Return the force function related to L2 projection math statement
    std::function<void(const VecDouble &co, VecDouble &result)> GetForceFunction() const
    {
        return forceFunction;
    }

    // Set the force function related to L2 projection math statement
    void SetForceFunction(const std::function<void(const VecDouble &co, VecDouble &result)> &f)
    {
        forceFunction = f;
    }
    
    // Set the exact solution related to L2 projection math statement
    void SetExactSolution(const std::function<void(const VecDouble &loc, VecDouble &result, Matrix &deriv)> &Exact)
    {
        SolutionExact = Exact;
    }
    
    // Return the number of state variables
    virtual int NState() const {
        return 1;
    };
    
    // Return the number of errors
    virtual int NEvalErrors() const;
    
    virtual int VariableIndex(const PostProcVar var) const;
    
    // Return the variable index associated with the name
    virtual PostProcVar VariableIndex(const std::string &name);
    
    // Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
    virtual int NSolutionVariables(const PostProcVar var);
    
    // Method to implement integral over element's volume
    virtual void Contribute(IntPointData &integrationpointdata, double weight, Matrix &EK, Matrix &EF) const;
    
    // Method to implement error over element's volume
    virtual void ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const;
    
    // Prepare and print post processing data
    virtual void PostProcessSolution(const IntPointData &integrationpointdata, const int var, VecDouble &sol) const;
    
    
};
#endif /* L2Projection_h */