//
//  CompElement.h
//  FEMCoursePZ
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef CompElement_h
#define CompElement_h

#include "DataTypes.h"
#include "IntRule.h"
#include "IntPointData.h"
#include "PostProcess.h"
#include <functional>

class CompMesh;
class GeoElement;
class MathStatement;

class CompElement
{
    
    // Computational mesh to which the element belongs
    CompMesh *compmesh = 0;
    
    // Element index into mesh element vector
    int64_t index = -1;
    
    // Geometrical element associated
    GeoElement *geoel = 0;
    
    // Integration rule object
    IntRule *intrule = 0;

    // Material object associated
    MathStatement *mat = 0;
    
public:
    
    // Default constructor of CompElement
    CompElement();
    
    // Constructor of CompElement
    CompElement(int64_t ind, CompMesh *cmesh, GeoElement *geo);
    
    // Copy constructor of CompElement
    CompElement(const CompElement &copy);
    
    // Operator of copy
    CompElement &operator=(const CompElement &);
    
    // Destructor of CompElement
    ~CompElement();
    
    // Method for creating a copy of the element
    virtual CompElement *Clone() const;
    
    // Return the material object associated with the element
    MathStatement *GetStatement() const;
    
    // Set the material objectdouble weight; associated with the element
    void SetStatement( MathStatement *statement);
    
    // Return integration rule established
    IntRule *GetIntRule() const;
    
    // Set integration rule established
    void SetIntRule(IntRule *intrule);
    
    // Return element index
    int64_t GetIndex() const
    {
        return index;
    }
    
    // Set element index
    void SetIndex(int64_t ind);
    
    // Return the geometric element associated
    GeoElement *GetGeoElement() const;
    // Set the geometric element associated
    void SetGeoElement(GeoElement *element);

    // Return a pointer to the element computational mesh
    CompMesh *GetCompMesh() const;
    
    // Set a pointer to the element computational mesh
    void SetCompMesh(CompMesh *mesh);

    // Initialize integration points data object
    void InitializeIntPointData(IntPointData &data) const;
    
    // Compute and fill integration points data object
    void ComputeRequiredData(IntPointData &data, VecDouble &intpoint) const;
    
    // Convert a shapefunction derivative in xi-eta to a function derivative in axes
    void Convert2Axes(const Matrix &dphi, const Matrix &jacinv, Matrix &dphidx) const;
    
    // Compute the element stifness matrix and force vector
    virtual void CalcStiff(Matrix &ek, Matrix &ef) const;
    
    // Compute error and exact solution
    virtual void EvaluateError(std::function<void(const VecDouble &loc,VecDouble &val,Matrix &deriv)> fp,
                               VecDouble &errors) const;

    // Compute shape functions set at point x
    virtual void ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi) const = 0;
    
    // Compute the solution and its gradient at a parametric point
    // for dsol the row indicates the direction, the column indicates the state variable
    virtual void Solution(VecDouble &intpoint, int var, VecDouble &sol) const;
    
    // Get Multiplying Coeficients
    virtual void GetMultiplyingCoeficients(VecDouble &coefs) const = 0;
    
    // Return the number of shape functions
    virtual int NShapeFunctions() const = 0;
    
    // Set number of DOF
    virtual void SetNDOF(int64_t ndof) = 0;
    
    // Set DOF index in vector position i
    virtual void SetDOFIndex(int i, int64_t dofindex) = 0;
    
    // Get DOF index in vector position i
    virtual int64_t GetDOFIndex(int i) = 0;
    
    // Return the number of degree of freedom
    virtual int NDOF() const = 0;
    
    // Return the number of shape functions stored in the DOF data structure
    virtual int NShapeFunctions(int doflocindex) const = 0;
    
    // Use the Shape template class to compute the number of shape functions
    virtual int ComputeNShapeFunctions(int doflocindex, int order) = 0;
    // Return the dimension of the element
    virtual int Dimension() const = 0;
    // Method to print computational element information
    virtual void Print(std::ostream &out) = 0;
};
#endif /* CompElement_h */
