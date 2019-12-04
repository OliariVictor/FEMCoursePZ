/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoElementTemplate.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "tpanic.h"
#include "GeoMesh.h"

using namespace std;

template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh, int index) : GeoElement(materialid,gmesh,index), Geom() {
    Geom.SetNodes(nodeindices);
    MaterialId = materialid;
}

template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const GeoElementTemplate &copy) : GeoElement(copy.MaterialId, copy.GMesh, copy.Index) {
    Reference = copy.Reference;
}

template<class TGeom>
GeoElementTemplate<TGeom> &GeoElementTemplate<TGeom>::operator=(const GeoElementTemplate &copy){
    GeoElement::MaterialId = copy.MaterialId;
    GeoElement::GMesh = copy.GMesh;
    GeoElement::Index = copy.Index;
    GeoElement::Reference = copy.Reference;
    return *this;
}

template<class TGeom>
ElementType GeoElementTemplate<TGeom>::Type() {
    return(TGeom::Type());
}

template<class TGeom>
void GeoElementTemplate<TGeom>::X(const VecDouble &xi, VecDouble &x) {
    //Computing the coordinates of each corner;
    VecInt ind;Geom.GetNodes(ind);
    GeoMesh *mesh = GetMesh();
    int dim = GetMesh()->Dimension();
    Matrix nodeCoord(ind.size(),dim,0);
    GeoNode node;

    for(int i = 0; i < ind.size() ; i++ ){
        node = mesh->Node(ind[i]);
        for(int j =0; j <dim; j++) {
            nodeCoord(i,j) = node.Coord(j);
        }
    } //std::cout << "\n Elem Ind" << GetIndex() << "\n";
    Geom.X(xi,nodeCoord,x);
}

template<class TGeom>
void GeoElementTemplate<TGeom>::GradX(const VecDouble &xi, VecDouble &x, Matrix &gradx) {
    VecInt ind;Geom.GetNodes(ind);
    GeoMesh *mesh = GetMesh();
    int dim = this->GetMesh()->Dimension();
    Matrix nodeCoord(ind.size(),dim,0);
    GeoNode node;

    int i =0;
    for(auto k: ind ){
        node = mesh->Node(k);
        for(int j =0; j <dim; j++) {
            nodeCoord(i,j) = node.Coord(j);
        }i++;
    }

    Geom.GradX(xi,nodeCoord,x,gradx);
}

template<class TGeom>
void GeoElementTemplate<TGeom>::Jacobian(const Matrix &gradx, Matrix &jac, Matrix &axes, double &detjac, Matrix &jacinv) {
    detjac = 0.0;
    int nrows = gradx.Rows();
    int ncols = gradx.Cols();
    int meshDim = this->GMesh->Dimension();//int meshDim   = ncols;
    int elemDim = this->Reference->Dimension(); //int elemDim = nrows;

    VecDouble v_1(meshDim,0),v_2(meshDim,0);
    VecDouble v_1_til(3,0),v_2_til(3,0);
    double norm_v_1 = 0.;
    double norm_v_1_til = 0.0;
    double norm_v_2_til = 0.0;
    double v_1_dot_v_2  = 0.0; Matrix res;

    switch(elemDim){
        case 0:
            jac.Resize(elemDim,elemDim);
            axes.Resize(elemDim,meshDim);
            jacinv.Resize(elemDim,elemDim);
            detjac = 1;
            break;
        case 1:
            jac.Resize(elemDim,elemDim);
            axes.Resize(meshDim,elemDim);
            jacinv.Resize(elemDim,elemDim);
            jac.Zero();

            /**  Definitions: v1 -> is the xi_direction of the Gradient*/
            for (int i = 0; i < meshDim; i++) {
                v_1[i]  = gradx.GetVal(0,i);
            }

            for(int i = 0; i < meshDim; i++) {
                norm_v_1 += v_1[i]*v_1[i];
            }

            norm_v_1    = sqrt(norm_v_1);
            jac(0,0)    = norm_v_1;
            detjac      = norm_v_1;

            if(detjac < 0.00001 && detjac > -0.000001) { std::cout << "Null Jacobian"; DebugStop();}

            jacinv(0,0) = 1.0/detjac;

            for(int i=0; i < meshDim; i++) {
                axes(i,0) = v_1[i]/norm_v_1;
            }
            axes.Multiply(jacinv,res,0); jacinv = res;
            break;

        case 2:
            jac.Resize(elemDim,elemDim);
            axes.Resize(meshDim,elemDim);
            jacinv.Resize(elemDim,elemDim);
            jac.Zero();

            /**  Definitions: v1 -> is the xi_direction of the Gradient, v2 -> is the eta_direction of the Gradient*/
            /**  Definitions: v_1_til and v_2_til -> asscoiated orthonormal vectors to v_1 and v_2*/

            for (int i = 0; i < nrows; i++) {
                v_1[i]  = gradx.GetVal(i,0);
                v_2[i]  = gradx.GetVal(i,1);
            }

            for(int i = 0; i < meshDim; i++) {
                norm_v_1_til    += v_1[i]*v_1[i];
                v_1_dot_v_2     += v_1[i]*v_2[i];
            }
            norm_v_1_til = sqrt(norm_v_1_til);

            for(int i=0 ; i < meshDim; i++) {
                v_1_til[i]          = v_1[i] / norm_v_1_til; // Normalizing
                v_2_til[i]          = v_2[i] - v_1_dot_v_2 * v_1_til[i] / norm_v_1_til;
                norm_v_2_til   += v_2_til[i]*v_2_til[i];
            }
            norm_v_2_til = sqrt(norm_v_2_til);

            jac(0,0) = norm_v_1_til;
            jac(0,1) = v_1_dot_v_2/norm_v_1_til;
            jac(1,1) = norm_v_2_til;

            detjac = jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);

            jacinv(0,0) = +jac(1,1)/detjac;// if(v_1)
            jacinv(1,1) = +jac(0,0)/detjac;
            jacinv(0,1) = -jac(0,1)/detjac;
            jacinv(1,0) = -jac(1,0)/detjac;

            if(detjac < 0.00001 && detjac > -0.000001) { std::cout << "Null Jacobian"; DebugStop();}

            for(int i=0; i < meshDim; i++) {
                v_2_til[i] /= norm_v_2_til; // Normalizing
                axes(i,0)  = v_1_til[i];
                axes(i,1)  = v_2_til[i];
            }
            axes.Multiply(jacinv,res,0); jacinv = res;
             //jacinv.Multiply(axes,res,0); jacinv = res;
            break;jacinv = res;

        case 3:
            jac.Resize(3,3);
            axes.Resize(3,3);
            jacinv.Resize(3,3);

            jac = gradx;
            detjac -= jac(0,2)*jac(1,1)*jac(2,0);//- a02 a11 a20
            detjac += jac(0,1)*jac(1,2)*jac(2,0);//+ a01 a12 a20
            detjac += jac(0,2)*jac(1,0)*jac(2,1);//+ a02 a10 a21
            detjac -= jac(0,0)*jac(1,2)*jac(2,1);//- a00 a12 a21
            detjac -= jac(0,1)*jac(1,0)*jac(2,2);//- a01 a10 a22
            detjac += jac(0,0)*jac(1,1)*jac(2,2);//+ a00 a11 a22
            //std::cout << "\n\nJacobian: \n";jac.Print(std::cout);
            if(detjac == 0){
                cout << "Singular Jacobian, determinant of jacobian = " << detjac << std::endl; DebugStop();
            }

            jacinv(0,0) = (-jac(1,2)*jac(2,1)+jac(1,1)*jac(2,2))/detjac;//-a12 a21 + a11 a22
            jacinv(0,1) = ( jac(0,2)*jac(2,1)-jac(0,1)*jac(2,2))/detjac;//a02 a21 - a01 a22
            jacinv(0,2) = (-jac(0,2)*jac(1,1)+jac(0,1)*jac(1,2))/detjac;//-a02 a11 + a01 a12
            jacinv(1,0) = ( jac(1,2)*jac(2,0)-jac(1,0)*jac(2,2))/detjac;//a12 a20 - a10 a22
            jacinv(1,1) = (-jac(0,2)*jac(2,0)+jac(0,0)*jac(2,2))/detjac;//-a02 a20 + a00 a22
            jacinv(1,2) = ( jac(0,2)*jac(1,0)-jac(0,0)*jac(1,2))/detjac;//a02 a10 - a00 a12
            jacinv(2,0) = (-jac(1,1)*jac(2,0)+jac(1,0)*jac(2,1))/detjac;//-a11 a20 + a10 a21
            jacinv(2,1) = ( jac(0,1)*jac(2,0)-jac(0,0)*jac(2,1))/detjac;//a01 a20 - a00 a21
            jacinv(2,2) = (-jac(0,1)*jac(1,0)+jac(0,0)*jac(1,1))/detjac;//-a01 a10 + a00 a11
            //std::cout <<"\n\nJacobianInv: \n";jacinv.Print(std::cout);
            for(int i =0; i<3; i++) for(int j = 0; j<3;j++) axes(i,j) = 0;
            axes(0,0) = axes(1,1) = axes(2,2) = 1;
            jacinv.Transpose(res);jacinv = res; //std::cout <<"\n\nJacobianInvTransposed: \n";jacinv.Print(std::cout);
            break;

        default:
            cout << "Please insert a valid dimension number";
            DebugStop();
    }
}
//The previously implemented Jacobian is only useful for elements aligned with XYZ coordinates. Hence, it doesn't suit many boundary elements. For example, if a tetrahedron face which does not belong to the planes XY, XZ,YZ is a boundary element, the following Jacobian is not valid.
/*template<class TGeom>
void GeoElementTemplate<TGeom>::Jacobian(const Matrix &gradx, Matrix &jac, Matrix &axes, double &detjac, Matrix &jacinv) {
    detjac = 0.0;
    int nrows = gradx.Rows();
    int ncols = gradx.Cols();
    int dim   = nrows;
    TMatrix temp;

    switch(dim){
        case 0:
            jac.Resize(0,0);
            axes.Resize(0,0);
            jacinv.Resize(0,0);
            detjac = 1;
        case 1:
            jac.Resize(1,1);
            axes.Resize(1,1);
            jacinv.Resize(1,1);

            temp = gradx;      // This operation is necessary because there is no const TMatrix (int, int) defined.
            jac(0,0) = temp(0,0);
            detjac = temp(0,0);
            axes(0,0) = 1;
            jacinv(0,0) = 1.0/jac(0,0);

            if(detjac == 0){
                cout << "Singular Jacobian, determinant of jacobian = " << detjac << std::endl; DebugStop();
            }break;
        case 2:
            jac.Resize(2,2);
            axes.Resize(2,2);
            jacinv.Resize(2,2);

            jac = gradx;
            axes(0,0) = axes(1,1) = 1;
            axes(0,1) = axes(1,0) = 0;
            detjac = jac(0,0)*jac(1,1) - jac(0,1)*jac(1,0);

            jacinv(0,0) = (1.0/detjac)*jac(1,1);
            jacinv(0,1) = -(1.0/detjac)*jac(0,1) ;
            jacinv(1,0) = -(1.0/detjac)*jac(1,0) ;
            jacinv(1,1) = (1.0/detjac)*jac(0,0) ;

            if(detjac == 0){
                cout << "Singular Jacobian, determinant of jacobian = " << detjac << std::endl; DebugStop();
            }break;
        case 3:
            jac.Resize(3,3);
            axes.Resize(3,3);
            jacinv.Resize(3,3);

            jac = gradx;
            detjac -= jac(0,2)*jac(1,1)*jac(2,0);//- a02 a11 a20
            detjac += jac(0,1)*jac(1,2)*jac(2,0);//+ a01 a12 a20
            detjac += jac(0,2)*jac(1,0)*jac(2,1);//+ a02 a10 a21
            detjac -= jac(0,0)*jac(1,2)*jac(2,1);//- a00 a12 a21
            detjac -= jac(0,1)*jac(1,0)*jac(2,2);//- a01 a10 a22
            detjac += jac(0,0)*jac(1,1)*jac(2,2);//+ a00 a11 a22

            if(detjac == 0){
                cout << "Singular Jacobian, determinant of jacobian = " << detjac << std::endl; DebugStop();
            }

            jacinv(0,0) = (-jac(1,2)*jac(2,1)+jac(1,1)*jac(2,2))/detjac;//-a12 a21 + a11 a22
            jacinv(0,1) = ( jac(0,2)*jac(2,1)-jac(0,1)*jac(2,2))/detjac;//a02 a21 - a01 a22
            jacinv(0,2) = (-jac(0,2)*jac(1,1)+jac(0,1)*jac(1,2))/detjac;//-a02 a11 + a01 a12
            jacinv(1,0) = ( jac(1,2)*jac(2,0)-jac(1,0)*jac(2,2))/detjac;//a12 a20 - a10 a22
            jacinv(1,1) = (-jac(0,2)*jac(2,0)+jac(0,0)*jac(2,2))/detjac;//-a02 a20 + a00 a22
            jacinv(1,2) = ( jac(0,2)*jac(1,0)-jac(0,0)*jac(1,2))/detjac;//a02 a10 - a00 a12
            jacinv(2,0) = (-jac(1,1)*jac(2,0)+jac(1,0)*jac(2,1))/detjac;//-a11 a20 + a10 a21
            jacinv(2,1) = ( jac(0,1)*jac(2,0)-jac(0,0)*jac(2,1))/detjac;//a01 a20 - a00 a21
            jacinv(2,2) = (-jac(0,1)*jac(1,0)+jac(0,0)*jac(1,1))/detjac;//-a01 a10 + a00 a11

            for(int i =0; i<3; i++) for(int j = 0; j<3;j++) axes(i,j) = 0;
            break;

        default:
            cout << "Please insert a valid dimension number";
            DebugStop();
    }
}*/

template<class TGeom>
int GeoElementTemplate<TGeom>::WhichSide(VecInt &SideNodeIds) {
    int numberNodes = SideNodeIds.size();
    int num = NSides();
    for(int side = 0; side <num ; side++){

        if(NSideNodes(side) == 2 && numberNodes ==2){
            int n1 = NodeIndex(SideNodeIndex(side,0)), n2 = NodeIndex(SideNodeIndex(side,1));
            if((n1 == SideNodeIds[0] && n2 == SideNodeIds[1]) ||
               (n2 == SideNodeIds[0] && n1 == SideNodeIds[1]))    return side;
        }
        if(NSideNodes(side) ==1 && numberNodes ==1){
            if(NodeIndex(SideNodeIndex(side,0)) == NodeIndex(SideNodeIds[0])) return side;
        }
        if(NSideNodes(side) ==3 && numberNodes ==3) {
            int sni[3], snx[3], k;
            for (k = 0; k < 3; k++) snx[k] = NodeIndex(SideNodeIndex(side, k));//el atual
            for (k = 0; k < 3; k++) sni[k] = SideNodeIds[k];//el viz
            for (k = 0; k < 3; k++) {
                if (snx[0] == sni[k] && snx[1] == sni[(k + 1) % 3] && snx[2] == sni[(k + 2) % 3]) return side;
                if (snx[0] == sni[k] && snx[1] == sni[(k + 2) % 3] && snx[2] == sni[(k + 1) % 3]) return side;
            }
        }
        else if(NSideNodes(side) == 4 && numberNodes == 4) {//face quadrilateral
            int64_t sni[4],snx[4],k;
            for(k=0;k<4;k++) snx[k] = NodeIndex(SideNodeIndex(side,k));//el atual
            for(k=0;k<4;k++) sni[k] = SideNodeIds[k];//vizinho
            if(snx[0]==sni[0]) {
                for(k=1;k<4;k++) {
                    if(snx[1]==sni[k] && snx[2]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
                    if(snx[1]==sni[k] && snx[2]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */
            } else if(snx[1]==sni[0]) {
                for(k=1;k<4;k++) {
                    if(snx[0]==sni[k] && snx[2]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
                    if(snx[0]==sni[k] && snx[2]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */                               /* 1023 1230 1302 , 1032 1203 1320 */
            } else if(snx[2]==sni[0]) {
                for(k=0;k<4;k++) {
                    if(snx[0]==sni[k] && snx[1]==sni[k%3+1]     && snx[3]==sni[(k+1)%3+1]) return side;
                    if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3+1] && snx[3]==sni[k%3+1])     return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */                               /* 2013 2130 2301 , 2031 2103 2310 */
            } else if(snx[3]==sni[0]) {
                for(k=0;k<4;k++) {
                    if(snx[0]==sni[k] && snx[1]==sni[k%3+1]     && snx[2]==sni[(k+1)%3+1]) return side;
                    if(snx[0]==sni[k] && snx[1]==sni[(k+1)%3+1] && snx[2]==sni[k%3+1])     return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */                               / * 3012 3120 3201 , 3021 3102 3210 * /
            }
        } else if(numberNodes<1 || numberNodes > 4) {
            int is;
            for (is=0; is<numberNodes; is++) {
                if (NSideNodes(is) == numberNodes) {
                    break;
                }
            }
            if (is != num) {
                std::cout << "TPZGeoEl::WhichSide must be extended\n";
                DebugStop();
            }
        }
    } return -1;
}

template<class TGeom>
void GeoElementTemplate<TGeom>::Print(std::ostream &out) {
    GeoElement::Print(out);
}

template<class TGeom>
int GeoElementTemplate<TGeom>::SideIsUndefined(int side) {
    if (side <0 || side >= NSides()){
        std::cout << "GeoElementTemplate<TGeom>::SideIsUndefined(int side): non-existent side"; DebugStop();
    }
    return(Neighbour(side).Side() == -1);
}


template class GeoElementTemplate<Geom1d>;
template class GeoElementTemplate<GeomQuad>;
template class GeoElementTemplate<GeomTriangle>;
template class GeoElementTemplate<GeomTetrahedron>;
