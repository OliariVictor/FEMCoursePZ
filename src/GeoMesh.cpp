/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoMesh.h"
#include "tpanic.h"
#include "GeoElementSide.h"
#include "GeoElement.h"
#include "GeoNode.h"

GeoMesh::GeoMesh(): Nodes(0),Elements(0), Reference(0), fDim(0) {
}

GeoMesh::GeoMesh(const GeoMesh &copy) {
    Nodes = copy.Nodes;
    Elements = copy.Elements;
    Reference = copy.Reference;
    fDim = copy.fDim;
}

GeoMesh &GeoMesh::operator=(const GeoMesh &copy) {
    Nodes = copy.Nodes;
    Elements = copy.Elements;
    Reference = copy.Reference;
    fDim = copy.fDim;
    return *this;
}

void GeoMesh::SetNumNodes(int nnodes) {
    Nodes.resize(nnodes);
}

void GeoMesh::SetNumElements(int numelements) {
    Elements.resize(numelements);
}

int GeoMesh::NumNodes() {
    return(Nodes.size());
}

int GeoMesh::NumElements() {
    return(Elements.size());
}

GeoNode &GeoMesh::Node(int node) {
    return Nodes[node];
}

void GeoMesh::SetElement(int elindex, GeoElement *gel) {
    Elements[elindex] = gel;
}

GeoElement *GeoMesh::Element(int elindex) {
    return(Elements[elindex]);
}

void GeoMesh::BuildConnectivity() {
    VecInt  SideNum(NumNodes(),-1);
    std::vector<GeoElement*> NeighNode(NumNodes(),0);
    int nelem = NumElements();
    int elem_index;
    for(elem_index = 0; elem_index < nelem ; elem_index++){ //Setting the neighbourhood of all nodes
        GeoElement *gel = Elements[elem_index];
        if (!gel) continue; //if gel points NULL -> go for the next loop iteration
        int ncor = gel->NCornerNodes();
        int in;
        for (in = 0 ; in < ncor ; in++) { // Iterate through all gel's corners
            int nod = gel->NodeIndex(in);
            if (SideNum[nod] == -1) {     //If this node index is not initialized -> set self as the node neighbour
                NeighNode[nod] = gel;
                SideNum[nod] = in;
                GeoElementSide gelSide(gel,in);
                if (gel->SideIsUndefined(in)) gel->SetNeighbour(in,gelSide);
            } else{ // If node index is already initialized -> set a neighbour connection
                GeoElementSide neigh(NeighNode[nod],SideNum[nod]);
                GeoElementSide gelSide(gel,in);
                if(!neigh.IsNeighbour(gelSide)){
                    neigh.IsertConnectivity(gelSide);
                }
            }
        }
    }
    //Insert neighbourhood for higher dimension sides...
    for(elem_index = 0; elem_index < nelem ; elem_index++){
        GeoElement *gel = Elements[elem_index];
        if (!gel) continue;
        int ncor = gel->NCornerNodes();
        int nside = gel->NSides();
        //Iterating through every side which it's dimension is higher then 1...
        for(int i = ncor; i < nside; i++){
            if (gel->SideIsUndefined(i)){ //if this side's neighbourhood is not yet defined...
                //Setting neighbour = self...
                GeoElementSide gelSide(gel,i);
                gel->SetNeighbour(i,gelSide);

                //Computing this side's neighbourhood...
                std::vector<GeoElementSide> neighbour;
                gelSide.ComputeNeighbours(neighbour);// if(gelSide.Side() == 5) {gelSide.Element()->Print(std::cout);}
                int nneigh = neighbour.size();

                //Iterating through each neighbour, checking its consistency and creating it's connectivity
                for(int j = 0; j < nneigh ; j++ ){
                    if(neighbour[j].Side() == -1){
                        std::cout << "Inconsistent mesh" << std::endl;
                        neighbour[j].Element()->Print(std::cout);
                        std::cout << std::endl << neighbour[j].Side(); continue;
                    }
                    if(neighbour[j].Element()->SideIsUndefined(neighbour[j].Side())){
                        neighbour[j].Element()->SetNeighbour(neighbour[j].Side(),neighbour[j]);
                    }
                    gelSide.IsertConnectivity(neighbour[j]);
                }
            }
        }
    }
}

void GeoMesh::Print(std::ostream &out) {
    out << "\n\t\t GEOMETRIC INFORMATIONS:\n\n";
    out << "number of nodes               = " << Nodes.size() << "\n";
    out << "number of elements            = " << Elements.size() << "\n";

    out << "\n\tGeometric Node Information:\n\n";
    int64_t i;
    int64_t nnodes = Nodes.size();
    for(i=0; i<nnodes; i++)
    {
        out << "Index: " << i << " ";
        Nodes[i].Print(out);
    }
    out << "\n\tGeometric Element Information:\n\n";
    int64_t nelem = Elements.size();
    for(i=0; i<nelem; i++)
    {
        if(Elements[i]) Elements[i]->Print(out);
        out << "\n";
    }
}
