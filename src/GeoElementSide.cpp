/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <algorithm>
#include "GeoElementSide.h"
#include "tpanic.h"
#include <bits/stdc++.h>
#include "GeoMesh.h"

GeoElementSide::GeoElementSide():  fSide(-1) { fElement = 0;
}

GeoElementSide::GeoElementSide(const GeoElementSide &copy) {
    fElement = copy.fElement;
    fSide = copy.fSide;
}

GeoElementSide &GeoElementSide::operator=(const GeoElementSide &copy) {
    fElement = copy.fElement;
    fSide = copy.fSide;
    return *this;
}

GeoElementSide GeoElementSide::Neighbour() const {
    return fElement ? fElement->Neighbour(fSide) : GeoElementSide();
}

void GeoElementSide::SetNeighbour(const GeoElementSide &neighbour) {
    if(fElement){
        fElement->SetNeighbour(fSide,neighbour);
    }else {std::cout << "Insert a valid GeoElement pointer"; DebugStop();}
}

bool GeoElementSide::IsNeighbour(const GeoElementSide &candidate) {
    if (Exists() || candidate.Exists()) {std::cout << "GeoElementSide::IsNeighbour: Inconsistent mesh\n"; return false;}
    if (candidate == *this) return true;
    GeoElementSide neighbour = Neighbour();
    if (neighbour.Exists()) return false;
    while(neighbour != *this){
        if (candidate == neighbour) return true;
        neighbour = neighbour.Neighbour();
    } return false;
}

void GeoElementSide::IsertConnectivity(GeoElementSide &candidate) {
    auto Exists = [](GeoElementSide trial) -> bool{
        return (trial.fElement != 0 && trial.fSide >-1);
    };

    if(!Exists(*this) || !Exists(candidate)) {std::cout << "GeoElementSide::IsertConnectivity: GeoElementSide variable improperly initialized\n"; DebugStop();}

    GeoElementSide neighneigh = candidate.Neighbour();
    GeoElementSide curentneigh= Neighbour();

    if (!Exists(neighneigh) && !Exists(curentneigh)){
        SetNeighbour(candidate);
        candidate.SetNeighbour(*this);
    } else if (Exists(neighneigh) && Exists(curentneigh)){
        bool a = IsNeighbour(candidate);
        bool b = candidate.IsNeighbour(*this);
        if ((!a && b) || (a && !b)) { std::cout <<"Inconsitent data structure";DebugStop;}
        else if (!a) {// "neither *this is neighbour of candidate nor candidate is neighbour of this"
            SetNeighbour(neighneigh);
            candidate.SetNeighbour(curentneigh);
         }
    } else if (Exists(neighneigh)){
        SetNeighbour(neighneigh);
        candidate.SetNeighbour(*this);
    } else if (Exists(curentneigh)){
        SetNeighbour(candidate);
        candidate.SetNeighbour(curentneigh);
    }
}

void GeoElementSide::AllNeighbours(std::vector<GeoElementSide> &allneigh) {
    if (!Neighbour().Exists() || !Exists()) {
        std::cout << "GeoElementSide::AllNeighbours: Warning! Undefined GeoElementSide\n";DebugStop();
    } GeoElementSide neighbour = Neighbour();
    while(neighbour !=*this)  {
        allneigh.push_back(neighbour);
        neighbour = neighbour.Neighbour();
    }
}

void GeoElementSide::ComputeNeighbours(std::vector<GeoElementSide> &compneigh) {
    if(fSide < fElement->NCornerNodes()){ //if side is a corner...
        AllNeighbours(compneigh);
    } if (!Exists()) {"Uninitialized GeoElementSide";DebugStop();}

    int numberNodes = fElement->NSideNodes(fSide);
    std::vector<GeoElementSide> GeoElSideSet;
    VecInt GeoElSet[15];

    VecInt nodeIndexes(numberNodes);
    //Iterating through each corner of the side and filling GeoElSet[i], 0<i<numberNodes, with the neighbour of the i_th node.
    for (int i =0; i<numberNodes;i++){
        // Computing the neighbouring geoEl of every node...
        //nodeIndexes[i] = fElement->NodeIndex(i);
        int locNod = fElement->SideNodeIndex(fSide,i);
        nodeIndexes[i] = fElement->NodeIndex(locNod);
        GeoElSideSet.resize(0);
        GeoElementSide locSide(fElement,locNod);
        locSide.AllNeighbours(GeoElSideSet);
        // Building each GeoElSet layer with the Index of every neighbour GeoEl of the current corner;
        int nel = GeoElSideSet.size();
        int el;
        for(el=0; el<nel; el++) {
            GeoElSet[i].push_back(GeoElSideSet[el].Element()->GetIndex());
        }
        //Sorting GeoElSet in a crescent order. This step is a must for using the intersect function.
        std::sort(GeoElSet[i].begin(),GeoElSet[i].end());
    }
    // Computing the intersection between the neighbouring nodes...
    VecInt result;
    switch(numberNodes){
        case 1:
            result = GeoElSet[0];
        case 2:
            Intersect(GeoElSet[0],GeoElSet[1],result);
            break;
        case 3:
            Intersect(GeoElSet[0],GeoElSet[1],GeoElSet[2],result);
            break;
        case 4:
        {
            VecInt inter1, inter2;
            Intersect(GeoElSet[0],GeoElSet[2],inter1);
            if(inter1.size()) break;
            Intersect(GeoElSet[1],GeoElSet[3],inter2);
            if(inter2.size()) break;
            Intersect(inter1,inter2,result);
        }
        default: // For sides with more then 4 nodes...
        {
            VecInt inter1, inter2;
            inter1 = GeoElSet[0];
            for(int in=0; in<numberNodes-1; in++) {
                inter2.resize(0);
                Intersect(inter1,GeoElSet[in+1],inter2);
                if(inter2.size() == 0) break;
                inter1 = inter2;
            }
            result = inter2;
        }
    } //
    int el,nel = result.size();
    GeoMesh *geoMesh = fElement->GetMesh();
    for(el=0; el<nel; el++) {
        GeoElement *gelResult = geoMesh->Element(result[el]);
        int whichSd = gelResult->WhichSide(nodeIndexes);
        if(whichSd > 0)
        {
            compneigh.push_back(GeoElementSide( gelResult, whichSd));
        }
    }
}
bool GeoElementSide::DataConsistency(GeoElementSide &candidate) {
    DebugStop();
}
