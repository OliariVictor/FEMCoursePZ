/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoNode.h"
#include "tpanic.h"

void GeoNode::Print(std::ostream &out) {
    int j = xco.size();
    out << "Node coordinate:\t";
    for (int i = 0; i < j; i++ ) std::cout  << xco[i] <<"\t";
    out << std::endl;
}
