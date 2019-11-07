//
//  GeoElementSide.h
//  FEMCoursePZ
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoElement.h"

#ifndef GeoElementSide_h
#define GeoElementSide_h

class GeoElement;

class GeoElementSide
{
    // Associated element
    GeoElement *fElement;
    
    // Associated side
    int fSide;
    
public:
    
    // Default Constructor of GeoElementSide
    GeoElementSide();
    
    // Constructor of GeoElementSide
    GeoElementSide(GeoElement *element, int side) : fElement(element), fSide(side)
    {
    }
    
    // Copy constructor of GeoElementSide
    GeoElementSide(const GeoElementSide &copy);
    
    // Operator of copy 
    GeoElementSide &operator=(const GeoElementSide &copy);
    
    int operator==(const GeoElementSide &other) const {
        return fElement == other.fElement && fSide == other.fSide;
    }
    int operator!=(const GeoElementSide &other) const {
        return ! operator==(other);
    }
    
    // Return the associated element
    GeoElement *Element() const
    {
        return fElement;
    }

    // Return the associated side
    int Side() const
    {
        return fSide;
    }

    // Return neighbour element of a given side
    GeoElementSide Neighbour() const;
    
    bool DataConsistency(GeoElementSide &candidate);
    
    int Exists() const {return (fElement != 0 && fSide > -1);}
    
    // Fill in the data structure for the neighbouring information
    void SetNeighbour(const GeoElementSide &neighbour);
    
    // Verifiy if an element is a neighbour
    bool IsNeighbour(const GeoElementSide &candidate);
    
    // Define elements neighbourhood
    void IsertConnectivity(GeoElementSide &connectivity);
    
    // Vector with all Neighbours
    void AllNeighbours(std::vector<GeoElementSide> &allneigh);
    
    // Compute all corner neighbours
    void ComputeNeighbours(std::vector<GeoElementSide> &neighbour);

    void Intersect(const std::vector< int > &one, const std::vector< int > &two, std::vector< int > &result) {
        int firstc, secondc, nfirst, nsecond;
        nfirst = one.size();
        nsecond = two.size();
        firstc = 0;
        secondc = 0;
        while (firstc < nfirst && secondc < nsecond) {
            while (firstc < nfirst && one[firstc] < two[secondc]) {
                firstc++;
            }
            if (firstc == nfirst) break;
            while (secondc < nsecond && two[secondc] < one[firstc]) {
                secondc++;
            }
            if (firstc < nfirst && secondc < nsecond && one[firstc] == two[secondc]) {
                result.push_back(one[firstc]);
                firstc++;
                secondc++;
            }
        }
    }

    void Intersect(const std::vector< int > &one, const std::vector< int > &two, const std::vector< int > &three, std::vector< int > &result) {
        int firstc, secondc, thirdc, nfirst, nsecond, nthird;
        nfirst = one.size();
        nsecond = two.size();
        nthird = three.size();
        firstc = 0;
        secondc = 0;
        thirdc = 0;
        while (firstc < nfirst && secondc < nsecond && thirdc < nthird) {
            while (firstc < nfirst && (one[firstc] < two[secondc] || one[firstc] < three[thirdc])) {
                firstc++;
            }
            if (firstc == nfirst)break;
            while (secondc < nsecond && (two[secondc] < one[firstc] || two[secondc] < three[thirdc])) {
                secondc++;
            }
            if (secondc == nsecond) break;
            while (thirdc < nthird && (three[thirdc] < one[firstc] || three[thirdc] < two[secondc])) {
                thirdc++;
            }
            if (firstc < nfirst && secondc < nsecond && thirdc < nthird && one[firstc] == two[secondc] && one[firstc] == three[thirdc]) {
                result.push_back(one[firstc]);
                firstc++;
                secondc++;
                thirdc++;
            }
        }
    }
    
    
};
#endif /* GeoElementSide_h */
