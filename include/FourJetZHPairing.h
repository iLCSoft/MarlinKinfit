#ifndef __FOURJETZHPAIRING_H
#define __FOURJETZHPAIRING_H

#include <iostream>
#include "BaseJetPairing.h"
#include "JetFitObject.h"

//  Class FourJetZHPairing:
/// Class to handle permutations of 4 jets into two different bosons
/**
 *
 * Author: Jenny List
 * Last update: $Date: 2016/06/15 10:19:05 $
 *          by: $Author: boehmej $
 *
 */

class FourJetZHPairing : public BaseJetPairing {
  public:
    /// constructor
    FourJetZHPairing (JetFitObject *jets_[]);
    
    /// Virtual destructor
    virtual ~FourJetZHPairing() {};    
        
    /// Number of permutaions
    virtual int getNPerm() const {return NPERM;};
    
    /// does the job
    virtual int nextPermutation (JetFitObject *permObjects[]);
    
  protected:
    enum {NPERM = 6};
    enum {NJETS = 4};
    JetFitObject *jets[NJETS]; 
    int permutations [NPERM][NJETS];

};
    
#endif // __FOURJETZHPAIRING_H

