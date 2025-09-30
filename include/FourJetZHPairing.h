#ifndef __FOURJETZHPAIRING_H
#define __FOURJETZHPAIRING_H

#include "BaseJetPairing.h"
#include "JetFitObject.h"

#include <array>

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
  protected:
    constexpr static int NPERM = 6;
    constexpr static int NJETS = 4;

  public:
    /// constructor
    FourJetZHPairing (std::array<JetFitObject*, NJETS> jets_);
    FourJetZHPairing(JetFitObject* jets_[]) : FourJetZHPairing({jets_[0], jets_[1], jets_[2], jets_[3]}) {}

    FourJetZHPairing(const FourJetZHPairing&) = delete;
    FourJetZHPairing& operator=(const FourJetZHPairing&) = delete;
    FourJetZHPairing(FourJetZHPairing&&) = delete;
    FourJetZHPairing& operator=(FourJetZHPairing&&) = delete;

    /// Virtual destructor
    virtual ~FourJetZHPairing() = default;
        
    /// Number of permutaions
    virtual int getNPerm() const {return NPERM;};
    
    /// does the job
    virtual int nextPermutation (JetFitObject *permObjects[]);
    
  protected:
    std::array<JetFitObject*, NJETS> jets;
    int permutations [NPERM][NJETS];

};
    
#endif // __FOURJETZHPAIRING_H

