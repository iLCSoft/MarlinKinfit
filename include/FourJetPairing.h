/*! \file 
 *  \brief Declares class FourJetPairing
 *
 * \b Changelog:
 *
 * \b CVS Log messages:
 * - $Log: FourJetPairing.h,v $
 * - Revision 1.1  2008/02/12 10:19:05  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.3  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 *
 */ 

#ifndef __FOURJETPAIRING_H
#define __FOURJETPAIRING_H

#include "BaseJetPairing.h"
#include "JetFitObject.h"

#include <iostream>
#include <array>

//  Class FourJetPairing:
/// Class to handle permutations of 2b jets and 4 light jets
/**
 *
 * Author: Jenny List, Anca Siebel
 * Last update: $Date: 2008/02/12 10:19:05 $
 *          by: $Author: blist $
 *
 */

class FourJetPairing : public BaseJetPairing {
  protected:
    constexpr static int NPERM = 3;
    constexpr static int NJETS = 4;

  public:
    /// constructor
    FourJetPairing(std::array<JetFitObject*, NJETS> jets_);
    FourJetPairing(const FourJetPairing&) = delete;
    FourJetPairing& operator=(const FourJetPairing&) = delete;
    FourJetPairing(FourJetPairing&&) = delete;
    FourJetPairing& operator=(FourJetPairing&&) = delete;
 
    /// Virtual destructor
    virtual ~FourJetPairing() = default;
        
    /// Number of permutaions
    virtual int getNPerm() const {return NPERM;};
    
    /// does the job
    virtual int nextPermutation (JetFitObject *permObjects[]);
    
  protected:
    std::array<JetFitObject*, NJETS> jets;
    int permutations [NPERM][NJETS];

};
    
#endif // __FOURJETPAIRING_H

