////////////////////////////////////////////////////////////////
// Class TwoB4JPairing
//
// Author: Jenny Boehme, Anca Siebel
// Last update: $Date: 2008/02/12 10:19:07 $
//          by: $Author: blist $
// 
// Description: handle permutations of 2b jets and 4 light jets
//               
////////////////////////////////////////////////////////////////

#ifndef __TWOB4JPAIRING_H
#define __TWOB4JPAIRING_H

#include "BaseJetPairing.h"
#include "JetFitObject.h"

#include<array>

class TwoB4JPairing : public BaseJetPairing {
  private:
    constexpr static int NPERM = 6;
    constexpr static int NJETS = 6;

  public:
    // constructor
    TwoB4JPairing (std::array<JetFitObject*, NJETS> jets_);
    
    // destructor
    virtual ~TwoB4JPairing() = default;
        
    // getters
    virtual int getNPerm() const {return NPERM;};
    
    // does the job
    virtual int nextPermutation (JetFitObject *permObjects[]);
    
  protected:
    std::array<JetFitObject*, NJETS> jets;
    int permutations [NPERM][NJETS];

};
    
#endif // __TWOB4JPAIRING_H

