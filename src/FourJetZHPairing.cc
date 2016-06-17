////////////////////////////////////////////////////////////////
// Class FourJetZHPairing
//
// Author: Jenny List
// Last update: $Date: 2016/06/15 10:19:08 $
//          by: $Author: boehmej $
// 
// Description: handle permutations of 4 jets into two different bosons
//               
////////////////////////////////////////////////////////////////

#include <iostream>
#include "FourJetZHPairing.h"
#include "JetFitObject.h"

FourJetZHPairing::FourJetZHPairing (JetFitObject *jets_[])  {
       
  for (int i = 0; i < NJETS; ++i) jets[i] = jets_[i];
  iperm = 0;
   
 
  int perms[NPERM][NJETS] = {{1, 2, 3, 4}, 
                             {1, 3, 2, 4}, 
                             {1, 4, 2, 3}, 
                             {3, 4, 1, 2}, 
                             {2, 4, 1, 3}, 
                             {2, 3, 1, 4}};
                             
  for (int i = 0; i < NPERM; ++i) 
     for (int j = 0; j < NJETS; ++j) 
                          permutations[i][j] = perms[i][j];    
                                                   
}

int FourJetZHPairing::nextPermutation (JetFitObject *permObjects[]) {
    
  for (int ijet = 0; ijet < NJETS; ++ijet) {
    permObjects[ijet] = jets[permutations[iperm][ijet]-1];
  } 
  
  ++iperm;
  return iperm;
}




