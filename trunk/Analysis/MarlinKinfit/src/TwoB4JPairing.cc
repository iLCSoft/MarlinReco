////////////////////////////////////////////////////////////////
// Class TwoB4JPairing
//
// Author: Jenny Boehme, Anca Siebel
// Last update: $Date: 2008/02/12 10:19:10 $
//          by: $Author: blist $
// 
// Description: handle permutations of 2b jets and 4 light jets
//               
////////////////////////////////////////////////////////////////

#include "TwoB4JPairing.h"

#include <iostream>
#include "JetFitObject.h"

TwoB4JPairing::TwoB4JPairing (JetFitObject *jets_[])  {
       
  for (int i = 0; i < NJETS; ++i) jets[i] = jets_[i];
  iperm = 0;
   
  // this assumes jet 5 & 6 to be the b-jets !!
 
  int perms[NPERM][NJETS] = {{1, 2, 3, 4, 5, 6}, 
                             {1, 3, 2, 4, 5, 6}, 
                             {1, 4, 2, 3, 5, 6}, 
                             {1, 2, 3, 4, 6, 5}, 
                             {1, 3, 2, 4, 6, 5}, 
                             {1, 4, 2, 3, 6, 5}};
                             
  for (int i = 0; i < NPERM; ++i) 
     for (int j = 0; j < NJETS; ++j) 
                          permutations[i][j] = perms[i][j];    
                                                   
}

int TwoB4JPairing::nextPermutation (JetFitObject *permObjects[]) {
    
  for (int ijet = 0; ijet < NJETS; ++ijet) {
    permObjects[ijet] = jets[permutations[iperm][ijet]-1];
  } 
  
  ++iperm;
  return iperm;
}




