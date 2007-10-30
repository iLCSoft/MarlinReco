////////////////////////////////////////////////////////////////
// Class FourJetPairing
//
// Author: Jenny Boehme, Anca Siebel
// Last update: $Date: 2007-10-30 15:51:14 $
//          by: $Author: gaede $
// 
// Description: handle permutations of 4 jets
//               
////////////////////////////////////////////////////////////////

#include <iostream>
#include "FourJetPairing.h"
#include "JetFitObject.h"

FourJetPairing::FourJetPairing (JetFitObject *jets_[])  {
       
  for (int i = 0; i < NJETS; ++i) jets[i] = jets_[i];
  iperm = 0;
   
  // this assumes jet 5 & 6 to be the b-jets !!
 
  int perms[NPERM][NJETS] = {{1, 2, 3, 4}, 
                             {1, 3, 2, 4}, 
                             {1, 4, 2, 3}};
                             
  for (int i = 0; i < NPERM; ++i) 
     for (int j = 0; j < NJETS; ++j) 
                          permutations[i][j] = perms[i][j];    
                                                   
}

int FourJetPairing::nextPermutation (JetFitObject *permObjects[]) {
    
  for (int ijet = 0; ijet < NJETS; ++ijet) {
    permObjects[ijet] = jets[permutations[iperm][ijet]-1];
  } 
  
  ++iperm;
  return iperm;
}




