////////////////////////////////////////////////////////////////
// Class FourVector
//
// Author: Benno List
// Last update: $Date: 2008-02-25 08:23:56 $
//          by: $Author: gaede $
// 
// Description: class for four-vectors
//               
////////////////////////////////////////////////////////////////

#include "FourVector.h"

#include <cmath>
#include <cassert>
#include "ftypes.h"
#include "cernlib.h"


FourVector& FourVector::boost (const FourVector& P) {
  // See CERNLIB U101 for a description
  
  double pP = -(p*P.p);
  double e = getE();
  double M = P.getM();
  
  E = (e*P.getE() - pP)/M;
  p = p - ((pP/(P.getE()+M)-e)/M)*P.p;
  
  return *this;
}

void FourVector::decayto (FourVector& d1, FourVector& d2) const {
  // Let this particle decay isotropically into 4-vectors d1 and d2;
  // d1 and d2 must have definite mass at beginning
  using std::abs;
  using std::sqrt; 
  using std::pow;

  double M2 = getM2();
  double M  = getM();;
  double m1 = d1.getM();
  double m2 = d2.getM();
  
  FReal randoms[2];
//  FInteger ilen = 2;
//  ranmar_ (randoms, &ilen);
  ranmar (randoms, 2);
  
  assert (m1+m2<=M);
  
  double pstar = 0.5*sqrt (abs((M2-pow(m1+m2,2))*(M2-pow(m1-m2,2))))/M;
  double phistar = 2*M_PI*randoms[0];
  double costhetastar = 2*randoms[1]-1;
  double sinthetastar = sqrt(abs (1-costhetastar*costhetastar));
  double E1 = sqrt(m1*m1+pstar*pstar);
  double E2 = sqrt(m2*m2+pstar*pstar);
  
//  cout << "pstar=" << pstar << ", E1=" << E1 << ", E2=" << E2 << endl;
  
  
          
  d1 = FourVector (E1, pstar*sinthetastar*cos(phistar), 
                       pstar*sinthetastar*sin(phistar),  
                       pstar*costhetastar);
  d2 = FourVector (E2, -pstar*sinthetastar*cos(phistar),  
                       -pstar*sinthetastar*sin(phistar),  
                       -pstar*costhetastar);
                       
//   cout << "d1 = " << d1 << "\nd2 = " << d2 << "\nsum= " << d1+d2 << ", mass: " << (d1+d2).getM() << endl;
  d1.boost (*this);
  d2.boost (*this);
  
//   std::cout << "Decay of " << mother 
//             << "\nto  " << d1 
//             << "\nand " << d2 
//             << "\n:   " << mother-(d1+d2) << endl;
}
