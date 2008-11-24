/*! \file 
 *  \brief Header file for some CERNLIB routines
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.1  2008/02/12 10:19:07  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.1  2007/10/30 15:51:14  gaede
 * - initial version of MarlinKinfit
 * -
 * - Revision 1.4  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 *
 */ 

#ifndef H_CERNLIB
#define H_CERNLIB
#include "ftypes.h"
extern "C" {
 // F010
  void reqn_  (FInteger *n, FReal   *a, FInteger *idim, FInteger *ir, FInteger *ifail,
               FInteger *k, FReal   *b);   
  void deqn_  (FInteger *n, FDouble *a, FInteger *idim, FInteger *ir, FInteger *ifail,
               FInteger *k, FDouble *b);   
  void rinv_  (FInteger *n, FReal   *a, FInteger *idim, FInteger *ir, FInteger *ifail);   
  void dinv_  (FInteger *n, FDouble *a, FInteger *idim, FInteger *ir, FInteger *ifail);   
  void reqinv_(FInteger *n, FReal   *a, FInteger *idim, FInteger *ir, FInteger *ifail,
               FInteger *k, FReal   *b);   
  void deqinv_(FInteger *n, FDouble *a, FInteger *idim, FInteger *ir, FInteger *ifail,
               FInteger *k, FDouble *b);   
 // F012
  void rsinv_  (FInteger *n, FReal   *a, FInteger *idim, FInteger *ifail);   
  void dsinv_  (FInteger *n, FDouble *a, FInteger *idim, FInteger *ifail);   
  void rseqn_  (FInteger *n, FReal   *a, FInteger *idim, FInteger *ifail,
               FInteger *k, FReal   *b);   
  void dseqn_  (FInteger *n, FDouble *a, FInteger *idim, FInteger *ifail,
                FInteger *k, FDouble *b);   
  void rsfact_ (FInteger *n, FReal   *a, FInteger *idim, FInteger *ifail,
                             FReal   *b, FInteger *jfail);   
  void dsfact_ (FInteger *n, FDouble *a, FInteger *idim, FInteger *ifail,
                             FDouble *b, FInteger *jfail);   
  void rsfeqn_ (FInteger *n, FReal   *a, FInteger *idim,
               FInteger *k, FReal   *b);   
  void dsfeqn_ (FInteger *n, FDouble *a, FInteger *idim,
                FInteger *k, FDouble *b);   
  void rsfinv_ (FInteger *n, FReal   *a, FInteger *idim);   
  void dsfinv_ (FInteger *n, FDouble *a, FInteger *idim);   
  // G100
  FReal prob_  (FReal   *x, FInteger *n); 
 
 
  void rannor_ (FReal *a, FReal *b);
  
  // V113
  void ranmar_ (FReal *vec, FInteger *len);
  void rmarin_ (FInteger *ijklin, FInteger *ntotin, FInteger *nto2in);
  void rmarut_ (FInteger *ijklin, FInteger *ntotin, FInteger *nto2in);
  void rmmar_  (FReal *vec, FInteger *len, FInteger *iseq);
  //void rmmaq_  (...);
  // V115
  void ranlux_ (FReal *rvec, FInteger *len);
  void rluxgo_ (FInteger *lux, FInteger *int_, FInteger *k1, FInteger *k2);
  void rluxat_ (FInteger *lux, FInteger *int_, FInteger *k1, FInteger *k2);
  void rluxin_ (FInteger *ivec);
  void rluxut_ (FInteger *ivec);
  // V120
  void rnorml_ (FReal *rvec, FInteger *len);
  void rnormx_ (FReal *rvec, FInteger *len, void urng (FReal *, FInteger *));
}

// F010
inline int  reqn   (FInteger n, FReal   a[], FInteger idim, FInteger ir[], FInteger& ifail,
                    FInteger k, FReal   b[])  {
  reqn_  (&n, a, &idim, ir, &ifail, &k, b);
  return ifail;
}  
inline int  reqn   (FInteger n, FReal   a[], FInteger idim, FInteger ir[],
                    FInteger k, FReal   b[])  {
  FInteger ifail = 0;
  reqn_  (&n, a, &idim, ir, &ifail, &k, b);
  return ifail;
}  
inline int  deqn   (FInteger n, FDouble a[], FInteger idim, FInteger ir[], FInteger& ifail,
                    FInteger k, FDouble b[]) {
  deqn_  (&n, a, &idim, ir, &ifail, &k, b);
  return ifail;
}    
inline int  deqn   (FInteger n, FDouble a[], FInteger idim, FInteger ir[], 
                    FInteger k, FDouble b[]) {
  FInteger ifail = 0;
  deqn_  (&n, a, &idim, ir, &ifail, &k, b);
  return ifail;
}   
 
inline int  rinv   (FInteger n, FReal   a[], FInteger idim, FInteger ir[], FInteger& ifail) {
  rinv_  (&n, a, &idim, ir, &ifail);
  return ifail;
}
inline int  rinv   (FInteger n, FReal   a[], FInteger idim, FInteger ir[]) {
  FInteger ifail = 0;
  rinv_  (&n, a, &idim, ir, &ifail);
  return ifail;
}
inline int  dinv   (FInteger n, FDouble a[], FInteger idim, FInteger ir[], FInteger& ifail) {
  dinv_  (&n, a, &idim, ir, &ifail);
  return ifail;
}   
inline int  dinv   (FInteger n, FDouble a[], FInteger idim, FInteger ir[]) {
  FInteger ifail = 0;
  dinv_  (&n, a, &idim, ir, &ifail);
  return ifail;
}

inline int  reqinv (FInteger n, FReal   a[], FInteger idim, FInteger ir[], FInteger& ifail,
                    FInteger k, FReal   b[])  {
  reqinv_(&n, a, &idim, ir, &ifail, &k, b);
  return ifail;
}  
inline int  reqinv (FInteger n, FReal   a[], FInteger idim, FInteger ir[],
                    FInteger k, FReal   b[])  {
  FInteger ifail = 0;
  reqinv_(&n, a, &idim, ir, &ifail, &k, b);
  return ifail;
}  
inline int  deqinv (FInteger n, FDouble a[], FInteger idim, FInteger ir[], FInteger& ifail,
                    FInteger k, FDouble b[]) {
  deqinv_(&n, a, &idim, ir, &ifail, &k, b);
  return ifail;
}    
inline int  deqinv (FInteger n, FDouble a[], FInteger idim, FInteger ir[], 
                    FInteger k, FDouble b[]) {
  FInteger ifail = 0;
  deqinv_(&n, a, &idim, ir, &ifail, &k, b);
  return ifail;
}   
 
// F012
inline int  rsinv  (FInteger n, FReal   a[], FInteger idim, FInteger& ifail) {
  rsinv_  (&n, a, &idim, &ifail);
  return ifail;
}
inline int  rsinv  (FInteger n, FReal   a[], FInteger idim) {
  FInteger ifail = 0;
  rsinv_  (&n, a, &idim, &ifail);
  return ifail;
}
inline int  dsinv  (FInteger n, FDouble a[], FInteger idim, FInteger& ifail) {
  dsinv_  (&n, a, &idim, &ifail);
  return ifail;
}   
inline int  dsinv  (FInteger n, FDouble a[], FInteger idim) {
  FInteger ifail = 0;
  dsinv_  (&n, a, &idim, &ifail);
  return ifail;
}

inline int  rseqn  (FInteger n, FReal   a[], FInteger idim, FInteger& ifail,
                    FInteger k, FReal   b[])  {
  rseqn_  (&n, a, &idim, &ifail, &k, b);
  return ifail;
}  
inline int  rseqn  (FInteger n, FReal   a[], FInteger idim,
                    FInteger k, FReal   b[])  {
  FInteger ifail = 0;
  rseqn_  (&n, a, &idim, &ifail, &k, b);
  return ifail;
}  
inline int  dseqn  (FInteger n, FDouble a[], FInteger idim, FInteger& ifail,
                    FInteger k, FDouble b[]) {
  dseqn_  (&n, a, &idim, &ifail, &k, b);
  return ifail;
}    
inline int  dseqn  (FInteger n, FDouble a[], FInteger idim, 
                    FInteger k, FDouble b[]) {
  FInteger ifail = 0;
  dseqn_  (&n, a, &idim, &ifail, &k, b);
  return ifail;
}   
 
inline int  rsfact (FInteger n, FReal   a[], FInteger idim, FInteger& ifail,
                                FReal   b[], FInteger& jfail) {
  rsfact_  (&n, a, &idim, &ifail, b, &jfail);
  return ifail;
}  
inline int  rsfact (FInteger n, FReal   a[], FInteger idim, 
                                FReal   b[], FInteger& jfail) {
  FInteger ifail = 0;
  rsfact_  (&n, a, &idim, &ifail, b, &jfail);
  return ifail;
}  
inline int  dsfact (FInteger n, FDouble a[], FInteger idim, FInteger& ifail,
                                FDouble b[], FInteger& jfail) {
  dsfact_  (&n, a, &idim, &ifail, b, &jfail);
  return ifail;
}    
inline int  dsfact (FInteger n, FDouble a[], FInteger idim, 
                                FDouble b[], FInteger& jfail) {
  FInteger ifail = 0;
  dsfact_  (&n, a, &idim, &ifail, b, &jfail);
  return ifail;
}    

inline void rsfeqn (FInteger n, FReal   a[], FInteger idim,
                    FInteger k, FReal   b[]) {
  rsfeqn_  (&n, a, &idim, &k, b);
}    
inline void dsfeqn (FInteger n, FDouble a[], FInteger idim,
                    FInteger k, FDouble b[]) {
  dsfeqn_  (&n, a, &idim, &k, b);
}       
inline void rsfinv (FInteger n, FReal   a[], FInteger idim) {
  rsfinv_  (&n, a, &idim);
}      
inline void dsfinv (FInteger n, FDouble a[], FInteger idim) {
  dsfinv_  (&n, a, &idim);
}
         
// G100
inline FReal prob  (FReal   x, FInteger n) {
  return prob_ (&x, &n);
}
  

inline void rannor (FReal& a, FReal& b) {
  rannor_ (&a, &b);
}
  // V113
inline void ranmar (FReal vec[], FInteger len) {
  ranmar_ (vec, &len);
}
inline void rmarin (FInteger ijklin, FInteger ntotin, FInteger nto2in) {
  rmarin_ (&ijklin, &ntotin, &nto2in);
}
inline void rmarut (FInteger& ijklin, FInteger& ntotin, FInteger& nto2in) {
  rmarut_ (&ijklin, &ntotin, &nto2in);
}
inline void rmmar  (FReal vec[], FInteger len, FInteger iseq) {
  rmmar_ (vec, &len, &iseq);
}
 
  // V115
inline void ranlux (FReal rvec[], FInteger len) {
  ranlux_ (rvec, &len); 
}
inline void rluxgo (FInteger lux, FInteger int_, FInteger k1, FInteger k2) {
  rluxat_ (&lux, &int_, &k1, &k2);
}
inline void rluxat (FInteger& lux, FInteger& int_, FInteger& k1, FInteger& k2) {
  rluxat_ (&lux, &int_, &k1, &k2);
}
inline void rluxin (FInteger ivec[25]) {
  rluxin_ (ivec);
}
inline void rluxut (FInteger ivec[25]) {
  rluxut_ (ivec);
}
  // V120
inline void rnorml (FReal rvec[], FInteger len) {
  rnorml_ (rvec, &len);
}

#endif
