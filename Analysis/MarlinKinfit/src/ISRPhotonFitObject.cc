/*! \file 
 *  \brief Implements class ISRPhotonFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: ISRPhotonFitObject.cc,v $
 * - Revision 1.4  2011/03/16 16:33:24  mbeckman
 * - Compatibility fixes with ILCSoft svn
 * -
 * - Revision 1.3  2011/03/03 15:32:32  boehmej
 * - removed obsolete PhotonFitObject.o and PhotonFitObjectPgxy.o from src/Makefile, activated NO_MARLIN in src/ISRPhotonFitObject.cc
 * -
 * - Revision 1.2  2010/07/05 20:08:43  mbeckman
 * - ISRPhotonFitObject.cc: Included flag for output via cout/marlin
 * -
 * - Revision 1.1  2010/06/11 20:32:51  mbeckman
 * - Renamed PhotonFitObjects, cleaned them up for usage
 * -
 * - Revision 1.10  2009/04/02 12:47:35  mbeckman
 * - PhotonFitObject.cc, PseudoMeasuredPhotonFitObjectPxyz.cc: bug fix (measured p = 0 instead of start value)
 * - PhotonFitObjectPxyg.cc: added assertion to catch up division by zero
 * -
 * - Revision 1.9  2009/04/01 09:00:15  mbeckman
 * - Corrected derivatives by factor sqrt(2)
 * -
 * - Revision 1.8  2009/03/26 08:47:30  mbeckman
 * - Bug fix (measured p = 0 instead of start value)
 * -
 * - Revision 1.7  2009/03/17 13:27:03  mbeckman
 * - fixed for compiling with gcc 3.2.3 (for compability with ILCSoft)
 * -
 * - Revision 1.6  2009/02/23 12:04:05  mbeckman
 * - - PhotonFitObject:     bug fix (1/0), removed dispensable variables
 * - - PhotonFitObjectPxyg: bug fixes (1/0, order of computing variables), modified parametrization
 * - - JetFitObject:        added start parameter check (inf, nan)
 * -
 * - Revision 1.5  2009/02/18 11:56:22  mbeckman
 * - PhotonFitObject*.cc: documentation, debug output
 * - NewtonFitterGSL.cc:  bug fix (Lagrange multipliers not initialized), debug output
 * - JetFitObject.cc:     bug fix: division by 0, if energy <= mass
 *
 */ 

#define NO_MARLIN		// if defined: all output via cout, Marlin inclusion not required
#include "ISRPhotonFitObject.h"
#include <cmath>
#include <cassert>
#include <iostream>
#ifndef NO_MARLIN
#include "marlin/Processor.h"
#endif

using std::sqrt;
using std::exp;
using std::pow;
using std::cout; 
using std::endl;
#ifndef NO_MARLIN
using namespace marlin;
#endif

static const double pi_ = 3.14159265358979323846264338328,
                    a   = 8./3./pi_*(pi_-3.)/(4.-pi_);    // = ca. 0.140012289

// constructor
ISRPhotonFitObject::ISRPhotonFitObject(double px, double py, double pz,
                                         double b_, double PzMaxB_, double PzMinB_) {
  initCov();
  b = b_;
  PzMinB = PzMinB_;
  PzMaxB = PzMaxB_;
  #ifdef DEBUG
    cout << "ISRPhotonFitObject:   b: " << b << "   PzMinB: " << PzMinB << "   PzMaxB: " << PzMaxB << endl;
  #endif

  if(b <= 0. || b >= 1.){
    cout << "ISRPhotonFitObject:   b must be from ]0,1[ "  << endl;
  }
  assert(b > 0. && b < 1.);
  if(PzMinB < 0. || PzMaxB <= PzMinB){
    cout << "ISRPhotonFitObject:   PzMinB and PzMaxB must be chosen such that 0 <= PzMinB < PzMaxB"  << endl;
  }
  assert(PzMinB >= 0.);
  assert(PzMaxB > PzMinB);
  dp2zFact = (PzMaxB-PzMinB)/b*sqrt(2./pi_);
  double pg = PgFromPz(pz);         // using internally Gauss-distributed parameter p_g instead of p_z
  setParam (0, px, true, true);
  setParam (1, py, true, true);
  setParam (2, pg, true);
  setMParam (0, 0.);                // all measured parameters
  setMParam (1, 0.);                // are assumed to be zero
  setMParam (2, 0.);                // in this photon parametrization
  #ifdef DEBUG
    cout << "ISRPhotonFitObject:   Initial pg: " << pg << endl;
  #endif
  setError (2, 1.);
  setMass (0.);
  invalidateCache();
}

// destructor
ISRPhotonFitObject::~ISRPhotonFitObject() {}


ISRPhotonFitObject *ISRPhotonFitObject::copy() const {
  return new ISRPhotonFitObject (*this);
}
    
ISRPhotonFitObject& ISRPhotonFitObject::assign (const BaseFitObject& source) {
  if (const ISRPhotonFitObject *psource = dynamic_cast<const ISRPhotonFitObject *>(&source)) {
    if (psource != this) *this = *psource;
  }
  else {
    assert (0);
  }
  return *this;
}


const char *ISRPhotonFitObject::getParamName (int ilocal) const {
  switch (ilocal) {
    case 0: return "P_x";
    case 1: return "P_y";
    case 2: return "P_g";
  }
  return "undefined";
}

// needed for constructor!
bool ISRPhotonFitObject::setParam (int ilocal, double par_, 
                                    bool measured_, bool fixed_) {
  assert (ilocal >= 0 && ilocal < 3);
  if (measured[ilocal] != measured_ || fixed[ilocal] != fixed_) invalidateCache();
  measured[ilocal] = measured_;
  fixed[ilocal] = fixed_;
// this doesn't work in constructor, since par[i], mass etc are not initialized yet!!!!  
//  return setParam (ilocal, par_);
// old version of bool ParticleFitObject::setParam (int ilocal, double par_ )
//  if (!isfinite(par_)) return false;  doesn't exist anymore?
//  assert (ilocal >= 0 && ilocal < NPAR);   done before
//  if (par[ilocal] == par_) return true;    
//  invalidateCache();                     done in constructor anyhow
  par[ilocal] = par_;
  return true;
} 

bool ISRPhotonFitObject::setParam (int i, double par_ ) {
  invalidateCache();
  bool result = (par_-par[i])*(par_-par[i]) > eps2*cov[i][i]; 
  switch (i) {
    // p_x
    case 0: par[0] = par_;
            break;
    // p_y
    case 1: par[1] = par_;
            break;          
    // p_g
    case 2: par[2] = par_;
            break;          
    default: std::cerr << "ISRPhotonFitObject::setParam:   Illegal i=" << i << std::endl;
  }
  return result;
}
 
bool ISRPhotonFitObject::updateParams (double p[], int idim) {
  invalidateCache();
  
  int i2 = getGlobalParNum(2);
  assert (i2 >= 0 && i2 < idim);
  
  double p2 = p[i2];
  #ifdef DEBUG
    std::cout << "ISRPhotonFitObject::updateParams:   p2(new) = " << p[i2] << "   par[2](old) = " << par[2] << endl;
  #endif
  bool result = ((p2-par[2])*(p2-par[2]) > eps2*cov[2][2]);

  par[2] = p2;
  p[i2] = par[2];
  return result;
}  

// these depend on actual parametrisation!
double ISRPhotonFitObject::getPx() const {
  return par[0];
}

double ISRPhotonFitObject::getPy() const {
  return par[1];
}

double ISRPhotonFitObject::getPz() const {
  if (!cachevalid) updateCache();
  return pz;
}

double ISRPhotonFitObject::getE() const {
  if (!cachevalid) updateCache();
  return p;
}
double ISRPhotonFitObject::getP() const {
  if (!cachevalid) updateCache();
  return p;
}
double ISRPhotonFitObject::getP2() const {
  if (!cachevalid) updateCache();
  return p2;
}
double ISRPhotonFitObject::getPt() const {
  if (!cachevalid) updateCache();
  return std::sqrt(pt2);
}
double ISRPhotonFitObject::getPt2() const {
  if (!cachevalid) updateCache();
  return pt2;
}

double ISRPhotonFitObject::getDPx(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpx0;
    case 1: return dpx1;
    case 2: return dpx2;
  }
  return 0; 
}

double ISRPhotonFitObject::getDPy(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpy0;
    case 1: return dpy1;
    case 2: return dpy2;
  }
  return 0; 
}

double ISRPhotonFitObject::getDPz(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dpz0;
    case 1: return dpz1;
    case 2: return dpz2;
  }
  return 0; 
}

double ISRPhotonFitObject::getDE(int ilocal) const {
  assert (ilocal >= 0 && ilocal < NPAR);
  if (!cachevalid) updateCache();
  switch (ilocal) {
    case 0: return dE0;
    case 1: return dE1;
    case 2: return dE2;
  }
  return 0; 
}
 
void   ISRPhotonFitObject::addToDerivatives (double der[], int idim, 
                                       double efact, double pxfact, 
                                       double pyfact, double pzfact) const {
  int i0 = globalParNum[0];
  int i1 = globalParNum[1];
  int i2 = globalParNum[2];
  assert (i0 >= 0 && i0 < idim);
  assert (i1 >= 0 && i1 < idim);
  assert (i2 >= 0 && i2 < idim);

  if (!cachevalid) updateCache();
  // for numerical accuracy, add up derivatives first,
  // then add them to global vector
  double der0 = efact*dE0;
  double der1 = efact*dE1;
  double der2 = efact*dE2;
  
  if (pxfact != 0) {
    der0 += pxfact*dpx0;
    der1 += pxfact*dpx1;
    der2 += pxfact*dpx2;
  }
  if (pyfact != 0) {
    der0 += pyfact*dpy0;
    der1 += pyfact*dpy1;
    der2 += pyfact*dpy2;
  }
  if (pzfact != 0) {
    der0 += pzfact*dpz0;
    der1 += pzfact*dpz1;
    der2 += pzfact*dpz2;
  }

  der[i0] += der0;
  der[i1] += der1;
  der[i2] += der2;
}
    
void   ISRPhotonFitObject::addTo2ndDerivatives (double der2[], int idim, 
                                          double efact, double pxfact, 
                                          double pyfact, double pzfact) const {
  int i2 = globalParNum[2];
  assert (i2 >= 0 && i2 < idim);
  
  if (!cachevalid) updateCache();
  
  double der22 = pzfact*d2pz22 + efact*d2E22;
    
  der2[idim*i2+i2] += der22;
}
    
void   ISRPhotonFitObject::addTo2ndDerivatives (double M[], int idim,  double lambda, double der[]) const {
  addTo2ndDerivatives (M, idim, lambda*der[0], lambda*der[1], lambda*der[2], lambda*der[3]);
}

void   ISRPhotonFitObject::addTo1stDerivatives (double M[], int idim, double der[], int kglobal) const {
  assert (kglobal >= 0 && kglobal < idim);
  int i2 = globalParNum[2];
  assert (i2 >= 0 && i2 < idim);
  
  if (!cachevalid) updateCache();
  
  double d2 = der[0]*dE2 + der[1]*dpx2 + der[2]*dpy2 + der[3]*dpz2;
  
  M[idim*kglobal + i2] += d2; 
  M[idim*i2 + kglobal] += d2; 
}

         
double ISRPhotonFitObject::PgFromPz(double pz){
  int sign = (pz>0.) - (pz<0.);
  double u = ( pow(fabs(pz),b) - PzMinB ) / (PzMaxB-PzMinB);

//   if(u<0.){
//   #ifdef NO_MARLIN
//     cout << 
//   #else
//     m_out(WARNING) << 
//   #endif
//     "ISRPhotonFitObject: Initial pz with abs(pz) < pzMin adjusted to zero." << std::endl;
//     u = 0.;
//   }
// 
//   if(u>=1.){
//   #ifdef NO_MARLIN
//     cout << 
//   #else
//     m_out(WARNING) << 
//   #endif
//     "ISRPhotonFitObject: Initial pz with abs(pz) >= pzMax adjusted." << std::endl;
//     u = 0.99999999;
//   }

  double g = std::log(1.-u*u);
  double g4pa = g + 4./pi_/a;
  return sign*sqrt( -g4pa+sqrt( g4pa*g4pa-4./a*g ) ) ;
}

void ISRPhotonFitObject::initCov() {
  for (int i = 0; i < NPAR; ++i) {
    for (int j = 0; j < NPAR; ++j) {
      cov[i][j] = static_cast<double>(i == j);
    }
  }    
}

void ISRPhotonFitObject::invalidateCache() const {
  cachevalid = false;
}

void ISRPhotonFitObject::addToGlobalChi2DerVector (double *y, int idim, 
                                             double lambda, double der[]) const {
  int i2 = globalParNum[2];
  assert (i2 >= 0 && i2 < idim);
  
  if (!cachevalid) updateCache();
  
  double d2 = der[0]*dE2 + der[1]*dpx2 + der[2]*dpy2 + der[3]*dpz2;

  y[i2] += lambda*d2;
}

void ISRPhotonFitObject::updateCache() const {
  double px = par[0];
  double py = par[1];
  double pg = par[2];

  int sign = (pg>0.) - (pg<0.);
  double pg2h = pg*pg/2.;
  double exponent = -pg2h*(4./pi_+a*pg2h)/(1.+a*pg2h);
  double u = sqrt( (exponent<-1.e-14) ? 1.-exp( exponent ) : -exponent );  // approximation to avoid numerical problem
  pz = sign*pow( ( PzMinB + (PzMaxB-PzMinB)*u ) , (1./b) );

  pt2 = px*px+py*py;
  p2  = pt2+pz*pz;
  p   = std::sqrt(p2);
  
  dpx0 = 1.;
  dpx1 = 0.;
  dpx2 = 0.;
  dpy0 = 0.;
  dpy1 = 1.;
  dpy2 = 0.;
  dpz0 = 0.;
  dpz1 = 0.;
  dpz2 = dp2zFact*pow(fabs(pz),(1.-b))*exp(-pg*pg/2.);
  
  // if p,pz==0, derivatives are zero (catch up 1/0)
  if(pz){
    d2pz22 = dpz2*( (1.-b)*dpz2/pz - par[2] );
  }
  else{
    d2pz22 = 0.;
  }
  if(p){
    dE0   = px/p;
    dE1   = py/p;
    dE2   = pz/p*dpz2;
    d2E22 = pz/p*d2pz22;
    if(pt2){
      d2E22 += pt2/p/p/p*dpz2*dpz2;   // NOT using /p2/p to avoid numerical problems
    }
  }
  else{
    dE0   = 0.;
    dE1   = 0.;
    dE2   = 0.;
    d2E22 = 0.;
  }

  #ifdef DEBUG
    cout << "ISRPhotonFitObject::updateCache:   pg: " << pg << "   pz: " << pz << "   p: " << p << "   p^2: " << p2 << "\n"
         << "                                    Dpz/Dpg: " << dpz2 << "   DE/Dpg: " << dE2 << "   D^2pz/Dpg^2: " << d2pz22
                                          << "   D^2E/Dpg^2: " << d2E22 << endl;
  #endif
  
  cachevalid = true;
}

double ISRPhotonFitObject::getError2 (double der[]) const {
  if (!cachevalid) updateCache();
  double cov4[4][4]; // covariance of E, px, py, px
  cov4[0][0] =                                       // E, E
        dE0*(dE0*cov[0][0]+2*dE1*cov[0][1]+2*dE2*cov[0][2])+dE1*(dE1*cov[1][1]+2*dE2*cov[1][2])+dE2*dE2*cov[2][2];
  cov4[0][1] = cov4[1][0] =                          // E, px 
        dE0*dpx0*cov[0][0]+(dE0*dpx1+dE1*dpx0)*cov[0][1]+(dE0*dpx2+dE2*dpx0)*cov[0][2]
        +dE1*dpx1*cov[1][1]+(dE1*dpx2+dE2*dpx1)*cov[1][2]+dE2*dpx2*cov[2][2];
  cov4[0][2] = cov4[2][0] =                          // E, py
        dE0*dpy0*cov[0][0]+(dE0*dpy1+dE1*dpy0)*cov[0][1]+(dE0*dpy2+dE2*dpy0)*cov[0][2]
        +dE1*dpy1*cov[1][1]+(dE1*dpy2+dE2*dpy1)*cov[1][2]+dE2*dpy2*cov[2][2];
  cov4[0][3] = cov4[3][0] =                          // E, pz
        dE0*dpz0*cov[0][0]+(dE0*dpz1+dE1*dpz0)*cov[0][1]+(dE0*dpz2+dE2*dpz0)*cov[0][2]
        +dE1*dpz1*cov[1][1]+(dE1*dpz2+dE2*dpz1)*cov[1][2]+dE2*dpz2*cov[2][2];
  cov4[1][1] =                                       // px, px
        dpx0*(dpx0*cov[0][0]+2*dpx1*cov[0][1]+2*dpx2*cov[0][2])+dpx1*(dpx1*cov[1][1]+2*dpx2*cov[1][2])+dpx2*dpx2*cov[2][2];
  cov4[1][2] = cov4[2][1] =                          // px, py
        dpx0*dpy0*cov[0][0]+(dpx0*dpy1+dpx1*dpy0)*cov[0][1]+(dpx0*dpy2+dpx2*dpy0)*cov[0][2]
        +dpx1*dpy1*cov[1][1]+(dpx1*dpy2+dpx2*dpy1)*cov[1][2]+dpx2*dpy2*cov[2][2];
  cov4[1][3] = cov4[3][1] =                          // px, pz
        dpx0*dpz0*cov[0][0]+(dpx0*dpz1+dpx1*dpz0)*cov[0][1]+(dpx0*dpz2+dpx2*dpz0)*cov[0][2]
        +dpx1*dpz1*cov[1][1]+(dpx1*dpz2+dpx2*dpz1)*cov[1][2]+dpx2*dpz2*cov[2][2];
  cov4[2][2] =                                       // py, py
        dpy0*(dpy0*cov[0][0]+2*dpy1*cov[0][1]+2*dpy2*cov[0][2])+dpy1*(dpy1*cov[1][1]+2*dpy2*cov[1][2])+dpy2*dpy2*cov[2][2];
  cov4[2][3] = cov4[3][2] =                          // py, pz
        dpy0*dpz0*cov[0][0]+(dpy0*dpz1+dpy1*dpz0)*cov[0][1]+(dpy0*dpz2+dpy2*dpz0)*cov[0][2]
        +dpy1*dpz1*cov[1][1]+(dpy1*dpz2+dpy2*dpz1)*cov[1][2]+dpy2*dpz2*cov[2][2];
  cov4[3][3] =                                       // pz, pz
        dpz0*(dpz0*cov[0][0]+2*dpz1*cov[0][1]+2*dpz2*cov[0][2])+dpz1*(dpz1*cov[1][1]+2*dpz2*cov[1][2])+dpz2*dpz2*cov[2][2];
  return der[0]*(der[0]*cov4[0][0] + 2*der[1]*cov4[0][1] + 2*der[2]*cov4[0][2] + 2*der[3]*cov4[0][3])
                             + der[1]*(der[1]*cov4[1][1] + 2*der[2]*cov4[1][2] + 2*der[3]*cov4[1][3])
                                                   + der[2]*(der[2]*cov4[2][2] + 2*der[3]*cov4[2][3])
                                                                          + der[3]*der[3]*cov4[3][3];
}
