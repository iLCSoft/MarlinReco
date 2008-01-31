#ifndef __NEUTRINOFITOBJECT_H
#define __NEUTRINOFITOBJECT_H

#include"ParticleFitObject.h"

// Class NeutrinoFitObject
/// Class for neutrinos with (E, eta, phi) in kinematic fits
/**
 *
 * Author: Jenny List, Benno List
 * $Date: 2008-01-31 13:01:55 $
 * $Author: listj $
 *
 * \b Changelog:
 * - 30.12.04 BL: addToGlobCov, getDChi2DParam, getDChi2DParam2, addToGlobalChi2DerMatrix
 *            moved up to ParticleFitObject,
 *            getParamName implemented
 */ 
class NeutrinoFitObject : public ParticleFitObject {
  public:
    NeutrinoFitObject(double E, double theta, double phi, 
                 double DE, double Dtheta, double Dphi, 
                 double m = 0);
    virtual ~NeutrinoFitObject();
    
    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                     ) const;
 
    
    
    // these depend on actual parametrisation!
    virtual double getPx() const;
    virtual double getPy() const;
    virtual double getPz() const;
    virtual double getE() const;
    
    virtual double getP() const;
    virtual double getP2() const;  
    virtual double getPt() const;
    virtual double getPt2() const;
     
    virtual double getDPx(int ilocal) const;
    virtual double getDPy(int ilocal) const;
    virtual double getDPz(int ilocal) const;
    virtual double getDE(int ilocal) const;
    
    virtual void   addToDerivatives (double der[],
                                     int idim,
                                     double pxfact=0, 
                                     double pyfact=0, 
                                     double pzfact=0, 
                                     double efact=0) const;

    // add to matrix der2 of size idim x idim
    // pxfact*d^2px/(dx_i dx_j) + pyfact...
    virtual void   addTo2ndDerivatives (double der2[],
                                        int idim,
                                        double pxfact, 
                                        double pyfact, 
                                        double pzfact, 
                                        double efact) const;
    
    virtual void addToGlobalDerMatrix (int idim, double c, double *M) const;

    virtual void invalidateCache() const;
  
  protected:
    
    virtual void initCov();
    
    void updateCache() const;
  
    mutable bool cachevalid;
    
    mutable double ctheta, stheta, cphi, sphi,
                   p2, p, pt, px, py, pz, dpdE, dptdE, 
                   dpxdE, dpydE, dpzdE, dpxdtheta, dpydtheta,
                   chi2;
                   // d2pdE2, d2ptsE2;
  
};



#endif // __NEUTRINOFITOBJECT_H

