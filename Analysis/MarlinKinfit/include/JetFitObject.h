/*! \file 
 *  \brief Declares class JetFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.5  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.4  2007/09/13 08:09:50  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 
 
#ifndef __JETFITOBJECT_H
#define __JETFITOBJECT_H

#include "ParticleFitObject.h"

// Class JetFitObject
/// Class for jets with (E, eta, phi) in kinematic fits
/**
 *
 * Author: Jenny List, Benno List
 * $Date: 2007-10-30 15:51:14 $
 * $Author: gaede $
 *
 * \b Changelog:
 * - 30.12.04 BL: addToGlobCov, getDChi2DParam, getDChi2DParam2, addToGlobalChi2DerMatrix
 *            moved up to ParticleFitObject,
 *            getParamName implemented
 */ 
class JetFitObject : public ParticleFitObject {
  public:
    JetFitObject(double E, double theta, double phi, 
                 double DE, double Dtheta, double Dphi, 
                 double m = 0);
    virtual ~JetFitObject();
    
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
    
    /// add derivatives to vector der of size idim
    /// pxfact*dpx/dx_i + pyfact*dpy/dx_i + pzfact*dpz/dx_i + efact*dE/dx_i 
    virtual void   addToDerivatives (double der[],      ///< Derivatives vector, length idim
                                     int idim,          ///< Length of derivatives vector
                                     double efact=0,    ///< Factor for dE/dx_i
                                     double pxfact=0,   ///< Factor for dpx/dx_i
                                     double pyfact=0,   ///< Factor for dpy/dx_i 
                                     double pzfact=0    ///< Factor for dpz/dx_i
                                     ) const;

    /// add second order derivatives to matrix der2 of size idim x idim
    /// pxfact*d^2px/(dx_i dx_j) + pyfact...
    virtual void   addTo2ndDerivatives (double der2[],  ///< Derivatives vector, size idim x idim
                                        int    idim,    ///< First dimension of derivatives matrix
                                        double efact,   ///< Factor for d^2E/dx_i dx_j
                                        double pxfact,  ///< Factor for d^2px/dx_i dx_j
                                        double pyfact,  ///< Factor for d^2py/dx_i dx_j
                                        double pzfact   ///< Factor for d^2pz/dx_i dx_j
                                       ) const;
    
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



#endif // __JETFITOBJECT_H

