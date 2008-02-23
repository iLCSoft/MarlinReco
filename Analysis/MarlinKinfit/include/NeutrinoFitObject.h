
#ifndef __NEUTRINOFITOBJECT_H
#define __NEUTRINOFITOBJECT_H

#include "ParticleFitObject.h"

#include <cmath>


// Class NeutrinoFitObject
/// Class for neutrinos with (E, eta, phi) in kinematic fits
/**
 *
 * Author: Jenny List, Benno List
 * $Date: 2008-02-23 11:18:39 $
 * $Author: listj $
 *
 * \b Changelog:
 * - 30.12.04 BL: addToGlobCov, getDChi2DParam, getDChi2DParam2,
 *            addToGlobalChi2DerMatrix moved up to ParticleFitObject,
 *            getParamName implemented
 */ 
class NeutrinoFitObject : public ParticleFitObject {
  public:
    NeutrinoFitObject(double E, double eta, double phi, 
                      double DE=1, double Dtheta=0.1, double Dphi=0.1);
    virtual ~NeutrinoFitObject();
    
    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                     ) const;
    
    /// Set value of parameter ilocal; return: significant change
    virtual bool setParam (int ilocal,    ///< Local parameter number
                           double par_    ///< New parameter value
                          );  
    /// Set value and measured flag of parameter ilocal; return=significant change
    virtual bool   setParam (int ilocal,         ///< Local parameter number
                             double par_,        ///< New parameter value
                             bool measured_,     ///< New "measured" flag
                             bool fixed_ = false ///< New "fixed" flag
                            );  
    /// Read values from global vector, readjust vector; return: significant change
    virtual bool   updateParams (double p[],   ///< The parameter vector
                                 int idim      ///< Length of the vector                         
                                );  
    
    // these depend on actual parametrisation!
    virtual double getPx() const;
    virtual double getPy() const;
    virtual double getPz() const;
    virtual double getE() const;
    virtual double getPt() const;
    virtual double getP2() const;
    virtual double getPt2() const;
    virtual double getDPx(int ilocal) const;
    virtual double getDPy(int ilocal) const;
    virtual double getDPz(int ilocal) const;
    virtual double getDE(int ilocal) const;
    
/*
 *     virtual double getD2Px(int ilocal1, int ilocal2) const;
 *     virtual double getD2Py(int ilocal1, int ilocal2) const;
 *     virtual double getD2Pz(int ilocal1, int ilocal2) const;
 *     virtual double getD2E (int ilocal1, int ilocal2) const;
 *     
 */
    virtual void invalidateCache() const;

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
                                       
    virtual void addToGlobalChi2DerVector (double *y,     ///< Vector of chi2 derivatives
                                           int idim,      ///< Vector size 
                                           double lambda, ///< The lambda value
                                           double der[]   ///< 4-vector with dg/dE, dg/dpx, dg/dpy, dg/dpz
                                           ) const;

    
    /// Add second order derivatives to matrix M of size idim x idim
    /// der[0]*d^2E/(dx_i dx_j) + der[1]*d^2px/(dx_i dx_j) + ...
    virtual void   addTo2ndDerivatives (double M[],     ///< Global derivatives matrix size idim x idim
                                        int    idim,    ///< First dimension of derivatives matrix
                                        double lamda,   ///< Global factor
                                        double der[]    ///< Factors for d^2(E,px,py,pz)/dx_i dx_j
                                       ) const;
    /// Add first order derivatives to matrix M of size idim x idim
    /// der[0]*dE/dx_i + der[1]*dpx/dx_i + ...
    virtual void   addTo1stDerivatives (double M[],     ///< Global derivatives matrix size idim x idim
                                        int    idim,    ///< First dimension of derivatives matrix
                                        double der[],   ///< Factors for d^2(E,px,py,pz)/dx_i dx_j
                                        int kglobal     ///< Global parameter number of constraint
                                       ) const;
    /// Calculates the squared error for a quantity with derivatives w.r.t. E, dx, dy, dz
    virtual double getError2 (double der[]    ///< Factors for d(E,px,py,pz)/dx_i
                             ) const;
      
  protected:
    inline  double getP() const;
    
    void updateCache() const;
    
  
    mutable bool cachevalid;
    
    mutable double ctheta, stheta, cphi, sphi,
                   pt, px, py, pz, dptdE, 
                   dpxdE, dpydE, dpxdtheta, dpydtheta,
                   chi2;

};

#endif // __NEUTRINOFITOBJECT_H

