/*! \file 
 *  \brief Declares class JetFitObject
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: JetFitObject.h,v $
 * - Revision 1.3  2011/03/03 15:03:02  blist
 * - Latest version, with NewFitterGSL
 * -
 * - Revision 1.2  2009/02/17 12:46:34  blist
 * - Improved version of NewtonFitterGSL, JetFitObject changed
 * -
 * - Revision 1.1  2008/02/12 10:19:05  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.9  2008/02/04 17:30:53  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.8  2008/01/30 21:48:02  blist
 * - Newton Fitter still doesnt work :-(
 * -
 * - Revision 1.7  2008/01/30 09:14:53  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.6  2008/01/29 17:23:00  blist
 * - new addTo2ndDerivatives and setParam
 * -
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
 * $Date: 2011/03/03 15:03:02 $
 * $Author: blist $
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

    
    /// Return a new copy of itself
    virtual JetFitObject *copy() const;
    
    /// Assign from anther object, if of same type
    virtual JetFitObject& assign (const BaseFitObject& source   ///< The source object
                                 );
    
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
    
    virtual void addToGlobalChi2DerVector (double *y,     ///< Vector of chi2 derivatives
                                           int idim,      ///< Vector size 
                                           double lambda, ///< The lambda value
                                           double der[]   ///< 4-vector with dg/dE, dg/dpx, dg/dpy, dg/dpz
                                           ) const;
    /// Calculates the squared error for a quantity with derivatives w.r.t. E, dx, dy, dz
    virtual double getError2 (double der[]    ///< Factors for d(E,px,py,pz)/dx_i
                             ) const;
    /// Get chi squared from measured and fitted parameters
    virtual double getChi2() const;

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
                   
    /// Adjust E, theta and phi such that E>=m, 0<=theta<=pi, -pi <= phi < pi; returns true if anything was changed             
    static bool adjustEThetaPhi (double& m, double &E, double& theta, double& phi);
    
    /// Calculate chi2 
    double calcChi2 () const;
    
};



#endif // __JETFITOBJECT_H

