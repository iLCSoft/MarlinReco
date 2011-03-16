 
#ifndef __IRSPHOTONFITOBJECT_H
#define __ISRPHOTONFITOBJECT_H

#include "ParticleFitObject.h"

/// Class ISRPhotonFitObject
/// Documention: arXiv:1006.0436 [hep-ex]

/// Class for ISR photons with (p_x,p_y (both fix),  p_z (free)) in kinematic fits
/// p_z is internally replaced by a parameter p_g
///
/// This class assumes a photon p_z distribution according to dN/d|p_z| = c*|p_z|^(b-1), where the total number of photons should be given by N = \int_{PzMin}^{PzMax} dN/d|p_z| d|p_z|
/// The parameters b, PzMaxB:=PzMax^b and PzMinB:=PzMin^b are required to describe the photon spectrum
/// Only |p_z| values in [PzMin,PzMax[ can be used as start values (assertion in PgFromPz(...)) and will occur in the fit
/// Recommended start value is the missing p_z (fitter will always find a minimum around p_z=0)

/**
 *
 * Author: Moritz Beckmann
 * $Date: 2011/03/16 16:33:24 $
 * $Author: mbeckman $
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: ISRPhotonFitObject.h,v $
 * - Revision 1.2  2011/03/16 16:33:24  mbeckman
 * - Compatibility fixes with ILCSoft svn
 * -
 * - Revision 1.1  2010/06/11 20:32:51  mbeckman
 * - Renamed PhotonFitObjects, cleaned them up for usage
 * -
 * - Revision 1.6  2009/03/26 08:47:16  mbeckman
 * - Bug fix (measured p = 0 instead of start value), extended documentation
 * -
 * - Revision 1.5  2009/02/23 12:03:18  mbeckman
 * - - PhotonFitObject:     bug fix (1/0), removed dispensable variables
 * - - PhotonFitObjectPxyg: bug fixes (1/0, order of computing variables), modified parametrization
 * -
 * - Revision 1.4  2009/02/18 11:53:42  mbeckman
 * - documentation, debug output
 *
 */ 
class ISRPhotonFitObject : public ParticleFitObject {
  public:
    ISRPhotonFitObject(double px, double py, double pz,                   /// initial values for photon (p_x,p_y fix)
                        double b_, double PzMaxB_, double PzMinB_ = 0.);  /// photon spectrum parametrization (see above)
    virtual ~ISRPhotonFitObject();

    
    /// Return a new copy of itself
    virtual ISRPhotonFitObject *copy() const;
    
    /// Assign from anther object, if of same type
    virtual ISRPhotonFitObject& assign (const BaseFitObject& source   ///< The source object
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

    virtual void invalidateCache() const;
  
  protected:
    
    virtual double PgFromPz(double pz);
    
    virtual void initCov();
    
    void updateCache() const;
  
    mutable bool cachevalid;
    
    mutable double pt2, p2, p, pz,
                   dpx0, dpy0, dpz0, dE0, dpx1, dpy1, dpz1, dE1,
                   dpx2, dpy2, dpz2, dE2, d2pz22, d2E22,
                   chi2,                   
                   b, PzMinB, PzMaxB, dp2zFact;
};

#endif // __ISRPHOTONFITOBJECT_H
