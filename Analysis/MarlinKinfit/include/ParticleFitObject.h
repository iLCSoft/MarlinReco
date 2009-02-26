/*! \file 
 *  \brief Declares class ParticleFitObject
 *
 * \b Changelog:
 * - 17.11.04 BL: First version (refactured from BaseFitObject)
 *
 */ 

#ifndef __PARTICLEFITOBJECT_H
#define __PARTICLEFITOBJECT_H

#include "BaseFitObject.h"


// Class ParticleFitObject
/// Abstract base class for particle objects of kinematic fits
/**
 * This class defines the minimal functionality any fit object must provide. 
 * The main task of a fit object is to keep parameters (and errors) that define 
 * the four-momentum of a particle and encapsulate the actually chosen 
 * parametrisation from the rest of the fitting machinery.
 *
 * Since for the fit a parametrisation distributed like a gaussian is most
 * favorable, different kinds of particles (implying different kinds of
 * measurements!) might require different parametrisations. For each 
 * desired parametrisation a concrete class should be derived from this 
 * abstract base class. It needs to be able to convert its parameters
 * to E, px, py, pz and to provide the derivatives of E, px, py, pz
 * with respect to the internal parameters.
 *
 * Depending on the type of particle, some or all parameters might
 * be unmeasured (neutrinos!), meaning that they come with a very large
 * and/or unknown error. They are treated differently by the 
 * fit algorithm and are thus flagged accordingly.
 *
 * In order to insert its derivatives into the global covariance matrix
 * of all FitObjects in the event, each FitObjects needs to know the
 * position of its parameters in the overall parameter list. 
 *
 * THIS iS JUNK!!!! It is done like this in JetFitObject.C,
 * but using measured[i] which is the bool giving the measured/unmeasured
 * status and NOT a bool containing the START VALUES!!!!!
 * From its stored initial parameters and the current fit parameters
 * the FitObject calculates its contribution to the $\chi^2$ of the fit.
 *
 * In its current state, a ParticleFitObject has a set of parameters, some
 * of them measured (i.e., they contribute to the \f$\chi^2\f$).
 * These parameters have a local numbering, running from 0 to n-1.
 * Global numbers can be assigned by the BaseFitter using
 * setGlobalParNum. 
 *
 * Author: Benno List, Jenny List
 * $Date: 2009-02-26 18:35:17 $
 * $Author: beckmann $
 *
 * \b Changelog:
 * - 30.12.04 BL: Added getCov, setCov
 *                addToGlobCov, getDChi2DParam, getDChi2DParam2, addToGlobalDerMatrix implemented
 */ 


class ParticleFitObject: public BaseFitObject {
  public:
    /// Default constructor
    ParticleFitObject();
    /// Virtual destructor
    virtual ~ParticleFitObject();
    
    /// Set value and measured flag of parameter ilocal; return=significant change
    virtual bool   setParam (int ilocal,         ///< Local parameter number
                             double par_,        ///< New parameter value
                             bool measured_,     ///< New "measured" flag
                             bool fixed_ = false ///< New "fixed" flag
                            );  
    /// Set value of parameter ilocal; return=significant change
    virtual bool   setParam (int ilocal,    ///< Local parameter number
                             double par_    ///< New parameter value
                             );  
                             
    /// Set measured value of parameter ilocal; return=success
    virtual bool   setMParam (int i, double mpar_ );  
    /// Set error of parameter ilocal; return=success
    virtual bool   setError (int ilocal,    ///< Local parameter number
                             double err_    ///< New error value
                             );
    /// Set covariance of parameters ilocal and jlocal; return=success
    virtual bool   setCov (int ilocal,    ///< Local parameter number
                           int jlocal,    ///< Local parameter number
                           double cov_    ///< New error value
                          );
    /// Set mass of particle; return=success
    virtual bool setMass (double mass_);
    /// Get mass of particle
    virtual double getMass () const;
    /// Set number of parameter ilocal in global list
    /// return true signals OK
    virtual bool setGlobalParNum (int ilocal, int iglobal); 
    /// Fix a parameter (fix=true), or release it (fix=false)
    virtual bool fixParam (int ilocal,    ///< Local parameter number
                           bool fix=true  ///< fix if true, release if false
                          );
                          
    
    /// Get covariance matrix
    virtual double *getCovMatrix() {return cov[0];};
    /// Get current value of parameter ilocal
    virtual double getParam (int ilocal ///< Local parameter number
                            ) const;
    /// Get measured value of parameter ilocal
    virtual double getMParam (int iocal   ///< Local parameter number
                            ) const;
    /// Get error of parameter ilocal
    virtual double getError (int ilocal     ///< Local parameter number
                            ) const;
    /// Get covariance between parameters ilocal and jlocal
    virtual double getCov (int ilocal,    ///< Local parameter number i
                           int jlocal     ///< Local parameter number j
                          ) const;
    /// Get measured flag for parameter i
    virtual bool isParamMeasured (int ilocal) const;
    /// Get fixed flag for parameter i
    virtual bool isParamFixed (int ilocal) const;
    /// Get global parameter number of parameter ilocal
    virtual int getGlobalParNum(int ilocal) const;
    /// Get number of parameters of this FitObject
    virtual int getNPar() const {return NPAR;};
    
    /// Add covariance matrix elements to global covariance matrix
    virtual void addToGlobCov(double *globCov,    ///< Global covariance matrix, size idim x idim
                              int idim            ///< First dimension of globCov
                             ) const; 
    
    /// print the four-momentum (E, px, py, pz)
    virtual std::ostream& print4Vector (std::ostream& os   ///< The output stream
                                       ) const;
        
    /// Return E
    virtual double getE() const = 0;
    /// Return px
    virtual double getPx() const = 0;
    /// Return py
    virtual double getPy() const = 0;
    /// Return pz
    virtual double getPz() const = 0;
    
    /// Return p (momentum) 
    virtual double getP() const = 0;
    /// Return p (momentum) squared
    virtual double getP2() const = 0;
    /// Return pt (transverse momentum) 
    virtual double getPt() const = 0;
    /// Return pt (transverse momentum)  squared
    virtual double getPt2() const = 0;
    
    /// Return d p_x / d par_ilocal (derivative of px w.r.t. local parameter ilocal)
    virtual double getDPx (int ilocal              ///< Local parameter number
                           ) const = 0;
    /// Return d p_y / d par_ilocal (derivative of py w.r.t. local parameter ilocal)
    virtual double getDPy (int ilocal              ///< Local parameter number
                           ) const = 0;
    /// Return d p_z / d par_ilocal (derivative of pz w.r.t. local parameter ilocal)
    virtual double getDPz (int ilocal              ///< Local parameter number
                           ) const = 0;
    /// Return d E / d par_ilocal (derivative of E w.r.t. local parameter ilocal)
    virtual double getDE (int ilocal               ///< Local parameter number
                          ) const = 0;
    
    
    /// Return all derivatives w.r.t. local parameters
    /**
     *  The result vector der will contain d E/d par1, d px/d par 1, d py/d par 1, d pz / d par1, d E/d par 2 ...
     */
    virtual void getDerivatives (double der[],    ///< Derivatives vector, length 4*idim
                                 int idim         ///< Length of derivatives vector
                                ) const;    
    
    /// add derivatives to vector der of size idim
    /// pxfact*dpx/dx_i + pyfact*dpy/dx_i + pzfact*dpz/dx_i + efact*dE/dx_i 
    virtual void   addToDerivatives (double der[],      ///< Derivatives vector, length idim
                                     int    idim,       ///< Length of derivatives vector
                                     double efact=0,    ///< Factor for dE/dx_i
                                     double pxfact=0,   ///< Factor for dpx/dx_i
                                     double pyfact=0,   ///< Factor for dpy/dx_i 
                                     double pzfact=0    ///< Factor for dpz/dx_i
                                     ) const = 0;
    /// add second order derivatives to matrix der2 of size idim x idim
    /// pxfact*d^2px/(dx_i dx_j) + pyfact...
    virtual void   addTo2ndDerivatives (double der2[],  ///< Derivatives vector, size idim x idim
                                        int    idim,    ///< First dimension of derivatives matrix
                                        double efact,   ///< Factor for d^2E/dx_i dx_j
                                        double pxfact,  ///< Factor for d^2px/dx_i dx_j
                                        double pyfact,  ///< Factor for d^2py/dx_i dx_j
                                        double pzfact   ///< Factor for d^2pz/dx_i dx_j
                                       ) const = 0;
    /// Add second order derivatives to matrix M of size idim x idim
    /// der[0]*d^2E/(dx_i dx_j) + der[1]*d^2px/(dx_i dx_j) + ...
    virtual void   addTo2ndDerivatives (double M[],     ///< Global derivatives matrix size idim x idim
                                        int    idim,    ///< First dimension of derivatives matrix
                                        double lamda,   ///< Global factor
                                        double der[]    ///< Factors for d^2(E,px,py,pz)/dx_i dx_j
                                       ) const = 0;
    /// Add first order derivatives to matrix M of size idim x idim
    /// der[0]*dE/dx_i + der[1]*dpx/dx_i + ...
    virtual void   addTo1stDerivatives (double M[],     ///< Global derivatives matrix size idim x idim
                                        int    idim,    ///< First dimension of derivatives matrix
                                        double der[],   ///< Factors for d^2(E,px,py,pz)/dx_i dx_j
                                        int kglobal     ///< Global parameter number of constraint
                                       ) const = 0;
    /// Calculates the squared error for a quantity with derivatives w.r.t. E, dx, dy, dz
    virtual double getError2 (double der[]    ///< Factors for d(E,px,py,pz)/dx_i
                             ) const = 0;
    
    /// Get chi squared from measured and fitted parameters
    virtual double getChi2() const;
    /// Get derivative of chi squared w.r.t. parameter ilocal
    virtual double getDChi2DParam (int ilocal           ///< Local parameter number
                                  ) const;
    /// Get second derivative of chi squared w.r.t. parameters ilocal and jlocal
    virtual double getD2Chi2DParam2(int ilocal,         ///< Local parameter number i
                                    int jlocal          ///< Local parameter number j
                                   ) const;
    
    /// Add derivatives of chi squared to global covariance matrix
    virtual void addToGlobalChi2DerMatrix (double *M,   ///< Global covariance matrix
                                           int idim     ///< First dimension of global covariance matrix
                                          ) const;
                                      
    /// Add  numerically determined derivatives of chi squared to global covariance matrix
    virtual void addToGlobalChi2DerMatrixNum (double *M,   ///< Global covariance matrix
                                           int idim,    ///< First dimension of global covariance matrix
                                           double eps   ///< Parameter variation
                                          );
                                      
    /// Add derivatives of chi squared to global derivative vector
    virtual void addToGlobalChi2DerVector (double *y,   ///< Vector of chi2 derivatives
                                           int idim     ///< Vector size 
                                           ) const;
    /// Add numerically determined derivatives of chi squared to global derivative vector
    virtual void addToGlobalChi2DerVectorNum (double *y,    ///< Vector of chi2 derivatives
                                              int idim,     ///< Vector size 
                                              double eps    ///< Parameter variation
                                             );


    /// Add derivatives of momentum vector to global derivative vector
    virtual void addToGlobalChi2DerVector (double *y,     ///< Vector of chi2 derivatives
                                           int idim,      ///< Vector size 
                                           double lambda, ///< The lambda value
                                           double der[]   ///< 4-vector with dg/dE, dg/dpx, dg/dpy, dg/dpz
                                           ) const = 0;
        
    
    /// print object to ostream
    virtual std::ostream& print (std::ostream& os    ///< The output stream
                                ) const;
                                
    void test1stDerivatives ();
    void test2ndDerivatives ();
                                
    /// Evaluates numerically the 1st derivative  of chi2 w.r.t. a parameter
    double num1stDerivative (int ilocal,  ///< Local parameter number 
                             double eps   ///< variation of  local parameter 
                            );
    /// Evaluates numerically the 2nd derivative of chi2 w.r.t. 2 parameters
    double num2ndDerivative (int ilocal1, ///< 1st local parameter number 
                             double eps1, ///< variation of 1st local parameter 
                             int ilocal2, ///< 1st local parameter number 
                             double eps2  ///< variation of 2nd local parameter 
                            );
  
  protected:
    /// Calculate the inverse of the covariance matrix
    virtual bool calculateCovInv() const;
        
    /// number of parameters
    enum {NPAR = 3};
    /// mass of particle
    double mass;
    /// fit parameters
    double par[NPAR];
    /// measured parameters
    double mpar[NPAR];
    /// measured flag
    bool measured[NPAR];
    /// fixed flag
    bool fixed[NPAR];
    /// global paramter number for each parameter
    int globalParNum [NPAR];
    /// local covariance matrix
    double cov [NPAR][NPAR];    
    /// inverse pf local covariance matrix
    mutable double covinv [NPAR][NPAR];    
    /// flag for valid inverse covariance matrix
    mutable bool covinvvalid; 

};

#endif // __PARTICLEFITOBJECT_H

