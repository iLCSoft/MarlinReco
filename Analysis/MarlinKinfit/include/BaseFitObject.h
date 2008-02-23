/*! \file 
 *  \brief Declares class BaseFitObject
 *
 * \b Changelog:
 * - 7.6.04 JB: First doxygen docu
 *
 * \b CVS Log messages:
 * - $Log: not supported by cvs2svn $
 * - Revision 1.10  2008/02/07 08:15:25  blist
 * - error calculation of constraints fixed
 * -
 * - Revision 1.9  2008/02/04 17:30:53  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.8  2008/01/30 09:14:53  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.7  2008/01/29 17:17:33  blist
 * - implemented setname
 * -
 * - Revision 1.6  2007/09/17 12:50:15  blist
 * - Some parameters reordered
 * -
 * - Revision 1.5  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 * - Revision 1.4  2007/09/13 13:33:06  blist
 * - Print methods return os
 * -
 */ 

#ifndef __BASEFITOBJECT_H
#define __BASEFITOBJECT_H

#include <iostream>


// Class BaseFitObject
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
 * In its current state, a BaseFitObject has a set of parameters, some
 * of them measured (i.e., they contribute to the \f$\chi^2\f$).
 * These parameters have a local numbering, running from 0 to n-1.
 * Global numbers can be assigned by the BaseFitter using
 * setGlobalParNum. 
 *
 * The class WWFitter needs the following routines from BaseFitObject:
 * - BaseFitObject::getNPar
 * - BaseFitObject::getMeasured
 * - BaseFitObject::getGlobalParNum
 * - BaseFitObject::getParam
 * - BaseFitObject::setParam
 * - BaseFitObject::addToGlobCov
 * - BaseFitObject::operator<<
 * - BaseFitObject::setError
 *
 *
 * 
 *
 * Author:  Benno List, Jenny Böhme
 * $Date: 2008-02-23 11:18:39 $
 * $Author: listj $
 *
 * \b Changelog:
 * - 30.12.04 BL: Added getCov, setCov
 * - 11.1.05 BL: return type of setParam changed to bool
 */ 


class BaseFitObject {
  public:
    BaseFitObject();
    /// Copy constructor
    BaseFitObject (const BaseFitObject& rhs              ///< right hand side
                   );
    /// Assignment               
    BaseFitObject& operator= (const BaseFitObject& rhs   ///< right hand side
                             );
    
    virtual ~BaseFitObject();
    
    /// Set value and measured flag of parameter i; return: significant change
    virtual bool setParam (int ilocal,         ///< Local parameter number
                           double par_,        ///< New parameter value
                           bool measured_,     ///< New "measured" flag
                           bool fixed_ = false ///< New "fixed" flag
                          ) = 0;  
    /// Set value of parameter ilocal; return: significant change
    virtual bool setParam (int ilocal,    ///< Local parameter number
                           double par_    ///< New parameter value
                          ) = 0;  
    /// Read values from global vector, readjust vector; return: significant change
    virtual bool updateParams (double p[],   ///< The parameter vector
                               int idim      ///< Length of the vector                         
                              );  
    /// Set measured value of parameter ilocal; return: success
    virtual bool setMParam (int ilocal,    ///< Local parameter number
                            double mpar_   ///< New measured parameter value
                           ) = 0;  
    /// Set error of parameter ilocal; return: success
    virtual bool setError (int ilocal,    ///< Local parameter number
                           double err_    ///< New error value
                           ) = 0;
    /// Set covariance of parameters ilocal and jlocal; return: success
    virtual bool setCov (int ilocal,    ///< Local parameter number
                         int jlocal,    ///< Local parameter number
                         double cov_    ///< New error value
                        ) = 0;
    /// Set number of parameter ilocal in global list
    /// return true signals OK
    virtual bool setGlobalParNum (int ilocal,  ///< Local parameter number
                                  int iglobal  ///< New global parameter number
                                  ) = 0; 
    
    /// Fix a parameter (fix=true), or release it (fix=false)
    virtual bool fixParam (int ilocal,    ///< Local parameter number
                           bool fix=true  ///< fix if true, release if false
                          ) = 0;
    /// Release a parameter 
    virtual bool releaseParam (int ilocal    ///< Local parameter number
                              ) 
    { return fixParam (ilocal, false); } 
    
    /// Returns whether parameter is fixed 
    virtual bool isParamFixed (int ilocal     ///< Local parameter number
                              ) const = 0; 
    
    /// Get current value of parameter ilocal
    virtual double getParam (int ilocal     ///< Local parameter number
                            ) const = 0;
    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                            ) const { return "???";}
    /// Get object's name
    virtual const char *getName () const { return name ? name : "???";}
    /// Set object's name
    virtual void setName (const char * name_);
    /// Get measured value of parameter ilocal
    virtual double getMParam (int ilocal     ///< Local parameter number
                             ) const = 0;
    /// Get error of parameter ilocal
    virtual double getError (int ilocal     ///< Local parameter number
                            ) const = 0;
    /// Get covariance between parameters ilocal and jlocal
    virtual double getCov (int ilocal,    ///< Local parameter number i
                           int jlocal     ///< Local parameter number j
                          ) const = 0;
    /// Get measured flag for parameter ilocal
    virtual bool isParamMeasured (int ilocal    ///< Local parameter number
                                   ) const = 0;
    /// Get global parameter number of parameter ilocal
    virtual int getGlobalParNum(int ilocal     ///< Local parameter number
                                ) const = 0;
    /// Get total number of parameters of this FitObject
    virtual int getNPar() const = 0;
    /// Get number of measured parameters of this FitObject
    virtual int getNMeasured() const;
    /// Get number of unmeasured parameters of this FitObject
    virtual int getNUnmeasured() const;
    /// Get number of free parameters of this FitObject
    virtual int getNFree() const;
    /// Get number of fixed parameters of this FitObject
    virtual int getNFixed() const;
    
    /// Add covariance matrix elements to 
    /// global covariance matrix of size idim x idim
    virtual void addToGlobCov(double *glcov, int idim) const = 0; 
            
    /// Get chi squared from measured and fitted parameters
    virtual double getChi2() const = 0;
    /// Get derivative of chi squared w.r.t. parameter ilocal
    virtual double getDChi2DParam(int ilocal   ///< Local parameter number
                                   ) const = 0;
    /// Get second derivative of chi squared w.r.t. parameters ilocal1 and ilocal2
    virtual double getD2Chi2DParam2(int ilocal,   ///< Local parameter number i
                                    int jlocal    ///< Local parameter number j
                                   ) const = 0;
    
    /// Add derivatives of chi squared to global covariance matrix
    virtual void addToGlobalChi2DerMatrix (double *M,   ///< Global covariance matrix
                                           int idim     ///< First dimension of global covariance matrix
                                           ) const = 0;
    /// Add derivatives of chi squared to global derivative matrix
    virtual void addToGlobalChi2DerVector (double *y,   ///< Vector of chi2 derivatives
                                           int idim     ///< Vector size 
                                           ) const = 0;
    /// print the parameters
    virtual std::ostream& printParams (std::ostream& os  ///< The output stream
                                      ) const;
    
    /// print object to ostream
    virtual std::ostream&  print (std::ostream& os       ///< The output stream
                                 ) const = 0;
    /// invalidate any cached quantities
    virtual void invalidateCache() const {};
                                 
    protected:
      char *name;  
      const static double eps2;                           
};

/** \relates BaseFitObject
 *  \brief Prints out a BaseFitObject, using its print method
 */
inline std::ostream& operator<< (std::ostream& os,         ///< The output stream
                                 const BaseFitObject& bfo  ///< The object to print
                                ) {
  return bfo.print(os);
}

#endif // __BASEFITOBJECT_H

