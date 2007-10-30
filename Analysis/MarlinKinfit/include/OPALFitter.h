////////////////////////////////////////////////////////////////
// Class OPALFitter
//
// Author: Jenny Boehme
// Last update: $Date: 2007-10-30 15:51:14 $
//          by: $Author: gaede $
// 
// Description: kinematic fit a la WWFGO
//               
////////////////////////////////////////////////////////////////

#ifndef __OPALFITTER_H
#define __OPALFITTER_H

#include<vector>
#include "BaseFitter.h"

class BaseFitObject;
class BaseConstraint;
/** Description of the fit algorithm and interface:
 *
 * The OPALFitter object holds a set (fitobjects) of BaseFitObject 
 * objects,
 * which represent the objects (particles, jets) whose fourvectors
 * shall be fitted.
 * It also holds a set (constraints) of BaseConstraint objects
 * that represent the constraints (momentum or mass constraints).
 *
 * OPALFitter::initialize goes over the list of fit objects and
 * determines the number of parameters (measured and unmeasured),
 * and assigns global numbers to them.
 *
 * OPALFitter::fit first initializes the global parameter numbering
 * using OPALFitter::initialize.
 *
 * Then it initializes vectors eta (\f$ \eta\f$)  and y (\f$ \vec y\f$)
 * with the current  parameter values.
 *
 * Next it initializes the matrix dfeta that represents 
 * \f$d \vec F / d \vec \eta\f$
 *
 * Used methods in OPALFitter::initialize:
 * - BaseFitObject::getNPar
 * - BaseFitObject::getMeasured
 * - BaseFitObject::setGlobalParNum
 * 
 * Used methods in OPALFitter::updateFitObjects:
 * - BaseFitObject::getGlobalParNum
 * - BaseFitObject::setParam
 *
 * Used methods in OPALFitter::fit:
 * - BaseFitObject::getNPar
 * - BaseFitObject::getGlobalParNum
 * - BaseFitObject::getParam
 * - BaseFitObject::addToGlobCov
 * - BaseFitObject::operator<<
 * - BaseFitObject::setError
 * - BaseConstraint::getDerivatives
 * - BaseConstraint::getValue
 *
 * Interaction between Constraints and Particles
 *
 * The Fitter does not care how the constraints and BaseFitObject objects
 * interact. The constraints are expressed in terms of four-vector
 * components of the BaseFitObjects.
 * A constraint object keeps a set of BaseFitObject objects, and, using
 * the chain rule, calculates its derivatives w.r.t. the global parameters.
 * 
 */

class OPALFitter : public BaseFitter {
  public:
    OPALFitter();
    virtual ~OPALFitter();
    virtual double fit();
    
    /// Return error code
    /** Error code meanings:
     *  - 0: No error
     *  - 1: out of iterations
     *  - 2: crazy chi^2
     *  - 3: minimum step size reached, chiK still increasing
     *  - 4: (step size decreased)
     *  - 5: (keep going)
     */
    virtual int getError() const;
    virtual double getProbability() const;
    virtual double getChi2() const;
    virtual int  getIterations() const;
  
  protected:
    virtual bool initialize();
    virtual bool updateFitObjects (double eta[]);
    enum {NPARMAX=50, NCONMAX=20, NUNMMAX=20};
    
    int npar;      // total number of parameters
    int nmea;      // total number of measured parameters
    int nunm;      // total number of unmeasured parameters
    int ncon;      // total number of constraints
    int ierr;      // Error status
    int nit;       // Number of iterations

    double fitprob;   // fit probability
    double chi2;      // final chi2
};

#endif // __OPALFITTER_H
