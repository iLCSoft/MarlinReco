
#ifndef __NEUTRINOFITOBJECT_H
#define __NEUTRINOFITOBJECT_H

#include "ParticleFitObject.h"

#include <cmath>


// Class NeutrinoFitObject
/// Class for neutrinos with (E, eta, phi) in kinematic fits
/**
 *
 * Author: Jenny List, Benno List
 * $Date: 2007-10-30 15:51:14 $
 * $Author: gaede $
 *
 * \b Changelog:
 * - 30.12.04 BL: addToGlobCov, getDChi2DParam, getDChi2DParam2,
 *            addToGlobalChi2DerMatrix moved up to ParticleFitObject,
 *            getParamName implemented
 */ 
class NeutrinoFitObject : public ParticleFitObject {
  public:
    NeutrinoFitObject(double E, double eta, double phi, 
                 double DE, double Deta, double Dphi, 
                 double m = 0);
    virtual ~NeutrinoFitObject();
    
    /// Get name of parameter ilocal
    virtual const char *getParamName (int ilocal     ///< Local parameter number
                                     ) const;
    
    // setters
    virtual bool setParam (int i, double par_, bool measured_ );  
    virtual bool setParam (int i, double par_ );  
    
    // these depend on actual parametrisation!
    virtual double getPx() const;
    virtual double getPy() const;
    virtual double getPz() const;
    virtual double getE() const;
    virtual double getDPx(int ilocal) const;
    virtual double getDPy(int ilocal) const;
    virtual double getDPz(int ilocal) const;
    virtual double getDE(int ilocal) const;
    
    virtual double getD2Px(int ilocal1, int ilocal2) const;
    virtual double getD2Py(int ilocal1, int ilocal2) const;
    virtual double getD2Pz(int ilocal1, int ilocal2) const;
    virtual double getD2E (int ilocal1, int ilocal2) const;
    
    virtual void addToGlobalDerMatrix (int idim, double c, double *M) const;
      
  protected:
    inline  double getP() const;
    

};

// Implementation of inline methods
double NeutrinoFitObject::getP() const {return sqrt(std::abs(par[0]*par[0]-mass*mass));}

#endif // __NEUTRINOFITOBJECT_H

