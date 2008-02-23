////////////////////////////////////////////////////////////////
// Class NewtonFitter
//
// Author: Benno List
// Last update: $Date: 2008-02-23 11:18:39 $
//          by: $Author: listj $
// 
// Description: kinematic fit using Newton method
//               
////////////////////////////////////////////////////////////////

#ifndef __NEWTONFITTER_H
#define __NEWTONFITTER_H

#include<vector>
#include "BaseFitter.h"

class NewtonFitter : public BaseFitter {
  public:
    NewtonFitter();
    virtual ~NewtonFitter();
    virtual double fit();
    
    virtual int getError() const;
    virtual double getProbability() const;
    virtual double getChi2() const;
    virtual int    getDoF() const;
    virtual int  getIterations() const;

    virtual bool initialize();
  
  protected:
    virtual double calcChi2();
    
    void printMy (double M[], double y[], int idim);
    
    enum {NPARMAX=50, NCONMAX=10, NUNMMAX=10};
    
    int npar;      // total number of parameters
    int ncon;      // total number of hard constraints
    int nsoft;     // total number of soft constraints
    int nunm;      // total number of unmeasured parameters
    int ierr;      // Error status
    int nit;       // Number of iterations

    double fitprob;   // fit probability
    double chi2;      // final chi2

};

#endif // __NEWTONFITTER_H
