////////////////////////////////////////////////////////////////
// Class PyConstraint
//
// Author: Jenny Boehme
// Last update: $Date: 2007-10-30 15:51:14 $
//          by: $Author: gaede $
// 
// Description: sum (py) = 0 constraint
//               
////////////////////////////////////////////////////////////////

#ifndef __PYCONSTRAINT_H
#define __PYCONSTRAINT_H

#include "ParticleConstraint.h"

class ParticleFitObject;

class PyConstraint : public ParticleConstraint {
  public:
    PyConstraint ();
    virtual ~PyConstraint();
    virtual double getValue() const;
    virtual void getDerivatives (int idim, double der[]) const;
  
  protected:

};

#endif // __PYCONSTRAINT_H
