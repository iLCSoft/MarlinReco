#ifndef ALGEBRAIMPLEMENTATION_H_1
#define ALGEBRAIMPLEMENTATION_H_1


#include "lcio.h"
#include <EVENT/ReconstructedParticle.h>
#include "TMatrixD.h"

int getCovMatrixMomenta(EVENT::ReconstructedParticle const *, TMatrixD &);

int getCovMatrixMomenta(EVENT::ReconstructedParticle const *, EVENT::FloatVec &);


#endif
