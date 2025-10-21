#ifndef ALGEBRAIMPLEMENTATION_H_1
#define ALGEBRAIMPLEMENTATION_H_1

#include "TMatrixD.h"
#include "lcio.h"
#include <EVENT/ReconstructedParticle.h>

int getCovMatrixMomenta(EVENT::ReconstructedParticle const*, TMatrixD&);

int getCovMatrixMomenta(EVENT::ReconstructedParticle const*, EVENT::FloatVec&);

#endif
