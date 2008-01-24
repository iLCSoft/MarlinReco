

#ifndef PADROW_INC
#define PADROW_INC

C     The maximum number of TPC Padrows that the tracking code 
C     can accept must be a multiple of 32
C     Note this is not the number of TPC Padrows in the dectector
C     but the maximum allow becauase of the predifined array sizes
C     and a set of 32 bit registers 

C     For example 512 as below

C     By increasing the value of N32BITREG you can increase the max 
C     number of pad rows

      INTEGER N32BITREG
      PARAMETER (N32BITREG=16)

      INTEGER LTPDRO
      PARAMETER (LTPDRO=N32BITREG*32) 
C     here the max number of TPC rows is set to N32BITREG*32=256


     
      INTEGER LMXCIR, LMXZCR
      PARAMETER (LMXCIR=LTPDRO,LMXZCR=LTPDRO)

#endif /* PADROW_INC */ es
