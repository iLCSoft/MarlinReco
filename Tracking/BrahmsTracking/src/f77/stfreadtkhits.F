
      REAL FUNCTION STFREADTKHITS(ATTRIBUTE,AHIT)
      IMPLICIT NONE
      integer ATTRIBUTE, AHIT
      external readtkhitscpp 
      REAL readtkhitscpp

      STFREADTKHITS = readtkhitscpp(ATTRIBUTE,AHIT)

      RETURN
      END


