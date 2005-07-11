
      REAL FUNCTION STFREADTKTE(ATTRIBUTE,ATE)
      IMPLICIT NONE
      integer ATTRIBUTE, ATE
      external readtktecpp 
      REAL readtktecpp

      STFREADTKTE = readtktecpp(ATTRIBUTE,ATE)

      RETURN
      END


