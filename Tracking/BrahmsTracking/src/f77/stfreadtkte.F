
      REAL FUNCTION RSTFREADTKTE(ATTRIBUTE,ATE)
      IMPLICIT NONE
      integer ATTRIBUTE, ATE
      external rreadtktecpp 
      REAL rreadtktecpp

      RSTFREADTKTE = rreadtktecpp(ATTRIBUTE,ATE)

      RETURN
      END

      INTEGER FUNCTION ISTFREADTKTE(ATTRIBUTE,ATE)
      IMPLICIT NONE
      integer ATTRIBUTE, ATE
      external ireadtktecpp 
      integer ireadtktecpp

      ISTFREADTKTE = ireadtktecpp(ATTRIBUTE,ATE)

      RETURN
      END
