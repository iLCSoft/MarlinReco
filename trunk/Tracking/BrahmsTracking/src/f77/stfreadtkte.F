
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
      integer  ireadtktecpp

      ISTFREADTKTE = ireadtktecpp(ATTRIBUTE,ATE)

      RETURN
      END

      INTEGER FUNCTION STFREADTKITEDAT(ATTRIBUTE,ATE)
      IMPLICIT NONE
      integer ATTRIBUTE, ATE
      external readtkitedatcpp 
      integer readtkitedatcpp

      STFREADTKITEDAT = readtkitedatcpp(ATTRIBUTE,ATE)
      
      RETURN
      END
