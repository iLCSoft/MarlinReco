
      REAL FUNCTION RSTFREADTKHITS(ATTRIBUTE,AHIT)
      IMPLICIT NONE
      integer ATTRIBUTE, AHIT
      external rreadtkhitscpp 
      REAL rreadtkhitscpp

      RSTFREADTKHITS = rreadtkhitscpp(ATTRIBUTE,AHIT)
c      write(*,*) "RSTFREADTKHITS = " , RSTFREADTKHITS 

      RETURN
      END

      INTEGER FUNCTION ISTFREADTKHITS(ATTRIBUTE,AHIT)
      IMPLICIT NONE
      integer ATTRIBUTE, AHIT
      external ireadtkhitscpp 
      INTEGER ireadtkhitscpp

      ISTFREADTKHITS = ireadtkhitscpp(ATTRIBUTE,AHIT)
c      write(*,*) "ISTFREADTKHITS = " , ISTFREADTKHITS 

      RETURN
      END


