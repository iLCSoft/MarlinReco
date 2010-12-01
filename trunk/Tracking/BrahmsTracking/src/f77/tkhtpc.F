*  
* $Id: tkhtpc.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKHTPC
C****************************************************************************
C* CALL TKGHIT FOR TPC                                                      *
C*                                                                          *
C* Kristian Harder, September 1999                                          *
C****************************************************************************
      IMPLICIT NONE
#include "include/detnums.inc"
#include "include/turnon.inc"
#include "include/tpcres.inc"
#include "include/tpcdims.inc"

      IF (LTPC) THEN
C.......the rphi resolution will be overwritten by parametrisation later
        if(TPRP.eq.0.0) then
C For simple studies, use a general r-phi resolution and z-resolution (cm^2):
          TPRP=0.016
          TPZ=0.10
          TPRPSQ=TPRP**2
          TPZSQ=TPZ**2
          TPRRSQ=TPCRRS**2
          write(*,*) ' +++++++ TPC RES ++++++ reinitialization'
       endif
C        CALL TKGHIT('CTPC','TPCH',4,ID_TPC,3,TPRP,TPZ)

      END IF

      RETURN
      END


