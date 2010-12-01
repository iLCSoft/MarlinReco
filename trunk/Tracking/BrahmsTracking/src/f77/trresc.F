*  
* $Id: trresc.F,v 1.2 2005-08-08 07:09:13 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
* Revision 1.2  2005/06/20 18:49:46  aplin
* change tpc geometery so that it is read from  same header file as cpp i.e. marlin_tpcgeom.h
*
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
      SUBROUTINE TRRESC
      IMPLICIT NONE
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C
C----------------------------------------------------------------*
C! Rescue a chain,
C! break it in half and refit
C!
C!  M.Mermikides  6/10/86   V 1.0
C!
C!  Called by TRKFND
C!
C----------------------------------------------------------------*
#include "include/fndparms.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/tchain.inc"
      INTEGER I,N1,N2,NC
      SAVE
C
      NC = NHCHAI
      N1 = NC/2
      N2 = NC - N1
      NHCHAI = N1
      CALL TRFTCH
      IF (ISTCHA.EQ.0) THEN
         DO 10 I=1,N2
            IRCHAI(I) = IRCHAI(N1+I)
            IHCHAI(I) = IHCHAI(N1+I)
            D2CHAI(I) = D2CHAI(N1+I)
 10      CONTINUE
         NHCHAI = N2
         CALL TRFTCH
      ENDIF
      RETURN
      END

c-----------------------------------------------------------------------
