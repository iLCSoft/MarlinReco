*  
* $Id: tkveto.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:51  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKVETO(IVETO)
**********************************************************************
*
*     subroutine TKVETO
*
*     book-keeping event veto requests for performance measurement
*     NB: This is not meant for cheating, but to avoid simple technical
*         problems like arrays that are too small in just a few cases out
*         of 1000s spoiling the performance measurement.
*
*     Kristian Harder, January 2001
*
*     INPUT: iveto = 0 : clear veto flag
*                    1 : set event veto. event will not be used in
*                        performance measurement
*                   -1 : get veto flag state
*     OUTPUT: (only if iveto was -1 on input)
*                    0 : event veto not set
*                    1 : event veto set
*
**********************************************************************
      IMPLICIT NONE
#include "include/tkffread.inc"

      INTEGER IVETO,IVSAVE
      SAVE    IVSAVE

      IF (IVETO.EQ.0) THEN
         IVSAVE=0
      ELSE IF (IVETO.GT.0) THEN
         IF (TKSTDBG.GT.0) THEN
            IF (IVSAVE.EQ.0) THEN
              WRITE(6,*) 'TKVETO: EVENT VETO, NO PERFORMANCE ANALYSIS'
            ELSE
              WRITE(6,*) 'TKVETO: same event vetoed again'
            END IF
         END IF
         IVSAVE=1
      ELSE
         IVETO=IVSAVE
      END IF

      RETURN
      END


