*  
* $Id: tkinit.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE TKINIT
**********************************************************************
*
*     subroutine TKINIT
*
*     track finder initialisation; called by UGINIT
*
*     Kristian Harder, September 1999
*
**********************************************************************
      IMPLICIT NONE
#include "include/tkffread.inc"
#include "include/brcntl.inc"
#include "include/fkdebug.inc"
#include "include/fsfmc.inc"
#include "include/tpcdbg.inc"
#include "include/tknoisepar.inc"

      INTEGER I

C.....init variables that are to be set by FFREAD
      DO I=1,3
        IPATRC(I)=1     ! default is realistic pattern recognition
      END DO
      MCMERGE=.FALSE.
      FSFDBG=0
      FSFHIST=0
      IDEB=0
      IHIS=0
      TKSTDBG=0
      TKSTHST=0
      DELADBG=0
      TPCDBG=0
      PFCOST=1.0
      PFPCUT=0.0
      PFKINE=.FALSE.
      PFMINH=1
      ITKPERF=0

      IFFORM=0
      IFBEST=0
      IFHIST=0
CSJA set the background to off (in brahms default is on)
      TNNBX=0
      TNENERGY=500

c      CALL FFKEY('PATREC',  IPATRC,3,'INTEGER')
c      CALL FFKEY('MCMERGE', MCMERGE,1,'LOGICAL')
c      CALL FFKEY('DESDBG', FSFDBG,1,'INTEGER')
c      CALL FFKEY('DESHST', FSFHIST,1,'INTEGER')
c      CALL FFKEY('DEFDBG', IDEB,1,'INTEGER')
c      CALL FFKEY('DEFHST', IHIS,1,'INTEGER')
c      CALL FFKEY('TKTDBG', TKSTDBG,1,'INTEGER')
c      CALL FFKEY('TKTHST', TKSTHST,1,'INTEGER')
c      CALL FFKEY('TPCDBG', TPCDBG,1,'INTEGER')
c      CALL FFKEY('DEADBG', DELADBG,1,'INTEGER')
c      CALL FFKEY('PFCOST',PFCOST,1,'REAL')
c      CALL FFKEY('PFPCUT',PFPCUT,1,'REAL')
c      CALL FFKEY('PFKINE',PFKINE,1,'LOGICAL')
c      CALL FFKEY('PFMINH',PFMINH,1,'INTEGER')
c      CALL FFKEY('TKPERF',ITKPERF,1,'INTEGER')
c      CALL FFKEY('BEST',IFBEST,1,'INTEGER')
c      CALL FFKEY('HIST',IFHIST,1,'INTEGER')
c      CALL FFKEY('NBX', TNNBX, 1, 'INTEGER')
c      CALL FFKEY('ENERGY', TNENERGY, 1, 'INTEGER')

C.....initialisation routine for GAB's TPC pattern recognition
      CALL TPCINI

C.....initialisation of RJH's silicon pattern recognition
c      CALL SIPINI

C.....create directory for tracking histograms


      RETURN
      END


