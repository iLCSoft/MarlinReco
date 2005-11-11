*  
* $Id: tpcini.F,v 1.2 2005-08-08 07:09:13 aplin Exp $
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
      SUBROUTINE TPCINI
C****************************************************************************
c* Version number: 0.1      Date of last modification   6.11.99             *
c*------------------------------------------------------------------------- *
c* Set up the pattern recognition constants from the cards file             *
c* List of modifications since last version:                                *
C* Main Authors:                                                            *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
      IMPLICIT NONE

#include "include/fndparms.inc"
#include "include/strcons.inc"
#include "include/tphtcns.inc"
#include "include/patrec.inc"
#include "include/trfdct.inc"
#include "include/rflags.inc"
#include "include/tchajj.inc"
#include "include/tpgeom.inc"
c#include "marlin_tpcgeom.h"
#include "include/tpgeop.inc"
#include "include/tpcdims.inc"
#include "include/my_ttrc.inc"

      CALL VZERO(IPRFLAG,10)
c      CALL FFKEY('PRFL',IPRFLAG,10,'INTEGER')

C  D(CHI**2)/KINK
      DCM=150.0
c      call ffkey('KNDCHI',DCM, 1, 'real')

C  MIN. POINTS FOR KINK
      MKINK=10
c      call ffkey('MINKINK',MKINK, 1, 'integer')

C MAX. POINTS TO THROW AWAY
      LOMAX=5
c      call ffkey('MAXTHROW',LOMAX, 1, 'integer')

C  MINIMAL NUMBER OF POINTS ON A TRACK
      MINPT=4
c      call ffkey('MNPN',MINPT, 1, 'integer')

C  MIN. FIT PROBABILITY FOR TRACK
      PMIN=0.00030
c      call ffkey('PRBM',PMIN, 1, 'real')

c      DATA CHIPM/300.0/
c    CHI2 DIFF OF LINKS
      CHIPM=300.0
c      call ffkey('MINCHI', CHIPM, 1, 'real')

C TPC occupancy:
      TPC_OCC=0.01
c      call ffkey('TPCOCC', TPC_OCC, 1, 'real')

C Set the threshold and saturation values for TPC hit coordinate inclusions
c the numbers are in ionizition electron number equivalence.
      TPSAT=999.
      TPTHRS=1.
c      call ffkey('TPCSAT', TPSAT, 1, 'real')
c      call ffkey('TPCTHR', TPTHRS, 1, 'real')

C     Max chisq for link helix (dz) test
      TDZTOL=16.
c      call ffkey('DCHIZ', TDZTOL, 1, 'real')

C Max value of Z0 for accepting link
      TZ0TOL=400.0
c      call ffkey('DZTOL', TDZTOL, 1, 'real')
C Max value of D0 for accepting link
      TD0TOL=160.
c      call ffkey('DZ0MX', TD0TOL, 1, 'real')
C Min radius of curvature for accepting link
      RHOMIN=10.
c      call ffkey('MINRHO', RHOMIN, 1, 'real')
C Max r-phi window for finding hits
      TDWMAX=1.
c      call ffkey('RPHWINDW', TDWMAX, 1, 'real')
C Max z window for finding hits
      TDZMAX=1.
c      call ffkey('ZWINDW', TDZMAX, 1, 'real')
C Nominal r-phi road width for finding hits
      TDWSTR=5.
c      call ffkey('RPRD', TDWSTR, 1, 'real')
C Nominal z road width for finding hits
      TDZSTR=5.
c      call ffkey('ZRDW', TDZSTR, 1, 'real')
C Max no of successive padrows with no found hit for local search
      MAXGAP= 3
c      call ffkey('MAXGAP', MAXGAP, 1, 'integer')
C Max difference in z between successive padrows for accepting hit in link
      DZNEXT= 45.
c      call ffkey('DZMX', DZNEXT, 1, 'real')
C Max diff. in phi between successive padrows for accepting hit in link
      DPNEXT= 0.3
c      call ffkey('DPMX', DPNEXT, 1, 'real')
C Max chisq of curvature diff. between successive links to accept in chain
      CURTOL= 16.0
c      call ffkey('DCHCURV', CURTOL, 1, 'real')
C Max number of hits that can be rejected in fitting chain
      MAXTRY= 4
c      call ffkey('RJMX', MAXTRY, 1, 'integer')
C Max value of chisq per degree of freedom for r-phi circle fit
      CHLIM1=18.
c      call ffkey('MAXRPCHI', CHLIM1, 1, 'real')
C Max value of chisq per degree of reedom for s-z line fit
      CHLIM2=18.
c      call ffkey('MAXZCHI', CHLIM2, 1, 'real')
C Obsolete variables
      DIPTOL=0.
      SIGMAW=0.

c----
C TTTRC Variables:
C Validity range
      MY_JTTRVR1=1
      MY_JTTRVR2=999999
c      call ffkey('VLD1', MY_JTTRVR1, 1, 'INTEGER')
c      call ffkey('VLD2', MY_JTTRVR2, 1, 'INTEGER')
c Number of attempts to find chains
      MY_JTTRNT=3
c      call ffkey('CHATMPT', MY_JTTRNT, 1, 'INTEGER')
c Number of pad rows a chain may jump over
      MY_JTTRNG=2
c      call ffkey('MAXJUMP', MY_JTTRNG, 1, 'INTEGER')
C minimal number of coordinates for chain
      MY_JTTRNM=4
c      call ffkey('MCHP', MY_JTTRNM, 1, 'INTEGER')
c Window width in phi for triplet search
      MY_JTTRDP=0.35
c      call ffkey('PHITRP', MY_JTTRDP, 1, 'real')
c Window width in z for triplet search
      MY_JTTRDZ=50.
c      call ffkey('ZTRP', MY_JTTRDZ, 1, 'real')
c Road width in (r*phi) for chaining
      MY_JTTRRR=2.
c      call ffkey('RPROADCH', MY_JTTRRR, 1, 'real')
c Road width in z for chaining
      MY_JTTRRZ=6.
c      call ffkey('ZROADCH', MY_JTTRRZ, 1, 'real')

c Max. chi**2/(deg.freedom) for a chain
      MY_JTTRCC=12.
c      call ffkey('MXCH', MY_JTTRCC, 1, 'real')
c Max. num. stand.dev. for up,down links
      MY_JTTRSU=20.
c      call ffkey('MXUD', MY_JTTRSU, 1, 'real')
c Max. num. stand.dev. for spiral links in TPC
      MY_JTTRST=30.
c      call ffkey('MXSP', MY_JTTRST, 1, 'real')
c Max. num. stand.dev. for links through ITC
      MY_JTTRSI=35.
c      call ffkey('MAXSDSP', MY_JTTRST, 1, 'real')


c Max. num. points to be thrown away in fit
      MY_JTTRLM=5
c      call ffkey('TKMP', MY_JTTRLM, 1, 'integer')
c Min. num. points for kink search
      MY_JTTRMK=10
c      call ffkey('MINP', MY_JTTRMK, 1, 'integer')
c Min. delta chi**2 to accept a kink
      MY_JTTRD2=20.
c      call ffkey('KNMN', MY_JTTRD2, 1, 'real')
c Minimal fit probability to accept a track
      MY_JTTRPM=0.001
c      call ffkey('KNPR', MY_JTTRPM, 1, 'real')

C TPC readout pixel size (cm) (from CDR)
       PIX_RP=0.33
c      call ffkey('TPXR', PIX_RP, 1, 'real')
       PIX_Z=1.4
c      call ffkey('TPXZ', PIX_Z, 1, 'real')


  999 RETURN
      END

C----------------------------------------------------------------------

