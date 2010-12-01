*  
* $Id: tfchek.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
        SUBROUTINE TFCHEK(VV0,LTASSO,IERR)
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C------------------------------------------------------------
C
C! Check TPC helix parameters for basic sanity
C
C  Input:       VV0(5)    helix fit(R*4)
C               LTASSO    helix convention(L)
C                         .TRUE. if TASSO convention used
C                         .FALSE. if ALEPH (standard) convention used
C  Output:      IERR      error code (I*4)
C                         0 = no problem
C                  non-zero = fit appears to be garbage
C
C               D. Casper  09/07/97
C
C-----------------------------------------------------------
#include "include/fndparms.inc"
#include "include/alcons.inc"
#include "include/tpcdims.inc"
#include "include/coildims.inc"
C -----------------------------------------------------------------------
      LOGICAL FIRST
C------------------------------------------------------------------------
C GAB:
C        parameter (DMAX = 180., ZMAX=225., PMAX=1000.)
        REAL DMAX,ZMAX,PMAX
        SAVE DMAX,ZMAX,PMAX
        DIMENSION VV0(*)
        LOGICAL LTASSO
        DIMENSION VV1(5)
        CHARACTER*40 REASON(6)
        DATA REASON/'(Field is zero)','(d0 too large)',
     &              '(1/Radius is zero)','(Momentum > 1 TeV)',
     &              '(Points away from Aleph)','(misses in z)'/
C
C  Get into a coordinate system we understand, if necessary
C
C GAB:
        SAVE

        DATA FIRST /.TRUE./
        IF (FIRST)THEN
         FIRST=.FALSE.
         DMAX=TPCACRO
         ZMAX=ZDRIFT
         PMAX=1000.
        ENDIF

        IF (LTASSO) THEN
            VV1(1) = - VV0(1)
            VV1(2) = VV0(2)
            VV1(3) = VV0(3)
            VV1(4) = VV0(4)*SIGN(1.,VV1(1))
            VV1(5) = VV0(5)
        ELSE
            CALL UCOPY(VV0,VV1,5)
        ENDIF
        IERR = 1
C  Did dummy call us with no field?
C GAB:
c        FIELD = ALFIEL(DUM)
        FIELD=BFIELD

        IF (FIELD.EQ.0) GOTO 999
        ROVERP = 1./(0.29979*FIELD/10.)*100.
C  Does the track pass within the TPC radius?
        IERR = 2
        IF (ABS(VV1(4)).GT.DMAX) GOTO 999
C  Does the track have an absurd momentum?
        IERR = 3
        IF (VV1(1).EQ.0.) GOTO 999
        IERR = 4
        PT = 1./(ROVERP*ABS(VV1(1)))
        XLAMBDA = ATAN(VV1(2))
        P = PT/COS(XLAMBDA)
        IF (P.GT.PMAX) GOTO 999
C  Does the track pass within the TPC endplates?
        IF (ABS(VV1(5)).GT.ZMAX) THEN
            IERR = 5
            IF (VV1(2)*VV1(5) .GE. 0.) GOTO 999
            IERR = 6
            RSQEND = VV1(4)**2 + 4.*((1.-VV1(4)*VV1(1))/VV1(1)**2)
     &          *SIN(VV1(1)*(ABS(VV1(5))-ZMAX)/(2.*VV1(2)))**2
            IF (RSQEND.GT.DMAX**2) GOTO 999
        ENDIF
C  Track is probably OK
        IERR = 0
        RETURN
C  Come here on trouble
 999    CONTINUE
c GAB:
c        CALL RERROR('TFCHEK',IERR,
c     &      'Track fails sanity check '//REASON(IERR))
C(KH)   WRITE(*,*) 'Track fails sanity check ',REASON(IERR)
        WRITE(*,*) 'TFCHEK: Track fails sanity check ',REASON(IERR)
        RETURN
        END
C-----------------------------------------------------
