*  
* $Id: fk3trk.F,v 1.4 2007-04-20 13:40:03 rasp Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.3  2006/11/02 11:49:00  rasp
*
* Modified files adopted to FullLDCTracking processor
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:18  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
************************************************************************
*                                                                      *
      SUBROUTINE FK3TRK(IOUTR,NDAT,
     $                  IDTER,DATTE,DATREF,DATTK1,NOUT,
     $                  IDOUT,IER1,IER2,IFLTE,IERX)
*                                                                      *
*   Steering routine to fit one track from input of TER data           *
*   R McNulty and J Wickens (after FK1TRK of P Billoir)                *
*                                                                      *
*   First release : 22.09.95                                           *
*   updated M.Feindt, R. Ehret, R. Fruehwirth   Oct/Nov 95             *
*                                                                      *
*   Input  :
*            IOUTR : steering for outlier logic :                      *
*                    last decimal digit:                               *
*                                  0  : skip outlier logic             *
*                                  1  : do not throw out ITC           *
*                                  2  : do not throw out TPC           *
*                                  3  : allow anything to be thrown out*
*                                  5  : keep either TPC or ITC         *
*                                  6  : allow anything to be thrown out*
*                                         if rank not decreased        *
*                                  8  : do not throw out FTD           *
*                    second last decimal digit:                        *
*                                  0  : called from search             *
*                                  i  : called from TKFB, i=iteration  *
*                                                            step      *
*                    third last decimal digit: 0 usually               *
*                                              1 chi2 definition       *
*                                                 without internal chi2*
*                    fourth last decimal digit: 0 usually              *
*                                               1 no extra/interpolations
*                                                                      *
*                                                                      *
*            NDAT  : Number of TER data input                          *
*            IDTER : Tanagra identifiers of the TER .Dimension NTEMX=15*
*                    Needed for communication with VD transform        *
*            DATTE : Array with TER data .                             *
*                    Must be dimensioned to DATTE(MXTE,NTEMX)  - at    *
*                    present (50,15). The data for the n th. TER should*
*                    be stored in DATTE(1 -> number of words in TER,n) *
*            DATREF: Parameters of the reference trajectory as in a TSR*
*                    Dimension is MXTS = 37                            *
*   Output : DATTK1: TKR data at beginning of the track.               *
*                    Dimension is MXTK = 37                            *
*            IER1  : error flag ( 0 if no error )                      *
*            NOUT  : the number of TER dropped from the fit            *
*            IDOUT : the Tanagra identifiers of the rejected TER       *
*                    Dimension is NTEMX=15                             *
*            IFLTE : flag for each input TER , Not 0 means TER dropped *
*                    from fit . Dimension is NTEMX = 15                *
*                                                                      *
************************************************************************
#include "include/fksave.inc"
#include "include/fkdebug.inc"
*
#include "include/fkclun.inc"
#include "include/fkrank.inc"
#include "include/fkfild.inc"
#include "include/fkcons.inc"
#include "include/fkparm.inc"
#include "include/fktkx.inc"
#include "include/fkcom.inc"
#include "include/fkexts.inc"
#include "include/fkpipar.inc"
      PARAMETER (RTOD=180./PI)
*
      DIMENSION DATTE(MXTE,NTEMX),IFLTE(NTEMX),DATREF(MXTS)
      DIMENSION DATS(MXTS),IDATS(MXTS),DATK(MXTK),IDATK(MXTK),
     &          DATE(MXTE),IDATE(MXTE),IDOUT(NTEMX)
      EQUIVALENCE (DATS,IDATS),(DATK,IDATK),(DATE,IDATE)
      DIMENSION DATTK1(MXTK),IDTER(NTEMX),IERX(NEXSMX)
*
      CHARACTER*30 ERRMES(3)
ccc      EXTERNAL FXPROB        ! MVL
*
      DATA ERRMES /'inconsistent TER data         ',
     &             'error matrix not positive     ',
     &             'error in conversion to Delphi '/
*

      IOUT=MOD(IOUTR,10)
      ITAN=MOD(IOUTR,100)/10
      ICHI=MOD(IOUTR/100,10)
      IEXT=IOUTR/1000
*

CSAJWRT
c      do i=1, MXTE
c      write(*,*) "DATTE(",i,1,") = ", DATTE(i,1) 
c      end do 


      IF(IDEB.GE.1) WRITE(IPRNT,
     &   '(/,'' *** Entering FK3TRK with IOUT= '',I8,/)') IOUT
      IF(IDEB.GE.1) WRITE(IPRNT,
     &   '(/,I3,'' TEs: ''/8(Z8,2X)/7(Z8,2X))')
     &   NDAT,(IDTER(III),III=1,NDAT)

CSAJWRT
c      WRITE(IPRNT,
c     &   '(/,'' *** Entering FK3TRK with IOUT= '',I8,/)') IOUT
CSAJWRT
c      WRITE(IPRNT,
c     &   '(/,I3,'' TEs: ''/8(Z8,2X)/7(Z8,2X))')
c     &   NDAT,(IDTER(III),III=1,NDAT)

*
      IERR=0
      IER1=0
*
      NOUT=0
      CALL VZERO(IFLTE,NTEMX)
      CALL VZERO(IDOUT,NTEMX)
      CALL VZERO(IERTE,NTEMX)
      CALL VZERO(IFLAG,NTEMX)
      CALL VZERO(DATTK1,MXTK)
*
      CALL UCOPY(DATREF,DATS,MXTS)
*
*     Transfer data for reference trajectory
*
      CALL FKSREF(DATS)
*
      IF(IDEB.GE.2) THEN
         WRITE(IPRNT,'(/,''    Reference trajectory '','//
     &       'F10.3,F10.5,F10.3,2F10.5,F10.6,/)') PARREF
      ENDIF
*
*     For each TER, store useful data
*
      IRNKXY=0
      IRNKZ=0
      NTE=NDAT
*
      DO 5 ITE=1,NTE
C(KH) ---begin
        IF (IDEB.GE.2) THEN
CSJAWRT
        WRITE(6,*) 'FK3TRK: TE reference coordinates ',
     >               DATTE(10,ITE),DATTE(11,ITE),DATTE(12,ITE),
     >               DATTE(13,ITE),DATTE(14,ITE)
        ENDIF
C(KH) ---end
        IBADCV=0
        CALL UCOPY(DATTE(1,ITE),DATE(1),MXTE)
*
*       store length of TE (projected onto xy)
*
        IF(IBADCV.EQ.0)THEN
          IF(ABS(DATE(13)).GT..001) THEN
            ALRFTE(ITE)=DATE(9)*SIN(DATE(13))
          ELSE
            ALRFTE(ITE)=DATE(9)
          ENDIF
*
*         store parameters and weight matrix from this TER
*
          CALL FKSTER(DATE,IDET(ITE),TYPE(ITE),MEASUR(ITE),NDTER(ITE),
     &    CHTER(ITE),PARAM(1,ITE),WGTER(1,ITE),IERTE(ITE),IRNKXY,IRNKZ)
*
        ELSE
          IERTE(ITE)=IBADCV
          IDET(ITE)=IDATE(1)
        ENDIF
*

        IF(IERTE(ITE).NE.0) THEN
          WRITE(IPRNT,'('' FKSTER(TKF)FK3TRK :TER'',Z8,''  det.'','//
     &         'I4,2X,A30)') IDTER(ITE),IDET(ITE),ERRMES(IERTE(ITE))
        ENDIF
*
    5 CONTINUE
*
*     do not try to fit if insufficient information in XY projection
      IF(IRNKXY.LT.3) THEN
        IF(IDEB.GE.2) THEN
          WRITE(IPRNT,
     &   '(/,'' FK3TRK(TKF): XY info insufficient : rank ='',I2)')
     &   IRNKXY
        ENDIF
*
*     mif putintocradle! was IER1=1
*
         IER1=2
         IER2=2
         GO TO 999
*
*     if sufficient information in XY projection, but not in Z :
*     increase slightly weights on z and theta for each TE
      ELSE IF(IRNKZ.LT.2) THEN
        DO 6 ITE=1,NTE
          WGTER(3,ITE)=WGTER(3,ITE)+.01
          WGTER(6,ITE)=WGTER(6,ITE)+1.
    6   CONTINUE
*
      ENDIF
*
*     Perform the fit
*
      CALL FKTFIT(IDEB,IHIS,IOUTR,IRNKXY,IRNKZ,
     & NTEFIT,IPATTB,IPATTF,ISREFB,ISREFF,IER1,IER2,IERX,IDTER)
*
*     Fill TKR array with the 'TS' data plus result of the backward fit
*     (if chi2 is acceptable)
*
      IF((IER1.EQ.1.AND.IRNKXY.LT.3).OR.IER1.GT.1) THEN
        IF(IDEB.GE.1) WRITE(IPRNT,'(''   IERB ='',I3,'' in FKTFIT'')')
     &                      IER1
      ELSE
*
*       Cut on chi2/ndeg
*
c MVL
c MVL      IF(CHBACK.GT.0.0.and.CHBACK.lT.300.0.and.NDBACK.GT.0)THEN
c MVL        PPP=max(PROB(CHBACK,NDBACK),1.e-25)
c MVL      ELSE
c MVL        write(*,*) ' Cut --- reason1',CHBACK,NDBACK
c MVL        PPP=1.e-25
c MVL      ENDIF

        IF(IOUT.NE.0.AND.NDBACK.GT.0 .AND.
c MVL     +        FXPROB(CHBACK,NDBACK).LT.PRTOT(IBAFO)) THEN
     +           PROB(CHBACK,NDBACK).LT.PRTOT(IBAFO)) THEN
          IER1=4
          IF (IDEB.GE.1)  WRITE (IPRNT,1234) PROB(CHBACK,NDBACK)
 1234     FORMAT(1X,'BACKWARD FIT PROB CUT ',E12.4)
        ELSE
          CALL UCOPY(DATS,DATK,16)
          IDATK(2)=2*IPATTB
          IF(TYPSUR(ISREFB).EQ.'CYLI') THEN
            IDATK(3)=1
          ELSE
            IDATK(3)=0
          ENDIF
          IF(IRNKZ.LT.2) IDATK(3)=IDATK(3)+2
          IDATK(5)=NTEFIT
          IDATK(8)=NDBACK
          DATK(9)=CHBACK
          CALL FKSTKR(ISREFB,-1,IHIS,DATK(17),IERS)
          IER1=MAX(IER1,IERS)
*
*         TK length = distance between ref. points of first and las
*                     TEs plus the length of the last TE
*
          DATK(10)=(ALRPHI(ISREFF)-ALRPHI(ISREFB))
          INUM=INUMB(ISREFF)
          IF(INUM.GT.0) DATK(10)=DATK(10)+ALRFTE(INUM)
          IF(ABS(DATK(20)).GT..001.AND.
     >       ABS(DATK(20)).LT.pi-0.01) DATK(10)=DATK(10)/SIN(DATK(20))
          IF(IER1.NE.0.AND.IDEB.GE.1)
     +        WRITE(IPRNT,'(A30,'' in FKSTKR'')') ERRMES(2)
          IF (IER1.LE.1) THEN
*
*           start point is redefined as that fitted on reference surface
            IF(TYPSUR(ISREFB).EQ.'CYLI') THEN
              CAPPHI=DATK(18)/DATK(17)
              DATK(11)=DATK(17)*COS(CAPPHI)
              DATK(12)=DATK(17)*SIN(CAPPHI)
              DATK(13)=DATK(19)
            ELSE
              CALL UCOPY(DATK(17),DATK(11),3)
            ENDIF
*
*           compute radius and centre of curvature, and coordinates
*           of last point with the length previously defined
*
            RTRK=SIN(DATK(20))/(CONSB*DATK(22))
            RDPHI=DATK(10)*SIN(DATK(20))
            PHIEND=DATK(21)+RDPHI/RTRK
            DATK(14)=DATK(11)+RTRK*(SIN(PHIEND)-SIN(DATK(21)))
            DATK(15)=DATK(12)-RTRK*(COS(PHIEND)-COS(DATK(21)))
            IF(ABS(DATK(20)).GT..001.AND.
     >         ABS(DATK(20)).LT.pi-0.01) THEN
              DATK(16)=DATK(13)+RDPHI/TAN(DATK(20))
            ELSE
              DATK(16)=DATK(13)+RDPHI
            ENDIF
*
*           the charge is the "preferred" one,
*           (unknown charge is not allowed)
            IF(-CONSB*DATK(22).GE.0.) THEN
              IDATK(6)=1
            ELSE
              IDATK(6)=2
            ENDIF
*
*           debug output
*
            IF(IDEB.GE.2) THEN
              WRITE(IPRNT,
     $             '(/,''   Fitted parameters at first point'','//
     $             '2X,A4,F10.3,F10.5,F10.3,2F10.5,F12.7)')
     $             TYPSUR(ISREFB),PARAMS(1,ISREFB),
     $             (PARAMS(K+1,ISREFB)+DQBCK2(K,INUMB(ISREFB)),K=1,5)
            ENDIF
            IF(IDEB.GE.1) WRITE(IPRNT,
     $          '(/,''   Successful backward fit - TKR data : '','//
     $         'I5,2(2X,Z8),5I5,/,8G14.6,/,6G14.6,/,10G12.4,/,5G12.4)')
     $          (IDATK(K),K=1,8),(DATK(K),K=9,37)
ccc            IF(IDEB.GT.1.AND.ITAN.GT.1) CALL FKCVPR(IPRNT,DATK(23))
*
*           Save TKR data in array DATTK1
*
            CALL UCOPY(DATK,DATTK1,MXTK)
          ENDIF
        ENDIF
      ENDIF
*
*     Computation of interpolations and extrapolations. Weighted means
*     on extrapolation surfaces if both backward and forward fits
*     succeeded. If only backward take these values.
*
*     IDEB = 3

      IF(IER1.EQ.0.AND.IEXT.EQ.0) THEN
*
        IF(IDEB.GE.2) WRITE(IPRNT,
     &  '(/,'' Computation of inter/extrapolations ---------'',/)')
*
        DO 10 IX=1,NEXS
        IF(IERX(IX).NE.0) THEN
          IF(IDEB.GE.2) WRITE(IPRNT,
     &    '(''    extr. surf.'',I4,2X,A4,''  IERX ='',I3)')
     &    IX,TYPX(IX),IERX(IX)
        ELSE
          IF(IER2.EQ.0)THEN
            CALL FKEXMN(DQBCKE(1,IX),WGBCKE(1,IX),DQFRWE(1,IX),
     &            WGFRWE(1,IX),DQX(1,IX),WGTX(1,IX),COVX(1,IX),IERE) ! this line was commented in the original version
            IF(IERE.NE.0) THEN
              IERX(IX)=-9999
              IF(IDEB.GT.0) WRITE(IPRNT,
     &        '(''    extr. surf.'',I4,2X,A4,'' : error in FKEXMN'')')
     &        IX,TYPX(IX)
            ELSE
               CALL FKDPVA(PARAMX(2,IX),DQX(1,IX),PARAMX(2,IX),5) ! this line was commented in the original version
            ENDIF
          ELSE
c            The following two lines were commented
             CALL FKDPVA(PARAMX(2,IX),DQBCKE(1,IX),PARAMX(2,IX),5) ! this line was commented in the original version
             CALL FKMI5D(WGBCKE(1,IX),COVX(1,IX),IERE) ! this line was commented in the original version
            IF(IERE.NE.0) THEN
              IERX(IX)=-9999
              IF(IDEB.GT.0) WRITE(IPRNT,
     &        '(''    extr. surf.'',I4,2X,A4,'' : error in FKMI5D'')')
     &        IX,TYPX(IX)
            ENDIF
          ENDIF
          IF(IDEB.GE.2) THEN
            WRITE(IPRNT,
     $           '(''    extr. surf.'',I4,2X,A4,'' fit  par.'','//
     $           'T39,3F11.4,2F11.5,F11.7)')
     $           IX,TYPX(IX),(PARAMX(J,IX),J=1,6)
            WRITE(IPRNT,
     $           '(30X,''errors'',T39,11X,2F11.6,2F11.5,F11.7)')
     $           DSQRT(COVX(1,IX)),DSQRT(COVX(3,IX)),DSQRT(COVX(6,IX)),
     $           DSQRT(COVX(10,IX)),DSQRT(COVX(15,IX))
          ENDIF
        ENDIF
   10   CONTINUE
*
*   transform into TANAGRA format
*
        DO 11 IX=1,NEXS
            IFLX(IX) = IERX(IX)
            IF(IERX(IX).EQ.0) THEN
              CALL FKSPAR(TYPX(IX),PARAMX(1,IX),PARAMX(2,IX),COVX(1,IX),
     &             DATX(1,IX)) ! this line was commented in the original version
C              CALL VZERO(DATX(22,IX),3) ! this line was uncommented in the original version
              IF(IDEB.GE.2) WRITE(IPRNT,
     $             '(/,6X,''TKX data ('',A4,'') : '','//
     $             '6G14.6,/,8X,10G12.4,/,8X,5G12.4)')
     $             TYPX(IX),(DATX(K,IX),K=1,21)
            ENDIF
            IF(IERX(IX).EQ.-9999)IERX(IX)=0
   11   CONTINUE
      ENDIF

      IDEB = 0
*
*           build the list of rejected TER's
*
            DO 7 ITE=1,NTE
              IF(IERTE(ITE).NE.0.OR.IFLAG(ITE).NE.0) THEN
*               NOTE: max outlayers are limited by Delana dimensions
                IF(NOUT.LT.NTEMX)THEN
                  NOUT=NOUT+1
                  IDOUT(NOUT)=IDTER(ITE)
                  IFLTE(ITE)=10*IERTE(ITE)+IFLAG(ITE)
                ELSE
                  WRITE(IPRNT,'(1X,''FK3TRK(TKF) : OL OVERFLOW'')')
                  CALL TKVETO(1)
                  IER1=5
                  IER2=5
                ENDIF
              ENDIF
    7       CONTINUE
*
  999 CONTINUE
*
      END
