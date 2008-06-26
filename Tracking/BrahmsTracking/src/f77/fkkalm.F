*     
*     $Id: fkkalm.F,v 1.2 2008-06-26 09:35:56 aplin Exp $
*     
*     $Log: not supported by cvs2svn $
*     Revision 1.1.1.1  2005/07/11 15:15:48  aplin
*     Initial version
*     
*     Revision 1.1.1.1  2005/06/18 10:05:34  gaede
*     Initial local version
*     
*     Revision 1.1  2005/06/13 16:14:18  aplin
*     *** empty log message ***
*     
*     Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
*     Brahms V308
*     
*     
************************************************************************
*     *
      SUBROUTINE FKKALM(IDIR,IESURF,IDEB,IHIS,
     &     NTEFIT,IPATT,ISUREF,IERR)
*     *
*     *
*     Input  :    IDIR  :  -1 for backward filter , 1 for forward filter *
*     IESURF  :  serial number of first surface in the filter  *
*     IDEB  :  level of debug output                         *
*     IHIS  :  flag for histograms                           *
*     *
*     Output :  NTEFIT  :  number of TE's included in the fit            *
*     IPATT  :  bit pattern of detectors used in the fit      *
*     ISUREF  :  reference surface for parameters              *
*     IERR  :  error flag :  0 if OK                         *
*     1 if insufficient information   *
*     2 if weight matrix not positive *
*     3 if abnormal topology          *
*     *
************************************************************************
*     
#include "include/fksave.inc"
#include "include/fkclun.inc"
#include "include/fkparm.inc"
#include "include/fkcom.inc"
#include "include/fkcons.inc"

*     
      DOUBLE PRECISION PI,TWOPI,DQKALM(5),WGKALM(15),PARSTA(6),
     /     PAREND(6),DER(8)
*     
      CHARACTER*4 TYPFIT
      PARAMETER (PI=3.14159265359D00,TWOPI=PI*2.0D00)
*     
      
      IERR=0
      IRANK=0
      IRNKRF=0
      IRNKZT=0
      ISUREF=0
      NTEFIT=0
      IPATT=0
      CALL FKDPV0(DQKALM,5)
      CALL FKDPV0(WGKALM,15)
      NDEG=0
      CHI2=0.
*     
*     loop from first surface to last one
*     
      TYPFIT=TYPE(INUMB(IORDS(IESURF)))
      IF(IDIR.EQ.1) THEN
         NSF=IESURF
         NSL=NSURF
      ELSE
         NSF=IESURF
         NSL=1
      ENDIF
*     
      DO 10 IS=NSF,NSL,IDIR
*     
         ISURF=IORDS(IS)
         INUM=INUMB(ISURF)
*     
*     if measurement, add information (weighted mean) and update chi2

*     
         IF(INUM.GT.0.AND.IFLAG(INUM).EQ.0) THEN
*     
*     if first information found : copy deviations, weight matrix, rank and chi2
*     
               IF(IRANK.EQ.0) THEN
*     
                  IF(IDEB.GE.2) WRITE
     &                 (IPRNT,'(''       First TE : det'','//
     &                 'I4,''  meas. code '',Z8)') 
     &                 IDET(INUM),MEASUR(INUM)
*     
                  CALL FKDPCO(DQTER(1,INUM),DQKALM,5)
                  CALL FKDPCO(WGTER(1,INUM),WGKALM,15)
                  NDEG=NDTER(INUM)
                  CHI2=CHTER(INUM)
*     
                  IRNKRF=JBIT(MEASUR(INUM),1)+JBIT(MEASUR(INUM),4)+
     &                 JBIT(MEASUR(INUM),5)
                  IRNKZT=JBIT(MEASUR(INUM),2)+JBIT(MEASUR(INUM),3)
                  IRANK=IRNKRF+IRNKZT
*     
                  IF(IDEB.GE.2) then
                     WRITE(IPRNT,*)
     &                    '(''       after 1st  : rank'',I2,''+'//
     &                    ',I1,6X,''''deviat.'',T50,5F11.6,'//
     &                    '  nd, chi2 '',I4,1X,G10.4)',
     &                    IRNKRF,IRNKZT,DQKALM,NDEG,CHI2

                     write(IPRNT,
     $                    '(a,t10,e17.8,/t10,2e17.8,/t10,3e17.8,'//
     $                    '/t10,4e17.8,/t10,5e17.8)') 'COV : ',wgkalm
                  end if
*     
*     if new information : add it to the fit
*     
               ELSE
*     
                  IF(IDEB.GE.2) 
     &                 WRITE(IPRNT,'(''       Next TE : det'',I4,'//
     &                 '''  meas. code '',Z8)') IDET(INUM),MEASUR(INUM)
*     
*     store deviations, weight matrix and ranks before adding information
*     
                  IF(IDIR.EQ.1) THEN
                     CALL FKDPCO(DQKALM,DQFRW1(1,INUM),5)
                     CALL FKDPCO(WGKALM,WGFRW1(1,INUM),15)
                     NDFRW1(INUM)=NDEG
                     CHFRW1(INUM)=CHI2
                     IRFFRW(INUM)=IRNKRF
                     IRZFRW(INUM)=IRNKZT
                  ELSE
                     CALL FKDPCO(DQKALM,DQBCK1(1,INUM),5)
                     CALL FKDPCO(WGKALM,WGBCK1(1,INUM),15)
                     NDBCK1(INUM)=NDEG
                     CHBCK1(INUM)=CHI2
                     IRFBCK(INUM)=IRNKRF
                     IRZBCK(INUM)=IRNKZT
                  ENDIF
*     
*     update the rank of information on  r*Phi,phi,1/R  and on  z,theta
*     
                  JRNKRF=JBIT(MEASUR(INUM),1)+JBIT(MEASUR(INUM),4)+
     &                 JBIT(MEASUR(INUM),5)
                  JRNKZT=JBIT(MEASUR(INUM),2)+JBIT(MEASUR(INUM),3)
                  JRANK=JRNKRF+JRNKZT
                  KRNKRF=MIN0(IRNKRF+JRNKRF,3)
                  KRNKZT=MIN0(IRNKZT+JRNKZT,2)
                  KRANK=KRNKRF+KRNKZT
*     
*     partial information : provisional solution : add small contribution to
*     diagonal elements, then compute weighted mean
*     
                  IF(KRANK.LT.5) THEN
*     
                     IF(IDEB.GE.2) 
     &                    WRITE(IPRNT,*)
     &                    '(''       Incomplete information '//
     &                    'in first TEs : '',''''increase '  //
     &                    'slightly the diagonal weights'')'
*     
                     IF(KRNKRF.LE.0) WGKALM( 1)=WGKALM( 1)+1.
                     IF(KRNKZT.LE.0) THEN
                        IF(TYPFIT.EQ.'CYLI') THEN
                           WGKALM( 3)=WGKALM( 3)+0.01
                        ELSE
                           WGKALM( 3)=WGKALM( 3)+1.
                        ENDIF
                     ENDIF
                     IF(KRNKZT.LE.1) WGKALM( 6)=WGKALM( 6)+1.
                     IF(KRNKRF.LE.1) WGKALM(10)=WGKALM(10)+1.
                     IF(KRNKRF.LE.2) WGKALM(15)=WGKALM(15)+100.
                  ENDIF
*     
*     compute weighted mean
*     
                  CALL FKWGMN(KRANK,DQKALM,WGKALM,NDEG,CHI2,IRANK,
     &                 DQTER(1,INUM),WGTER(1,INUM),NDTER(INUM),
     &                 CHTER(INUM),JRANK,DQKALM,WGKALM,NDEG,CHI2,IERR)

                  IF(KRANK.LT.5) THEN
*     
                     IF(IDEB.GE.2) WRITE(IPRNT,
     &                    '(''       Incomplete information in '//
     &                    'first TEs : '',' //
     &                    '''subtract again the diagonal weights'')')
*     
                     IF(KRNKRF.LE.0) WGKALM( 1)=WGKALM( 1)-1.
                     IF(KRNKZT.LE.0) THEN
                        IF(TYPFIT.EQ.'CYLI') THEN
                           WGKALM( 3)=WGKALM( 3)-0.01
                        ELSE
                           WGKALM( 3)=WGKALM( 3)-1.
                        ENDIF
                     ENDIF
                     IF(KRNKZT.LE.1) WGKALM( 6)=WGKALM( 6)-1.
                     IF(KRNKRF.LE.1) WGKALM(10)=WGKALM(10)-1.
                     IF(KRNKRF.LE.2) WGKALM(15)=WGKALM(15)-100.
                  ENDIF
*     
                  IF(IERR.EQ.0) THEN
                     IRNKRF=KRNKRF
                     IRNKZT=KRNKZT
                     IRANK=KRANK
*     
                     IF(IDEB.GE.2) WRITE(IPRNT,*)
     &                    '(''       after add. : rank'',I2,'//
     &                    '+'',I1,6X,'//
     &                    '''deviat.'',T50,5F11.6,''   nd, chi2 '//
     &                    ',I4,1X,G10.4)',
     &                    IRNKRF,IRNKZT,DQKALM,NDEG,CHI2
                  ELSE
                     IF(IDEB.GE.1) THEN
                        WRITE(IPRNT,*) ' FKKALM(TKF) : error in FKWGMN'
                        WRITE(IPRNT,'(''       FKWGMN failed'')')
                     END IF
                     RETURN
                  ENDIF
               ENDIF
               if (ideb.ge.2) then
                  write(IPRNT,
     $                 '(a,t10,e17.8,/t10,2e17.8,/t10,3e17.8,'//
     $                 '/t10,4e17.8,/t10,5e17.8)') 'COV : ',wgkalm
               end if
*     
               NTEFIT=NTEFIT+1
               ISUREF=ISURF
               CALL SBIT1(IPATT,IDET(INUM))
*     
*     store deviations and weight matrix after adding information
*     
               IF(IDIR.EQ.1) THEN
                  CALL FKDPCO(DQKALM,DQFRW2(1,INUM),5)
                  CALL FKDPCO(WGKALM,WGFRW2(1,INUM),15)
                  NDFRW2(INUM)=NDEG
                  CHFRW2(INUM)=CHI2
               ELSE
                  CALL FKDPCO(DQKALM,DQBCK2(1,INUM),5)
                  CALL FKDPCO(WGKALM,WGBCK2(1,INUM),15)
                  NDBCK2(INUM)=NDEG
                  CHBCK2(INUM)=CHI2
               ENDIF
c            ENDIF
*     
*     if interpolation/extrapolation requested, store deviations from
*     reference trajectory and weight matrix
*     
         ELSE IF(INUM.LT.0) THEN
            IF(IDIR.EQ.-1) THEN
               CALL FKDPCO(DQKALM,DQBCKE(1,-INUM),5)
               CALL FKDPCO(WGKALM,WGBCKE(1,-INUM),15)
            ELSE
               CALL FKDPCO(DQKALM,DQFRWE(1,-INUM),5)
               CALL FKDPCO(WGKALM,WGFRWE(1,-INUM),15)
            ENDIF
*     
*     if material surface, add contribution of multiple scattering
*     
         ELSE IF(INUM.EQ.0)THEN
            IF(IRANK.NE.0) THEN
               CALL FKMSCT(SINTH,COSTH,PINV,XRLS(ISURF),WGKALM)
            ENDIF
         end if
*     
*     propagate to next surface (unless the first one is reached)
*     
         IF(IS.EQ.NSL) GO TO 10
         INEXT = IORDS(IS+IDIR)
         CALL FKDPCO(PARAMS(1,ISURF),PARSTA,6)
         CALL FKDPCO(PARAMS(1,INEXT),PAREND,6)

         if (inum.eq.0 .and. idir.eq.1)
     $        call fkdpco(paraft(1,isurf),parsta,6)
         if (inumb(inext).eq.0 .and. idir.eq.-1)
     $        call fkdpco(paraft(1,inext),parend,6)
*     
*     if next surface has a type different from TYPFIT, switch TYPFIT
*     and transform parameters, deviations and weight matrix
*     
         IF(TYPFIT.NE.TYPSUR(INEXT)) THEN
            CALL FKSWCH(TYPFIT,PARSTA,DQKALM,WGKALM)
            IF(IDEB.GE.2) then
               WRITE(IPRNT,
     &              '(7X,A4,'' to '',A4,'' : deviations'',T50,5F11.6)')
     &              TYPFIT,TYPSUR(INEXT),DQKALM
               write(IPRNT,
     $              '(a,t10,e17.8,/t10,2e17.8,/t10,3e17.8,'//
     $              '/t10,4e17.8,/t10,5e17.8)') 'COV : ',wgkalm
            end if
            TYPFIT=TYPSUR(INEXT)
         ENDIF
*     
*     compute derivatives and propagate deviations
*     then compute inverse derivatives and propagate weight matrix
*     
         IF(TYPFIT.EQ.'CYLI') THEN
            CALL FKXDCY(PARSTA,PAREND,DER)
            CALL FKPRCD(DQKALM,DER,DQKALM)
            CALL FKXDCY(PAREND,PARSTA,DER)
            CALL FKPRCW(WGKALM,DER,WGKALM)
         ELSE
            CALL FKXDPL(PARSTA,PAREND,DER)
            CALL FKPRPD(DQKALM,DER,DQKALM)
            CALL FKXDPL(PAREND,PARSTA,DER)
            CALL FKPRPW(WGKALM,DER,WGKALM)
         ENDIF
*     
         IF(IDEB.GE.2) then
            WRITE(IPRNT,
     &           '(''       propag. to'',8X,F8.2,'' : deviat.'',T50,'//
     &           '5F11.6)') PAREND(1),DQKALM
            write(IPRNT,
     $           '(a,t10,e17.8,/t10,2e17.8,/t10,3e17.8,'//
     $           '/t10,4e17.8,/t10,5e17.8)') 'COV : ',wgkalm
         end if
*     
 10   CONTINUE
*     
      IF(IDIR.EQ.-1) THEN
         NDBACK=NDEG
         CHBACK=CHI2
         CALL FKDPCO(DQKALM,DQBACK,5)
         CALL FKDPCO(WGKALM,WGBACK,15)
      ELSE
         NDFORW=NDEG
         CHFORW=CHI2
         CALL FKDPCO(DQKALM,DQFORW,5)
         CALL FKDPCO(WGKALM,WGFORW,15)
      ENDIF
*     
*     error exit if the rank at the end is not 5
*     
      IF(IRANK.LT.5) THEN
         IERR=1
         IF(IDEB.GE.1) WRITE(IPRNT,
     &        '(''       Insufficient information in the TEs'')')
      ENDIF
*     
      IF(IDEB.GE.1) WRITE(IPRNT,
     &     '(7X,''Kalman Filter end  -  NDEG, CHI2 :'',I4,F12.2)') 
     &     NDEG,CHI2
*     
      END

