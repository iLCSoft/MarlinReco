*  
* $Id: fkster.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
************************************************************************
*                                                                      *
      SUBROUTINE FKSTER(DATE,IDET,TYPE,MEASUR,ND,CH2,PARAM,WEIGHT,
     &                  IERR,IRNKXY,IRNKZ)
*                                                                      *
*   Decodes TER bank and determines parameters and weight matrix       *
*   in fixed format : 5-vector and (5x5) matrix                        *
*                                                                      *
*   Input  :  DATE  :  TER data array                                  *
*                                                                      *
*   Output :  IDET  :  Identifier of detector                          *
*             TYPE  :  Type of surface ('CYLI' or 'PLAN')              *
*             MEASUR : Measurement code on 5 bits                      *
*             ND    :  Number of degrees of freedom                    *
*             CH2   :  chi2                                            *
*             PARAM :  parameters in "internal" format :               *
*                        r(fixed) ; Phi,z,theta,beta,1/R for 'CYLI'    *
*                        z(fixed) ; x,y,theta,phi,1/R    for 'PLAN'    *
*             WEIGHT(15) : (5x5) weight matrix on parameters           *
*             IERR  :  error flag  :  0 if OK                          *
*                                     1 if inconsistency in TER data   *
*                                     2 if error matrix not positive   *
*                                                                      *
*   Input/Output (incremented if no error) :                           *
*             IRNKXY : rank of XY information (0 to 3)                 *
*             IRNKZ  : rank of Z information  (0 to 2)                 *
*                                                                      *
************************************************************************
*
#include "include/fksave.inc"
#include "include/fkfild.inc"
#include "include/fsfparams.inc"
#include "include/fkcons.inc"
*
      DOUBLE PRECISION PARAM,COVI,WEIGHT
      DIMENSION DATE(MXTE)
*
      DIMENSION PARAM(6),COVI(15),WEIGHT(15)
*
      CHARACTER*4 TYPE
      DIMENSION CHMAX(40),CHMIN(40)
      DATA CHMAX/
     +   6.6349,  9.2103, 11.3417, 13.2757, 15.0855,
     +  16.8112, 18.4748, 20.0900, 21.6661, 23.2097,
     +  24.7249, 26.2169, 27.6882, 29.1412, 30.5779,
     +  31.9999, 33.4086, 34.8053, 36.1908, 37.5662,
     +  38.9322, 40.2893, 41.6384, 42.9798, 44.3141,
     +  45.6417, 46.9629, 48.2782, 49.5879, 50.8922,
     +  52.1914, 53.4858, 54.7755, 56.0609, 57.3421,
     +  58.6192, 59.8925, 61.1621, 62.4281, 63.6907/

      DATA CHMIN/
     +   0.0002,  0.0201,  0.1149,  0.2972,  0.5543,
     +   0.8721,  1.2391,  1.6465,  2.0879,  2.5583,
     +   3.0536,  3.5707,  4.1070,  4.6605,  5.2294,
     +   5.8122,  6.4078,  7.0149,  7.6328,  8.2604,
     +   8.8972,  9.5425, 10.1957, 10.8564, 11.5240,
     +  12.1982, 12.8785, 13.5647, 14.2565, 14.9535,
     +  15.6555, 16.3622, 17.0735, 17.7892, 18.5089,
     +  19.2327, 19.9602, 20.6914, 21.4262, 22.1643/
*
      IERR=0
*
*   decode and check data descriptor and measurement code
*
      IDET= JBYT(DATE(1),1,5)
      IDD=  JBYT(DATE(4),1,5)
      IMEAS=JBYT(DATE(4),6,6)
      ICOORD=MOD(IDD,2)
      IMOM=IMEAS/32
      MEASUR=IMEAS
      IF(IMOM.EQ.1) MEASUR=MEASUR-16
      NCOV=JBYT(DATE(5),1,16)-17
*  If dE/dx measured, removes its error
      IF(IDET.EQ.detid(i_TPC).AND.NCOV.EQ.16) NCOV=NCOV-1
      NQXY=JBIT(MEASUR,1)+JBIT(MEASUR,4)+JBIT(MEASUR,5)
      NQZ=JBIT(MEASUR,2)+JBIT(MEASUR,3)
      NQ=NQXY+NQZ
*
*   check consistency of data
*
      IF((NQ*(NQ+1))/2.NE.NCOV) THEN
        print*, idet,nq,ncov,imeas
        IERR=1
        RETURN
      ENDIF
*
*   store chi2 and number of degrees of freedom
*
       ND=JBYT(DATE(7),1,32)
       CH2=DATE(8)
*
*    restrict internal chi^2s to values corresponding between 1 and 99 %     mif
*
      IF (ND.GT.0.AND.ND.LE.40) THEN
        CH2=MIN(CH2,CHMAX(ND))
        CH2=MAX(CH2,CHMIN(ND))
      END IF
*
*   covariance matrix
*
      DO 1 I=1,NCOV
    1 COVI(I)=DATE(16+I)
*
*  input coordinates are r, r*Phi, z
*
      IF(ICOORD.EQ.1) THEN
        TYPE='CYLI'
        PARAM(1)=DATE(10)
        PARAM(2)=DATE(11)/DATE(10)
        PARAM(3)=DATE(12)
        PARAM(4)=DATE(13)
        PARAM(5)=DATE(14)-PARAM(2)
        PARAM(6)=CONSB*DATE(15)
*
*  input coordinates are x, y, z
*
      ELSE
        TYPE='PLAN'
        PARAM(1)=DATE(12)
        PARAM(2)=DATE(10)
        PARAM(3)=DATE(11)
        PARAM(4)=DATE(13)
        PARAM(5)=DATE(14)
        PARAM(6)=CONSB*DATE(15)
*   if plane at negative z, change geometrical sign of 1/R
        IF(PARAM(1).LT.0.) PARAM(6)=-PARAM(6)
      ENDIF
*
*   if last parameter is 1/p, divide 1/R by sin(theta)
      IF(MEASUR.GE.16.AND.IMOM.EQ.0) PARAM(6)=PARAM(6)/DSIN(PARAM(4))
*
*   invert covariance matrix and rearrange it to get the weight matrix
*   on 5 parameters in all cases
*
      CALL FKCTOW(PARAM,IMOM,TYPE,COVI,NQ,MEASUR,WEIGHT,IERR)
*
      IF(IERR.EQ.0) THEN
        IRNKXY=MIN0(3,IRNKXY+NQXY)
        IRNKZ=MIN0(2,IRNKZ+NQZ)
      else
         print*,'FKCTOW error',ierr
      ENDIF
*
      END
