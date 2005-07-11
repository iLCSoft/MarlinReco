*  
* $Id: fkspar.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE FKSPAR(TYPE,RZREF,PAR,COV,DATSP)
*                                                                      *
*   transform internal parameters and covariance into TANAGRA format   *
*                                                                      *
************************************************************************
*
#include "include/fksave.inc"
#include "include/fkfild.inc"
*
      DOUBLE PRECISION PI,TWOPI,PAR,COV,CAPPHI,DAT,RZREF,DER1,DER2
      DIMENSION PAR(5),COV(15),DAT(21),DATSP(21)
      CHARACTER*4 TYPE
      PARAMETER (PI=3.14159265359D00,TWOPI=PI*2.0D00)
*
*   first case : reference surface is a cylinder ----------------
*
      IF(TYPE.EQ.'CYLI') THEN
        DAT(1)=RZREF
CKMO        CAPPHI=DMOD(PAR(1)+TWOPI,TWOPI)
        CAPPHI=PAR(1)
        IF (CAPPHI.LT.0.) CAPPHI = CAPPHI+TWOPI
        IF (CAPPHI.GT.TWOPI) CAPPHI = CAPPHI-TWOPI
        DAT(2)=DAT(1)*CAPPHI
        DAT(3)=PAR(2)
        DAT(4)=PAR(3)
CKMO        DAT(5)=DMOD(CAPPHI+PAR(4)+TWOPI,TWOPI)
        DAT(5)=CAPPHI+PAR(4)
        IF (DAT(5).LT.0.) DAT(5) = DAT(5)+TWOPI
        IF (DAT(5).GT.TWOPI) DAT(5) = DAT(5)-TWOPI
*
*  change from beta to phi = Phi + beta
*
        COV(10)=COV(10)+2.*COV(7)+COV(1)
        COV(7)=COV(7)+COV(1)
        COV(8)=COV(8)+COV(2)
        COV(9)=COV(9)+COV(4)
        COV(14)=COV(14)+COV(11)
*
*  change from Phi to r*Phi
*
        COV(1)=COV(1)*DAT(1)**2
        COV(2)=COV(2)*DAT(1)
        COV(4)=COV(4)*DAT(1)
        COV(7)=COV(7)*DAT(1)
        COV(11)=COV(11)*DAT(1)
*
*  second case : reference surface is a plane  ----------------
*
      ELSE
        DAT(1)=PAR(1)
        DAT(2)=PAR(2)
        DAT(3)=RZREF
        DAT(4)=PAR(3)
c --- -ME- protect phi
CKMO        DAT(5)=DMOD(PAR(4)+TWOPI,TWOPI)
        DAT(5)=PAR(4)
        IF (DAT(5).LT.0.) DAT(5) = DAT(5)+TWOPI
        IF (DAT(5).GT.TWOPI) DAT(5) = DAT(5)-TWOPI
      ENDIF
*
*   in both cases : change from 1/R to 1/p = sin(theta)/(R*CONSB) ------
*
      DAT(6)=PAR(5)*DSIN(DAT(4))/CONSB
*
      DER1=DAT(6)/DTAN(DAT(4))
      DER2=DAT(6)/PAR(5)
      COV(15)=COV(6)*DER1**2+2.*COV(13)*DER1*DER2+COV(15)*DER2**2
      COV(11)=COV(4)*DER1+COV(11)*DER2
      COV(12)=COV(5)*DER1+COV(12)*DER2
      COV(13)=COV(6)*DER1+COV(13)*DER2
      COV(14)=COV(9)*DER1+COV(14)*DER2
*
*   store covariance matrix in TANAGRA array
*
      CALL FKDPCO(COV,DAT(7),15)
*
      DO 1 I=1,21
C --- -KH- protect against underflows
      IF (DABS(DAT(I)).LT.1.D-30) DAT(I)=0.0D0
    1 DATSP(I)=DAT(I)
*
      END

