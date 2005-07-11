*  
* $Id: fkctow.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
      SUBROUTINE FKCTOW(PAR,IMOM,TYPE,COV,NQ,MEASUR,WGT,IERR)
*                                                                      *
*   transformation of covariance matrix in TANAGRA format into a       *
*   fixed format (5x5) weight matrix. Missing elements are set to zero *
*                                                                      *
*   Input  :  PAR  :  parameters                                       *
*             IMOM :  1 if word 15 of TER is 1/pt ; 0 if 1/p           *
*             TYPE :  type ('CYLI' or 'PLAN')                          *
*             COV  :  words 17 and following of TER                    *
*             NQ   :  number of quantities measured in this TE         *
*                                                                      *
*   Output :  WGT  :  weight matrix                                    *
*             IERR :  error flag : 0 if OK                             *
*                                  2 if error matrix not positive      *
*                                                                      *
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "include/fksave.inc"
#include "include/fkfild.inc"
*
      DIMENSION PAR(6),COV(15),WGP(15),WGT(15),KM(5)
*
      CHARACTER*4 TYPE
*
*   invert packed covariance matrix and check positivity
*   NQ = number of measured quantities (dimension of TANAGRA covariance matrix)
*
      IF(NQ.EQ.1) THEN
        IF(COV(1).LE.0.) THEN
          IERR=2
        ELSE
          WGP(1)=1./COV(1)
        ENDIF
      ELSE IF(NQ.EQ.2) THEN
        DET=COV(1)*COV(3)-COV(2)**2
        IF(DET.LE.0.) THEN
          IERR=2
        ELSE
          WGP(1)=COV(3)/DET
          WGP(2)=-COV(2)/DET
          WGP(3)=COV(1)/DET
        ENDIF
      ELSE IF(NQ.EQ.3) THEN
        CALL FKMI3D(COV,WGP,IERR)
      ELSE IF(NQ.EQ.4) THEN
        CALL FKMI4D(COV,WGP,IERR)
      ELSE
        CALL FKMI5D(COV,WGP,IERR)
      ENDIF
      IF(IERR.NE.0) THEN
         PRINT*,NQ,COV(1),COV(2),COV(3)
         RETURN
      ENDIF
*
*   assign to each measured quantity a standard numbering from 1 to 5
*
      K=0
      DO 1 I=1,5
        IF(JBIT(MEASUR,I).EQ.1) THEN
          K=K+1
          KM(K)=I
        ENDIF
    1 CONTINUE
*
*   expand weight matrix from (NQxNQ) to (5x5) :
*   fill WGT with zeroes and then copy each term of WGP to its standard
*   location in WGT
*
      CALL FKDPV0(WGT,15)
      L=0
      DO 2 I=1,NQ
        II=(KM(I)*(KM(I)-1))/2
        DO 2 J=1,I
          L=L+1
          WGT(II+KM(J))=WGP(L)
    2 CONTINUE
*
*   if barrel parameters : change parameters  r*Phi,phi  to  Phi,beta
*
      IF(TYPE.EQ.'CYLI') THEN
*   change r*Phi to Phi
        WGT( 1)=WGT( 1)*PAR(1)**2
        WGT( 2)=WGT( 2)*PAR(1)
        WGT( 4)=WGT( 4)*PAR(1)
        WGT( 7)=WGT( 7)*PAR(1)
        WGT(11)=WGT(11)*PAR(1)
*
*   change phi to beta=phi-Phi
        WGT( 1)=WGT( 1)+2*WGT( 7)+WGT(10)
        WGT( 2)=WGT( 2)+WGT( 8)
        WGT( 4)=WGT( 4)+WGT( 9)
        WGT( 7)=WGT( 7)+WGT(10)
        WGT(11)=WGT(11)+WGT(14)
*
*   if plane parameters and z negative, change geometrical sign of 1/p or 1/pt
*
      ELSE IF(MEASUR.GE.16.AND.PAR(1).LT.0.) THEN
        WGT(11)=-WGT(11)
        WGT(12)=-WGT(12)
        WGT(13)=-WGT(13)
        WGT(14)=-WGT(14)
*
      ENDIF
*
*   change 1/p or 1/pt to 1/R
*
      IF(MEASUR.GE.16) THEN
*
*   input is 1/pt = 1/CONSB * 1/R
*
        IF(IMOM.EQ.1) THEN
          WGT(11)=WGT(11)/CONSB
          WGT(12)=WGT(12)/CONSB
          WGT(13)=WGT(13)/CONSB
          WGT(14)=WGT(14)/CONSB
          WGT(15)=WGT(15)/CONSB**2
*
*   input is 1/p = 1/CONSB * sin(theta) * 1/R
*
        ELSE
          D1=DCOS(PAR(4))*PAR(6)/CONSB
          D2=DSIN(PAR(4))/CONSB
*
          WGT( 6)= WGT( 6)+ 2.*D1*WGT(13)+ D1**2 *WGT(15)
          WGT(13)= D2*WGT(13)+ D1*D2*WGT(15)
          WGT(15)= D2**2 *WGT(15)
*
          WGT( 4)= WGT( 4)+ D1*WGT(11)
          WGT( 5)= WGT( 5)+ D1*WGT(12)
          WGT( 9)= WGT( 9)+ D1*WGT(14)
          WGT(11)= D2*WGT(11)
          WGT(12)= D2*WGT(12)
          WGT(14)= D2*WGT(14)
        ENDIF
      ENDIF
*
      END
