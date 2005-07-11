*  
* $Id: fkswch.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE FKSWCH(TYP,PAR,DQ,WG)
*                                                                      *
*   change from "plane" parameters to "barrel" parameters              *
*   at a given point (without propagation), or conversely,             *
*   and tranform the deviations from ref. traj. and the weight matrix  *
*                                                                      *
*   Input   :   TYP       :  type to be changed ('PLAN' or 'CYLI')     *
*               PAR(1-6)  :  plane  parameters                         *
*                             z(fixed),x,y,theta,phi,1/R               *
*                            or cylinder parameters                    *
*                             r(fixed),Phi,z,theta,beta=phi-Phi,1/R    *
*               DQ(1-5)   :  deviations from reference trajectory      *
*                             on x,y,theta,phi,1/R       (plane)       *
*                             or Phi,z,theta,beta,1/R    (cylinder)    *
*               WG(1-15)  :  weight matrix on these parameters         *
*                                                                      *
*   Output  :   PAR(1-6)  :  parameters with the new type              *
*               DQ(1-5)   :  deviations from reference trajectory      *
*               WG(1-15)  :  weight matrix on these parameters         *
*                                                                      *
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.14159265359D00,TWOPI=PI*2.0D00)
*
      CHARACTER*4 TYP
      DIMENSION PAR(6),DQ(5),WG(15)
*
*   transformation from 'PLAN' to 'CYLI' ------------------------------
*
      IF(TYP.EQ.'PLAN') THEN
*
*   transformation of parameters
*
        Z=PAR(1)
        X=PAR(2)
        Y=PAR(3)
        PHI=PAR(5)
*
        PAR(1)=DSQRT(X**2+Y**2)
CKMO        PAR(2)=DMOD(DATAN2(Y,X)+TWOPI,TWOPI)
        PAR(2)=DATAN2(Y,X)
        IF (PAR(2).LT.0.) PAR(2) = PAR(2)+TWOPI
        PAR(3)=Z
CKMO        PAR(5)=DMOD(PHI-PAR(2)+TWOPI+PI,TWOPI)-PI
        PAR(5)=PHI-PAR(2)
        IF (PAR(5).LT.-PI) PAR(5) = PAR(5)+TWOPI
        IF (PAR(5).GT.PI) PAR(5) = PAR(5)-TWOPI
*
*   transformation of deviations
*
        RCOSB=PAR(1)*DCOS(PAR(5))
        TANTH=DTAN(PAR(4))
        COSF=DCOS(PHI)
        SINF=DSIN(PHI)
        DQQ=(X*DQ(1)+Y*DQ(2))/RCOSB
        DQ(1)= (-SINF*DQ(1)+COSF*DQ(2))/RCOSB
        DQ(2)=-DQQ/TANTH
        DQ(4)=DQ(4)-PAR(6)*DQQ-DQ(1)
*
*   transformation of weight matrix  x,y,phi -> Phi,z,phi
*
        D11=-Y
        D21=X
        D12=-TANTH*COSF
        D22=-TANTH*SINF
        D42=-TANTH*PAR(6)
*
        W11   = D11*WG( 1) + D21*WG( 2)
        W12   = D11*WG( 2) + D21*WG( 3)
        W13   = D11*WG( 4) + D21*WG( 5)
        W14   = D11*WG( 7) + D21*WG( 8)
        W15   = D11*WG(11) + D21*WG(12)
        W21   = D12*WG( 1) + D22*WG( 2) + D42*WG( 7)
        W22   = D12*WG( 2) + D22*WG( 3) + D42*WG( 8)
        WG( 5)= D12*WG( 4) + D22*WG( 5) + D42*WG( 9)
        WG( 8)= D12*WG( 7) + D22*WG( 8) + D42*WG(10)
        WG(12)= D12*WG(11) + D22*WG(12) + D42*WG(14)
*
        WG( 1)= W11*D11 + W12*D21
        WG( 2)= W21*D11 + W22*D21
        WG( 3)= W21*D12 + W22*D22 + WG(8)*D42
*
*   transformation  Phi,phi -> Phi,beta
*
        WG( 1)= WG( 1) +2.*W14   + WG(10)
        WG( 2)= WG( 2) +   WG( 8)
        WG( 4)= W13    +   WG( 9)
        WG( 7)= W14    +   WG(10)
        WG(11)= W15    +   WG(14)
*
*   transformation from 'CYLI' to 'PLAN' ------------------------------
*
      ELSE
*
*   transformation of parameters
*
        RAD=PAR(1)
        CAPPHI=PAR(2)
        COSB=DCOS(PAR(5))
*
        PAR(1)=PAR(3)
        PAR(2)=RAD*DCOS(CAPPHI)
        PAR(3)=RAD*DSIN(CAPPHI)
CKMO        PAR(5)=DMOD(CAPPHI+PAR(5)+2.*TWOPI,TWOPI)
        PAR(5)=CAPPHI+PAR(5)
        IF (PAR(5).LT.0.) PAR(5) = PAR(5)+TWOPI
        IF (PAR(5).LT.0.) PAR(5) = PAR(5)+TWOPI
        IF (PAR(5).GT.TWOPI) PAR(5) = PAR(5)-TWOPI
        COSF=DCOS(PAR(5))
        SINF=DSIN(PAR(5))
        TANTH=DTAN(PAR(4))
*
*   transformation of deviations
*
        DQ(4)=DQ(4)-TANTH*PAR(6)*DQ(2)+DQ(1)
        DQ1= -PAR(3)*DQ(1)-TANTH*COSF*DQ(2)
        DQ(2)=PAR(2)*DQ(1)-TANTH*SINF*DQ(2)
        DQ(1)=DQ1
*
*   transformation of weight matrix  Phi,beta -> Phi,phi
*
        WG( 1)= WG( 1) -2.*WG( 7) + WG(10)
        WG( 2)= WG( 2) -   WG( 8)
        WG( 4)= WG( 4) -   WG( 9)
        WG( 7)= WG( 7) -   WG(10)
        WG(11)= WG(11) -   WG(14)
*
*   transformation of weight matrix  Phi,z,phi -> x,y,phi
*
        RCOSB=RAD*COSB
        D11=-SINF/RCOSB
        D12= COSF/RCOSB
        D21=-PAR(2)/(TANTH*RCOSB)
        D22=-PAR(3)/(TANTH*RCOSB)
        D41=-PAR(2)*PAR(6)/RCOSB
        D42=-PAR(3)*PAR(6)/RCOSB
*
        W11   = D11*WG( 1) + D21*WG( 2) + D41*WG( 7)
        W12   = D11*WG( 2) + D21*WG( 3) + D41*WG( 8)
        W13   = D11*WG( 4) + D21*WG( 5) + D41*WG( 9)
        W14   = D11*WG( 7) + D21*WG( 8) + D41*WG(10)
        W15   = D11*WG(11) + D21*WG(12) + D41*WG(14)
        W21   = D12*WG( 1) + D22*WG( 2) + D42*WG( 7)
        W22   = D12*WG( 2) + D22*WG( 3) + D42*WG( 8)
        WG( 5)= D12*WG( 4) + D22*WG( 5) + D42*WG( 9)
        WG( 8)= D12*WG( 7) + D22*WG( 8) + D42*WG(10)
        WG(12)= D12*WG(11) + D22*WG(12) + D42*WG(14)
*
        WG( 1)= W11*D11 + W12*D21 + W14  *D41
        WG( 2)= W21*D11 + W22*D21 + WG(8)*D41
        WG( 3)= W21*D12 + W22*D22 + WG(8)*D42
        WG( 4)= W13
        WG( 7)= W14
        WG(11)= W15
*
      ENDIF
*
      END
