*  
* $Id: fkxdcy.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:56  hvogt
* Brahms V308
*
*  
************************************************************************
*                                                                      *
      SUBROUTINE FKXDCY(PARAMI,PARAMF,DER)
*                                                                      *
*   compute derivatives for a propagation cylinder -> cylinder         *
*   knowing the initial and final parameters                           *
*                                                                      *
*   see routine FKXCYL for the definition of PARAMI,PARAMF,DER         *
*   and comments                                                       *
*                                                                      *
************************************************************************
*                                                                      *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION PARAMI(6),PARAMF(6),DER(8)
      PARAMETER (PI=3.141592654D00,TWOPI=PI*2.0D00,DPHIMN=1.0D-02)
*
*   radius of the helix, cot(theta) and sin(beta),cos(beta) at both ends
*
      RADI=PARAMI(1)
      RADF=PARAMF(1)
      RTRK=1./PARAMI(6)
      COTTH=1./DTAN(PARAMI(4))
      SINB=DSIN(PARAMI(5))
      COSB=DCOS(PARAMI(5))
      SINBF=DSIN(PARAMF(5))
      COSBF=DCOS(PARAMF(5))
*
*   variation DPHI of phi angle
*
      DPHI=PARAMF(2)+PARAMF(5)-PARAMI(2)-PARAMI(5)
CKMO      DPHI=DMOD(DPHI+2.*TWOPI+PI,TWOPI)-PI
      IF (DPHI.LT.-PI) DPHI = DPHI+TWOPI
      IF (DPHI.LT.-PI) DPHI = DPHI+TWOPI
      IF (DPHI.GT.PI) DPHI = DPHI-TWOPI
      IF (DPHI.GT.PI) DPHI = DPHI-TWOPI
*
*   if |DPHI| >  DPHIMN , use "exact" formulae
*
      IF(DABS(DPHI).GE.DPHIMN) THEN
        XC=RADI-RTRK*SINB
        YC=RTRK*COSB
        RC2=XC**2+YC**2
*
        CCPSI=RADI-RTRK*SINB
        SCPSI=RTRK*COSB
        CCPSF=RADF-RTRK*SINBF
        SCPSF=RTRK*COSBF
*
        CPSII=RTRK-RADI*SINB
        SPSII=-RADI*COSB
        CPSIF=RTRK-RADF*SINBF
        SPSIF=-RADF*COSBF
*
        SDPHI=DSIN(DPHI)
        CDPHI=DCOS(DPHI)
*
        DER(1)=0.
        FACT=-RTRK/SPSIF
        DER(2)=SDPHI*FACT
        DER(3)=FACT*RTRK*(1.-CDPHI)
        DER(4)=-RTRK*DPHI*(1.+COTTH**2)
        DER(5)=RTRK*COTTH*(RADF*CCPSF*SPSII/SPSIF-RADI*CCPSI)
     &         /RC2
        DER(6)=RTRK**2*COTTH*(-DPHI+SINBF/COSBF-
     &         (RADI*SCPSI+RADF*CCPSF*CPSII/SPSIF)/RC2)
        DER(7)=SPSII/SPSIF
        DER(8)=RTRK*(CPSIF-CPSII)/SPSIF
*
*   if |DPHI| < DPHIMN , use first order in 1/R
*
      ELSE
        RDPHI=RTRK*DPHI
        DER(1)=0.
        DER(2)=RDPHI/(RADF*COSBF)
        DER(3)=.5*RDPHI*DER(2)
        DER(4)=-RDPHI*(1.+COTTH**2)
        DER(5)=RADI*COTTH*DSIN(PARAMF(2)-PARAMI(2))/COSBF
        DER(6)=.5*RDPHI*DER(5)
        DER(7)=(RADI*COSB)/(RADF*COSBF)
        DER(8)=.5*RDPHI*(1.+DER(7))
      ENDIF
*
      END
