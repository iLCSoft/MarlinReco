*  
* $Id: fkxdpl.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE FKXDPL(PARAMI,PARAMF,DER)
*                                                                      *
*   compute derivatives for a propagation plane -> plane               *
*   knowing the initial and final parameters                           *
*                                                                      *
*   see routine FKXPLA for the definition of PARAMI,PARAMF,DER         *
*   and comments                                                       *
*                                                                      *
************************************************************************
*                                                                      *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PARAMI(6),PARAMF(6),DER(8)
*
      PARAMETER (PI=3.141592654D00,TWOPI=PI*2.0D00,DPHIMN=1.0D-02)
*
      COSF0=DCOS(PARAMI(5))
      COSF1=DCOS(PARAMF(5))
      SINF0=DSIN(PARAMI(5))
      SINF1=DSIN(PARAMF(5))
      TANTH=DTAN(PARAMI(4))
      DZ=PARAMF(1)-PARAMI(1)
      RDPHI=TANTH*DZ
      DER(1)=(1.0+TANTH**2)*DZ*COSF1
      DER(4)=(1.0+TANTH**2)*DZ*SINF1
      DER(7)=DZ*PARAMI(6)*(1.0+TANTH**2)
      DER(8)=RDPHI
*
      DPHI=RDPHI*PARAMI(6)
*
*   "exact" formulae if |DPHI| > DPHIMN
*
      IF(DABS(DPHI).GE.DPHIMN) THEN
        RTRK=1.0/PARAMI(6)
        DCOSF=COSF1-COSF0
        DSINF=SINF1-SINF0
        DER(2)=RTRK*DCOSF
        DER(3)=RTRK**2*(DPHI*COSF1-DSINF)
        DER(5)=RTRK*DSINF
        DER(6)=RTRK**2*(DPHI*SINF1+DCOSF)
*
*   first order in 1/R if |DPHI| < DPHIMN
      ELSE
        DER(2)=-RDPHI*SINF0
        DER(3)=0.5*RDPHI*DER(2)
        DER(5)=RDPHI*COSF0
        DER(6)=0.5*RDPHI*DER(5)
      ENDIF
*
      END
