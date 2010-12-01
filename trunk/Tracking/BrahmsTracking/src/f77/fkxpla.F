*  
* $Id: fkxpla.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE FKXPLA (PARAMI,IDIR,ZF,RMIN,RMAX,IOPT,
     &                   PARAMF,DER,ALRPHI,IERR)
*                                                                      *
************************************************************************
*                                                                      *
*    AIM :                                                             *
*    -----                                                             *
*    extrapolate a helix defined by the initial parameters PARAMI      *
*    up to a given z-plane, and compute the derivatives of the         *
*    final parameters w.r.t. the initial ones                          *
*                                                                      *
*    The computation uses double precision on intermediate variables   *
*    If the variation of phi angle is less than DPHIMN (.0001 in this  *
*    version) the computation is done at first order in 1/R in order   *
*    to avoid rounding errors in the derivatives                       *
*                                                                      *
*    Input  :  PARAMI(1-6)   : initial parameters                      *
*                              (z,x,y,theta,phi,1/R)                   *
*              ZF            : z of the final plane                    *
*              IDIR    :  if  1 : positive extrapolation only          *
*                         if -1 : negative         "                   *
*                         if  0 : extrapolation on both sides          *
*              RMIN          : lower limit of r on the plane           *
*              RMAX          : upper limit of r on the plane           *
*              IOPT          : 0 if derivatives not requested          *
*                              1 if derivatives requested              *
*                                                                      *
*    Output :  IERR          : 0 if OK                                 *
*                              1 if no intersection found              *
*                              3 if intersection outside of limits     *
*              PARAMF(1-6)   : final parameters                        *
*              DER(1-8)      : deriv. of final w.r.t. initial param.   *
*                              DER(1) = d(x)/d(theta)                  *
*                              DER(2) = d(x)/d(phi)                    *
*                              DER(3) = d(x)/d(1/R)                    *
*                              DER(4) = d(y)/d(theta)                  *
*                              DER(5) = d(y)/d(phi)                    *
*                              DER(6) = d(y)/d(1/R)                    *
*                              DER(7) = d(phi)/d(theta)                *
*                              DER(8) = d(phi)/d(1/R)                  *
*              ALRPHI        : length (in r-phi projection) from start *
*                              to extrapolation, with a sign (positive *
*                              if the extrapolation is towards the     *
*                              direction defined by theta,phi)         *
*                                                                      *
*    Author  :  P. Billoir                                             *
*                                                                      *
*    First version : 26-01-88                                          *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL ZF,RMIN,RMAX,ALRPHI
*
      DIMENSION PARAMI(6),PARAMF(6),DER(8)
*
      PARAMETER (PI=3.141592654D00,TWOPI=PI*2.0D00,DPHIMN=1.0D-04)
*
      IERR=0
*
      DZ=ZF-PARAMI(1)
      TANTH=DTAN(PARAMI(4))
      RDPHI=DZ*TANTH
*
*   select positive or negative extrapolation, or both
      IF(RDPHI*IDIR.LT.0.) THEN
        IERR=1
      ELSE
        ALRPHI=RDPHI
*
*   x,y coordinates of intersection
        RTRK=1./PARAMI(6)
        COSF0=DCOS(PARAMI(5))
        SINF0=DSIN(PARAMI(5))
        XC=PARAMI(2)-RTRK*SINF0
        YC=PARAMI(3)+RTRK*COSF0
        DPHI=PARAMI(6)*RDPHI
CKMO        PHI1=DMOD(PARAMI(5)+DPHI,TWOPI)
        PHI1=PARAMI(5)+DPHI
        COSF1=DCOS(PHI1)
        SINF1=DSIN(PHI1)
        X1=XC+RTRK*SINF1
        Y1=YC-RTRK*COSF1
        R1=DSQRT(X1**2+Y1**2)
*
*   intersection outside of limits in r
        IF(R1.LT.RMIN.OR.R1.GT.RMAX) THEN
          IERR=3
*
*   parameters at the intersection
        ELSE
          PARAMF(1)=ZF
          PARAMF(2)=X1
          PARAMF(3)=Y1
          PARAMF(4)=PARAMI(4)
          PARAMF(5)=PHI1
          IF(PARAMF(5).LT.0.) PARAMF(5)=PARAMF(5)+TWOPI
          IF(PARAMF(5).GT.TWOPI) PARAMF(5)=PARAMF(5)-TWOPI
          PARAMF(6)=PARAMI(6)
*
*   computation of derivatives -----------------------------------
*
          IF(IOPT.EQ.1) THEN
            CT2INV=1.+TANTH**2
            DER(1)=CT2INV*DZ*COSF1
            DER(4)=CT2INV*DZ*SINF1
            DER(7)=DZ*PARAMI(6)*CT2INV
            DER(8)=RDPHI
*
*   "exact" formulae if |DPHI| > DPHIMN
            IF(DABS(DPHI).GE.DPHIMN) THEN
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
              DER(3)=.5*RDPHI*DER(2)
              DER(5)=RDPHI*COSF0
              DER(6)=.5*RDPHI*DER(5)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
*
      END
