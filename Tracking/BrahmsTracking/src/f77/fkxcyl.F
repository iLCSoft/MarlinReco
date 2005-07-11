*  
* $Id: fkxcyl.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE FKXCYL(PARAMI,IDIR,RADF,ZMIN,ZMAX,SINBMX,IOPT,
     &                    PARAMF,DER,ALRPHI,IERR)
*                                                                      *
*    AIM :                                                             *
*    -----                                                             *
*    Extrapolate a helix defined by the initial parameters PARAMI      *
*    up to a given cylinder, and compute if requested the derivatives  *
*    of the final parameters w.r.t. the initial ones                   *
*                                                                      *
*    The computation uses double precision on intermediate variables   *
*    If the variation of phi angle is less than DPHIMN (.0001 in this  *
*    version) the computation is done at first order in 1/R in order   *
*    to avoid rounding errors, especially in the derivatives           *
*                                                                      *
*    Input  :  PARAMI(1-6)   : initial parameters                      *
*                              (r,Phi,z,theta,beta,1/R)                *
*                              with beta = phi-Phi                     *
*                                       geometrical sign)              *
*              IDIR    :  if  1 : positive extrapolation only          *
*                         if -1 : negative         "                   *
*                         if  0 : extrapolation on both sides          *
*              RADF          : radius of the cylinder                  *
*              ZMIN          : lower z limit of the cylinder           *
*              ZMAX          : upper z limit of the cylinder           *
*              SINBMX        : maximum allowed for |sin(beta)| at the  *
*                              intersection                            *
*              IOPT          : 0 if derivatives not requested          *
*                              1 if derivatives requested              *
*                                                                      *
*    Output :  IERR          : 0 if intersection found                 *
*                              1 if no intersection with the cylinder  *
*                              2 if SINBMX exceeded                    *
*                              3 if intersection outside of limits     *
*              PARAMF(1-6)   : final parameters                        *
*              DER(1-8)      : deriv. of final w.r.t. initial param.   *
*                              DER(1) = d(Phi)/d(theta)                *
*                              DER(2) = d(Phi)/d(beta)                 *
*                              DER(3) = d(Phi)/d(1/R)                  *
*                              DER(4) = d(z)/d(theta)                  *
*                              DER(5) = d(z)/d(beta)                   *
*                              DER(6) = d(z)/d(1/R)                    *
*                              DER(7) = d(beta)/d(beta)                *
*                              DER(8) = d(beta)/d(1/R)                 *
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
      REAL RADF,ZMIN,ZMAX,SINBMX,ALRPHI
*
      DIMENSION PARAMI(6),PARAMF(6),DER(8)
*
      PARAMETER (PI=3.141592654D00,TWOPI=PI*2.0D00,DPHIMN=1.0D-04)
*
      IERR=0
*
      RADI=PARAMI(1)
      PHII=PARAMI(2)
      BETA=PARAMI(5)
      COSB=DCOS(BETA)
      SINB=DSIN(BETA)
      COTTH=1./DTAN(PARAMI(4))
      RTRK=1./PARAMI(6)
*
*  center and squared radius of the projected helix, in a rotated frame
*  (x-axis through the starting point)
      XC=RADI-RTRK*SINB
      YC=RTRK*COSB
      RC2=XC**2+YC**2
*
*  resolution of a second order equation
      RRR=(RADF**2-RTRK**2-RC2)/(2.*RTRK)
      DELT=RC2-RRR**2
*
*   intersection exists if DELT > 0
*
      IF(DELT.LE.0.) THEN
        IERR=1
      ELSE
        DELT=DSQRT(DELT)
*
*   choose intersection on the same side as the starting point
*   w.r.t. the plane containing the z axis and the axis of the helix
*
        SINF=(XC*RRR+YC*DELT)/RC2
        COSF=(XC*DELT-YC*RRR)/RC2
        XF=XC+RTRK*SINF
        YF=YC-RTRK*COSF
        SINBF=(SINF*XF-COSF*YF)/RADF
*
*   exit if beta too large at the intersection
        IF(DABS(SINBF).GT.SINBMX) THEN
          IERR=2
        ELSE
          ALFF=DATAN2(SINF,COSF)
          DPHI=ALFF-BETA
          ALRPHI=RTRK*DPHI
*
*   select positive or negative extrapolations, or both
*
          IF(ALRPHI*IDIR.LT.0.) THEN
            IERR=1
            RETURN
          ENDIF
*
*   switch to approximate expressions if the variation of phi
*   is less than DPHIMN
*
*   "exact" expressions ---------------------------------------
*
          IF(DABS(DPHI).GE.DPHIMN) THEN
*
            ZF=PARAMI(3)+COTTH*RTRK*DPHI
*
*   exit if outside of limits in z,theta,Phi
            TTH=PI/2.-DATAN(ZF/RADF)
            PHIF=DATAN2(YF,XF)
            PPH=PHIF
            IF(PPH.LT.0.) PPH=PPH+TWOPI
            IF(ZF.LT.ZMIN.OR.ZF.GT.ZMAX) THEN
              IERR=3
            ELSE
*
*   final parameters
              PARAMF(1)=RADF
CKMO              PARAMF(2)=DMOD(PHII+PHIF+TWOPI,TWOPI)
              PARAMF(2)=PHII+PHIF
              IF (PARAMF(2).LT.0.) PARAMF(2) = PARAMF(2)+TWOPI
              IF (PARAMF(2).GT.TWOPI) PARAMF(2) = PARAMF(2)-TWOPI
              PARAMF(3)=ZF
              PARAMF(4)=PARAMI(4)
              PARAMF(5)=ALFF-PHIF
              PARAMF(6)=PARAMI(6)
*
*   computation of derivatives ---------
*
              IF(IOPT.EQ.1) THEN
                COSBF=DSQRT(1.-SINBF**2)
*
*  CCPSI = RC*cos(cap.psi) ; SCPSI = RC*sin(cap.psi)     (initial point)
*  CCPSF = RC*cos(cap.psi) ; SCPSF = RC*sin(cap.psi)     (final point)
                CCPSI=RADI-RTRK*SINB
                SCPSI=RTRK*COSB
                CCPSF=RADF-RTRK*SINBF
                SCPSF=RTRK*COSBF
*
*  CPSII = sgn*RC*cos(psi) ; SPSII = sgn*RC*sin(psi)     (initial point)
*  CPSIF = sgn*RC*cos(psi) ; SPSIF = sgn*RC*sin(psi)     (final point)
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
     &                 /RC2
                DER(6)=RTRK**2*COTTH*(-DPHI+SINBF/COSBF-
     &                 (RADI*SCPSI+RADF*CCPSF*CPSII/SPSIF)/RC2)
                DER(7)=SPSII/SPSIF
                DER(8)=RTRK*(CPSIF-CPSII)/SPSIF
              ENDIF
            ENDIF
*
*   approximation at first order in 1/R --------------------------
*
          ELSE
            DR2=RADF**2-RADI**2
            RCOSB=RADI*COSB
            AA=1.-RADI*SINB/RTRK
            DELT=RCOSB**2+AA*DR2
*   exit if no solution
            IF(DELT.LE.0.) THEN
             IERR=1
            ELSE
              RDPHI=(DSQRT(DELT)-RCOSB)/AA
              DPHI=RDPHI/RTRK
              DCOSF=-SINB-.5*COSB*DPHI
              COSF=COSB+DCOSF*DPHI
              YF=-RDPHI*DCOSF
              DSINF=COSB-.5*SINB*DPHI
              SINF=SINB+DSINF*DPHI
              XF=RADI+RDPHI*DSINF
              SINBF=(SINF*XF-COSF*YF)/RADF
              ZF=PARAMI(3)+COTTH*RDPHI
*
*
*   exit if outside of limits in z,theta,Phi
              TTH=PI/2.-DATAN(ZF/RADF)
              PHIF=DATAN2(YF,XF)
              PPH=PHIF
              IF(PPH.LT.0.) PPH=PPH+TWOPI
              IF(ZF.LT.ZMIN.OR.ZF.GT.ZMAX) THEN
                IERR=3
              ELSE
                PARAMF(1)=RADF
                PHIF=DATAN2(YF,XF)
CKMO                PARAMF(2)=DMOD(PARAMI(2)+PHIF+TWOPI,TWOPI)
                PARAMF(2)=PARAMI(2)+PHIF
                IF (PARAMF(2).LT.0.) PARAMF(2) = PARAMF(2)+TWOPI
                IF (PARAMF(2).GT.TWOPI) PARAMF(2) = PARAMF(2)-TWOPI
                PARAMF(3)=ZF
                PARAMF(4)=PARAMI(4)
                PARAMF(5)=PARAMI(5)+DPHI-PHIF
                PARAMF(6)=PARAMI(6)
*
*   Computation of derivatives --------------
*
                IF(IOPT.EQ.1) THEN
                  COSBF=DSQRT(1.-SINBF**2)
                  SPHIF=YF/RADF
*
                  DER(1)=0.
                  DER(2)=RDPHI/(RADF*COSBF)
                  DER(3)=.5*RDPHI*DER(2)
                  DER(4)=-RDPHI*(1.+COTTH**2)
                  DER(5)=RADI*COTTH*SPHIF/COSBF
                  DER(6)=.5*RDPHI*DER(5)
                  DER(7)=(RADI*COSB)/(RADF*COSBF)
                  DER(8)=.5*RDPHI*(1.+DER(7))
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
*
      END
