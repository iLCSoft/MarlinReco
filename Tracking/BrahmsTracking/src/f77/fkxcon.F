*  
* $Id: fkxcon.F,v 1.1 2006-11-02 12:03:41 rasp Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
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
      SUBROUTINE FKXCON(PARAMI,IDIR,ZMIN,ZMAX,RADMIN,RADMAX,
     &                  SINBMX,IOPT,
     &                  PARAMF,DER,ALRPHI,IERR)
*                                                                      *
*    AIM :                                                             *
*    -----                                                             *
*    Extrapolate a helix defined by the initial parameters PARAMI      *
*    up to a given cone, and compute if requested the derivatives      *
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
*              ZMIN          : lower z limit of the cylinder           *
*              ZMAX          : upper z limit of the cylinder           *
*              RADMIN        : radius defined at ZMIN                  * 
*              RADMAX        : radius defined at ZMAX                  *
*              SINBMX        : maximum allowed for |sin(beta)| at the  *
*                              intersection                            *
*              IOPT          : 0 if derivatives not requested          *
*                              1 if derivatives requested              *
*                                                                      *
*    Output :  IERR          : 0 if intersection found                 *
*                              1 if no intersection with the cone      *
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
*    Author  :  A. Raspereza                                           *
*                                                                      *
*    First version : 26-01-88                                          *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL RADMIN,RADMAX,ZMIN,ZMAX,SINBMX,ALRPHI
      REAL ZMINCYL,ZMAXCYL,RADCYL,TEMP
*
      DIMENSION PARAMI(6),PARAMF(6),DER(8)
*
      PARAMETER (PI=3.141592654D00,TWOPI=PI*2.0D00,DPHIMN=1.0D-04)
*
      COMMON /CONICAL_EQUATION/ APHI,BPHI,CPHI,DPHI,PHIDOT

      DOUBLE PRECISION FCONEQ
      EXTERNAL FCONEQ

      IERR=0
*
c$$$      WRITE(*,*) 'PARAMETERS = ',PARAMI(1),' ',PARAMI(2),' ',PARAMI(3),
c$$$     &                           PARAMI(4),' ',PARAMI(5),' ',PARAMI(6)

      RADI=PARAMI(1)
      PHII=PARAMI(2)
      BETA=PARAMI(5)
      COSB=DCOS(BETA)
      SINB=DSIN(BETA)
      TANT=DTAN(PARAMI(4))
      COTTH=1./DTAN(PARAMI(4))
      RTRK=1./PARAMI(6)
*
*  XREF and YREF are x and y coordinates at reference point  
*  PHIREF - phi angle of track momentum at reference point
*  PHII = ATAN2(YREF,XREF)
*
      PHIREF = BETA + PHII
      XREF = RADI * DCOS(PHII)
      YREF = RADI * DSIN(PHII)
      ZREF = PARAMI(3)

*
*  helix parameterization used :
*  X = X0 + R*COS(P)
*  Y = Y0 + R*SIN(P)
*  Z = Z0 + BZ*P
*
*  find helix parameters
*

      BZ = COTTH * RTRK
      X0 = XREF - RTRK * DSIN(PHIREF)
      Y0 = YREF + RTRK * DCOS(PHIREF)
      PHI_REF = DATAN2(YREF-Y0,XREF-X0)
      Z0 = ZREF - BZ*PHI_REF

*
* transcendental equation to be solved
* APHI + BPHI*P + CPHI*P**2 + DPHI*COS(P-PHIDOT) = 0 
* 
* find coefficients 
*

      A = RADMIN - ZMIN*(RADMAX-RADMIN)/(ZMAX-ZMIN)
      B = (RADMAX-RADMIN)/(ZMAX-ZMIN)
      R02 = X0*X0+Y0*Y0
      R0 = DSQRT(R02)
      PHIDOT = DATAN2(Y0,X0)

      TWOAB = 2*A*B
      BSQUARE = B*B

      APHI = A*A + TWOAB*Z0 + BSQUARE*Z0*Z0 - R02 - RTRK*RTRK
      BPHI = TWOAB * BZ + 2*BSQUARE*Z0*BZ
      CPHI = BSQUARE * BZ*BZ
      DPHI = -2.0D0*DABS(RTRK)*R0
      BZINV = 1.0D0/BZ

      PMIN = (ZMIN - Z0)*BZINV
      PMAX = (ZMAX - Z0)*BZINV

      PRODUCT = FCONEQ(PMIN)*FCONEQ(PMAX)
      
      IF (PRODUCT>0.0) THEN
C         WRITE(*,*) 'INTERSECTION NOT FOUND '
         IERR = 1
         RETURN
      ENDIF

      MAXF = 1000
      EPS  = 1.0E-8
      MODE = 1

C      WRITE(*,*) 'X0,Y0,Z0 = ',X0,' ',Y0,' ',Z0
C      WRITE(*,*) 'BZ = ',BZ
C      WRITE(*,*) 'RTRK = ',RTRK
C      WRITE(*,*) 'APHI,BPHI,CPHI,DPHI,PHIDOT = ',
C     &            APHI,' ',BPHI,' ',CPHI,' ',DPHI,' ',PHIDOT

      P0 = DZEROX(PMIN,PMAX,EPS,MAXF,FCONEQ,MODE)
*
* SOME DIAGNOSTICS
*

*      WRITE(*,*) 

      XHELIX = X0 + DABS(RTRK)*DCOS(P0)
      YHELIX = Y0 + DABS(RTRK)*DSIN(P0)
      RHELIX = DSQRT(XHELIX*XHELIX+YHELIX*YHELIX)
      ZHELIX = Z0 + BZ*P0
      ZCONE = ZHELIX
      RCONE = A + B*ZCONE
C      WRITE(*,*) 'Z, R, X, Y ON HELIX = ',ZHELIX,RHELIX,
C     &                                    XHELIX,YHELIX
C      WRITE(*,*) 'Z AND R ON CONE     = ',ZCONE,RCONE

      ZMINCYL = SNGL(ZCONE) - 1.0E-3
      ZMAXCYL = SNGL(ZCONE) + 1.0E-3
      RADCYL  = SNGL(RCONE)

      CALL FKXCYL(PARAMI,IDIR,RADCYL,ZMINCYL,ZMAXCYL,SINBMX,IOPT,
     &            PARAMF,DER,ALPHI,IERR)


C      WRITE(*,*) 'AFTER CALLING FKXCYL -----> IERR = ',IERR
      XFIN = PARAMF(1)*DCOS(PARAMF(2))
      YFIN = PARAMF(1)*DSIN(PARAMF(2))
      ZFIN = PARAMF(3)
      RFIN = DSQRT(XFIN*XFIN+YFIN*YFIN)
      
C      WRITE(*,*) 'Z, R, X, Y FINAL    = ',ZFIN,RFIN,XFIN,YFIN


*
*  
*       

*
      END



      DOUBLE PRECISION FUNCTION FCONEQ(X)
*
*  NON-LINEAR FUNCTION TO BE SLOVED IN ORDER TO FIND
*  HELIX INTERSECTION WITH CONE
*
      DOUBLE PRECISION X
      DOUBLE PRECISION APHI,BPHI,CPHI,DPHI,PHIDOT
      COMMON /CONICAL_EQUATION/ APHI,BPHI,CPHI,DPHI,PHIDOT

      FCONEQ = APHI + BPHI*X + CPHI*X*X + DPHI*DCOS(X-PHIDOT)


      RETURN
      END
