*  
* $Id: fkxpln.F,v 1.1 2007-11-23 08:50:09 rasp Exp $
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
*************************************************************************
*                                                                       *
      SUBROUTINE FKXPLN(PARAMI,IDIR,XPOS,YPOS,ZPOS,WIDTH,LENGTH,
     &                  LPHI,SINBMX,IOPT,PARAMF,DER,ALRPHI,COSNOR,IERR) *
*    AIM :                                                              *
*    -----                                                              *
*    Extrapolate a helix defined by the initial parameters PARAMI       *
*    up to a given ladder, and compute if requested the derivatives     *
*    of the final parameters w.r.t. the initial ones                    *
*                                                                       *
*    The computation uses double precision on intermediate variables    *
*    If the variation of phi angle is less than DPHIMN (.0001 in this   *
*    version) the computation is done at first order in 1/R in order    *
*    to avoid rounding errors, especially in the derivatives            *
*                                                                       *
*    Input  :  PARAMI(1-6)   : initial parameters                       *
*                              (r,Phi,z,theta,beta,1/R)                 *
*                              with beta = phi-Phi                      *
*                                       geometrical sign)               *
*              IDIR          :  if  1 : positive extrapolation only     *
*                               if -1 : negative         "              *
*                               if  0 : extrapolation on both sides     *
*              XPOS          : x coordinate of the center of the ladder *
*              YPOS          : y               "                        *
*              ZPOS          : z               "                        *
*              WIDTH         : width of the ladder                      *
*              LENGTH        : length of the ladder (parallel to z axis)*
*              LPHI          : angle between normal of ladder and x axis*
*              SINBMX        : maximum allowed for |sin(beta)| at the   *
*                              intersection                             *
*              IOPT          : 0 if derivatives not requested           *
*                              1 if derivatives requested               *
*                                                                       *
*    Output :  IERR          : 0 if intersection found                  *
*                              1 if no intersection with the ladder     *
*                              2 if SINBMX exceeded                     *
*                              3 if intersection outside of limits      *
*              PARAMF(1-6)   : final parameters                         *
*              DER(1-8)      : deriv. of final w.r.t. initial param.    *
*                              DER(1) = d(Phi)/d(theta)                 *
*                              DER(2) = d(Phi)/d(beta)                  *
*                              DER(3) = d(Phi)/d(1/R)                   *
*                              DER(4) = d(z)/d(theta)                   *
*                              DER(5) = d(z)/d(beta)                    *
*                              DER(6) = d(z)/d(1/R)                     *
*                              DER(7) = d(beta)/d(beta)                 *
*                              DER(8) = d(beta)/d(1/R)                  *
*              ALRPHI        : length (in r-phi projection) from start  *
*                              to extrapolation, with a sign (positive  *
*                              if the extrapolation is towards the      *
*                              direction defined by theta,phi)          *
*              COSNOR        : cosine of the angle between the normal   *
*                              of the ladder and the tangent            *
*                              of the helix at the intersection point   *
*                                                                       *
*    Author  :  P. Billoir                                              *
*                                                                       *
*    First version : 26-01-88                                           *
*                                                                       *
*************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL XPOS,YPOS,ZPOS,WIDTH,LENGTH,LPHI,SINBMX,ALRPHI,COSNOR
*
      DIMENSION PARAMI(6),PARAMF(6),DER(8)
*
      PARAMETER (PI=3.141592654D00,TWOPI=PI*2.0D00,DPHIMN=1.0D-04)
*
      IERR=0
      COSNOR=1.
*
      RADI=PARAMI(1)
      PHII=PARAMI(2)
      BETA=PARAMI(5)
      SPHI=BETA+PHII
      COSP=DCOS(PHII)
      SINP=DSIN(PHII)
      COSSP=DCOS(SPHI)
      SINSP=DSIN(SPHI)
      COSTH=DCOS(PARAMI(4))
      SINTH=DSIN(PARAMI(4))
      COTTH=1./DTAN(PARAMI(4))
      RTRK=1./PARAMI(6)

*     center and squared radius of the projected helix
      CX=RADI*COSP-RTRK*SINSP
      CY=RADI*SINP+RTRK*COSSP
      RC2=RTRK**2

*     half width and half length
      HWID=.5*WIDTH
      HLEN=.5*LENGTH
      ZLMIN=ZPOS-HLEN
      ZLMAX=ZPOS+HLEN

*     cos and sin of the phi angle of the ladder
      SINLP=SIN(LPHI)
      COSLP=COS(LPHI)

*     start point and direction vector of plane (ladder segment)
      PX=XPOS+HWID*SINLP
      PY=YPOS-HWID*COSLP
      DX=-WIDTH*SINLP
      DY=WIDTH*COSLP
      DLEN2=WIDTH**2

*     the delta vector between the start point of the plane
*     and the helix center
      DELTAX=PX-CX
      DELTAY=PY-CY
      DELTAL=DELTAX**2+DELTAY**2

*     evaluate the discriminant of the second order equation
      QEB=2*(DX*DELTAX+DY*DELTAY)
      QEC=DELTAL-RC2
      DISC=QEB**2-4*DLEN2*QEC

*     intersection exists if DISC > 0
      IF(DISC.LE.0.) THEN
        IERR=1
      ELSE
*       compute both intersection points (the t value of the line X=P+t*D)
*       and check their values (a modified version of the quadratic solution
*       formula is used here, see e.g. Numerical Recipes)
        DISC=DSQRT(DISC)

        QEQ=-.5*(QEB+DSIGN(DISC,QEB))
        TVAL1=QEQ/DLEN2

        IF(TVAL1.GE.0..AND.TVAL1.LE.1.) THEN
          XF=PX+TVAL1*DX
          YF=PY+TVAL1*DY
          SINPNW=(XF-CX)/RTRK
          COSPNW=(CY-YF)/RTRK

*         cos of the angle between the normal of the ladder and the
*         tangent of the helix at the intersection point
          COSNOR=(COSLP*COSPNW*SINTH)+(SINLP*SINPNW*SINTH)
          IERR=0

*         choose only particles with outward direction (w.r.t. the ladder)
          IF(COSNOR.LT.0.) THEN
            IERR=1
          ENDIF
        ELSE
          IERR=3
        ENDIF

*       if the first intersection was not valid, try the second
        IF(IERR.EQ.3) THEN
          TVAL2=QEC/QEQ

          IF(TVAL2.GE.0..AND.TVAL2.LE.1.) THEN
            XF=PX+TVAL2*DX
            YF=PY+TVAL2*DY
            SINPNW=(XF-CX)/RTRK
            COSPNW=(CY-YF)/RTRK

            COSNOR=(COSLP*COSPNW*SINTH)+(SINLP*SINPNW*SINTH)
            IERR=0

            IF(COSNOR.LT.0.) THEN
              IERR=1
            ENDIF
          ELSE
            IERR=3
          ENDIF
        ENDIF

*       if an intersection was found, calculate the remaining values
        IF(IERR.EQ.0) THEN
          RADIF2=XF**2+YF**2
          RADIF=DSQRT(RADIF2)
          SINPHF=YF/RADIF
          COSPHF=XF/RADIF
          SINBF=(SINPNW*COSPHF)-(COSPNW*SINPHF)

*         exit if beta too large at the intersection
          IF(DABS(SINBF).GT.SINBMX) THEN
            IERR=2
          ELSE
            SPHINW=DATAN2(SINPNW,COSPNW)
            IF (SPHINW.LT.0.) SPHINW = SPHINW+TWOPI

            DPHI=SPHINW-SPHI

            ALRPHI=RTRK*DPHI

            STAPX=RADI*COSP
            STAPY=RADI*SINP
            TEST=DSQRT((STAPX-XF)**2+(STAPY-YF)**2)

*           select positive or negative extrapolations, or both
            IF(ALRPHI*IDIR.LT.0.) THEN
              IERR=1
              RETURN
            ENDIF

            ZF=PARAMI(3)+(COTTH*RTRK*DPHI)

*           exit if outside of limits
            IF(ZF.LT.ZLMIN.OR.ZF.GT.ZLMAX) THEN
              IERR=3
            ELSE
*             final parameters
              PARAMF(1)=RADIF
              PARAMF(2)=DATAN2(YF,XF)
              IF (PARAMF(2).LT.0.) PARAMF(2) = PARAMF(2)+TWOPI

              PARAMF(3)=ZF
              PARAMF(4)=PARAMI(4)
              PARAMF(5)=SPHINW-PARAMF(2)
              PARAMF(6)=PARAMI(6)

*             computation of derivatives (TODO: To be finished)
              IF(IOPT.EQ.1) THEN
                DER(1)=0.

                DER(4)=RTRK*DPHI*(1.+COTTH**2)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
*
      END
