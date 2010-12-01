
      subroutine bcalhit(pin,q,vin,B,EB,zbcal,pout,vout)
*
*    INPUT:
*      pin = momentum (px,py,pz, in GeV/c)
*      q = charge of particle
*      vin = point where pin is given (x,y,z, in mm)
*      B = B-field (in Tesla, assumed to be in the same
*       direction as the z-axis)
*      EB = beam energy (GeV)
*      zbcal = z where output is requested (in mm)
*
*    OUTPUT:
*      pout = momentum at  zbcal
*      vout = track cooridnates at at zbcal
*
* AUTHORS:
*
*    M. Berggren, DESY
*
*
* CREATION DATE:
*
*    20-JUL-2007
*
*

      real pin(*), q, vin(*),B,EB,zbcal,pout(*),vout(*)
      real p(7), v(10), zpla, iprp, ipz, phi0,pend(10),dst,vend(3),
     .       track(5)
      logical first
      integer error
*+CDE,PXCONS.
      REAL BMAG, AMPI, PI, EBEAM, EBMDEF, EPS, CONSA
* General Constants.
      COMMON /PXCONS/ BMAG, AMPI, PI, EBEAM, EBMDEF
      PARAMETER( EPS=1.0E-20)
      PARAMETER( CONSA = .0029979246 )
*-----
*-----
*      print*,"helloworld, pin = ", pin(1), pin(2), pin(3)    
      BMAG=B
      PI=4.*ATAN(1.)
      EBEAM=EB
*      EBMDEF=
*      AMPI=
*-----
      call ucopy(pin,p,3)
      p(7) = q
      call vscale(vin,0.1,v,3)
      first=.true.
      zpla=zbcal/10.

*      print*,"helloworld, p = ", p(1), p(2), p(3)

      call ZTTPLA(P,V,ZPLA,FIRST,IPRP,IPZ,PHI0,
     .                  PEND,VEND,DST,TRACK,ERROR)
      call ucopy(pend,pout,3)
      call vscale(vend,10.0,vout,3)
      if ( pin(3) .lt. 0. ) then
         vout(3)=-vout(3)
      endif
*        (also inpact parameters (IPRP,IPZ,PHI0),
*         distance traveled (DST), and helix-
*         parameters (TRACK, Delphi convention)
*         are calculated, if they should be of
*         use ...)
      
      end
   
      SUBROUTINE ZTTPLA(P,V,ZPLA,FIRST,IPRP,IPZ,PHI0,
     .                  PEND,VEND,DST,TRACK,ERROR)
*
*      SUBROUTINE ZTTPLA(P, V, ZPLA, FIRST, IPRP, IPZ, PHI0,
*                  PEND, VEND, DST, TRACK, ERROR) : Track intersect with PLAne
*      --------------------------------------------------
*
*
*        This routine tracks a particle with momentum p from v to the
*     plane surface at z coordinate ZPLA. It returns the momentum (pend)
*     and position (vend) at this point. It also returns the TKR parameters
*     at the point.
*
* ARGUMENTS:
*
*  INPUT :
*
*    P : Momentum of track
*    V  : Position where P is given
*    ZPLA : z at which track parameters are requested.
*           If ZPLA=0, the perigee parameters will be returned.
*    FIRST : TRUE if first call for this track.
*    IPRP : true IP in R-Phi (input if FIRST = FALSE)
*    IPZ : true IP in Z (input if FIRST = FALSE)
*    PHI0 : Phi angle to perigee point in R-phi plane
*            (input if FIRST = FALSE)
*
*  OUTPUT :
*
*    IPRP : true IP in R-Phi (output if FIRST = TRUE)
*    IPZ : true IP in Z (output if FIRST = TRUE)
*    PHI0 : Phi angle to perigee point in R-phi plane
*            (output if FIRST = TRUE)
*    PEND: Momentum of track (if ERROR <>0 then PEND=P)
*    VEND  : Position where PEND is given (if ERROR <>0 then VEND=V)
*    DST : Distance along track helix from V to VEND
*    TRACK : TKR parameters for the track.
*                 for the track
*    ERROR : error flag :
*             0 - no error
*             1 - no intersection in barrel.
*             2 - track turns back.
*             3 - other error.
*
* IMPLICIT INPUTS:
*
*    From /PXCONS/ : BMAG, the magnetic field; CONSA, Speed of light.
*
* IMPLICIT OUTPUTS:
*
*    none
*
* CALLED BY :
*
*    ZTCALO, Track to CALOrimeter
*    ZTTISP, TRacK InterSection with Plane
* CALLED ROUTINE(S):
*
*         ZTREFL (To REFLect track in x-y plane.)
*         ZTTRIP (To get TRue Impact Parameters.)
*         ZUI2PI (To bring angle In the range 0 to 2PI.)
*         ZUINPI (To bring angle IN the range 0 to PI.)
*         PXTPPX (To transform Track parameters on Plan to PXpypz)
*         ZTTCYL (To get Track intersection with CYLinder.)
*
*
* SIDE EFFECTS:
*
*    none
*
*
* AUTHORS:
*
*    M. Berggren, UFRJ
*
*
* CREATION DATE:
*
*    10-dec-1995
*
*
*%            C H A N G E   L O G
*%
*%       Date     | Name  | Description
*%----------------+-------+-----------------------------------------------------
*%   02-04-2001   | K.H.  | crash-prot for track along beam
*%   22-06-2001   | M.B.  | fix error in geometry giving wrong vend
*%[change_entry]
*
      IMPLICIT NONE
*+CDE,PXCONS.
      REAL BMAG, AMPI, PI, EBEAM, EBMDEF, EPS, CONSA
* General Constants.
      COMMON /PXCONS/ BMAG, AMPI, PI, EBEAM, EBMDEF
      PARAMETER( EPS=1.0E-20)
      PARAMETER( CONSA = .0029979246 )
*-----
      REAL P(*),V(*), ZPLA, IPRP,IPZ,PHI0,
     .     PEND(*),VEND(*), DST, TRACK(*)
      LOGICAL FIRST
      INTEGER ERROR

      REAL COTTH,ZDIST,COSTH,SINTH,TANTH,THETA,PINV,PMOM,
     .     BCONST,RINV,RTRK,PDST,PDSTXY,PHIROT,X,Y,XPRI,YPRI,
     .     XBIS,YBIS,Q,PHIPRI,ZSTART
      LOGICAL OTHSID

      REAL VMOD
      integer iii
      ERROR = 0
      IF ( FIRST ) THEN

*          Reflect tracks with p-z < 0 to positive side

         CALL ZTREFL(P,THETA,OTHSID)

         COSTH=COS(THETA)
         SINTH=SIN(THETA)
         IF(SINTH.EQ.0) THEN
            ERROR = 3
            RETURN
         ELSE
            COTTH=COSTH/SINTH
            TANTH=1.0/COTTH
         ENDIF

*          signed inverse momentum,signed radius and signed curvature
         PMOM=-VMOD(P,3)*SIGN(1.0,P(7))

         BCONST=CONSA*BMAG
         PINV=1./PMOM
         RINV=BCONST*PINV/SINTH
         RTRK=1.0/RINV

*          get true ip implied by p and v
         CALL ZTTRIP(P,V,IPRP,IPZ,PHI0)

         ZSTART = V(3)
      ENDIF


      Q=P(7)
      IF ( ZPLA .NE. 0.0) THEN
*
*          Distance in z from xy-perigee to plane
         IF ( OTHSID ) THEN
            ZDIST = -ZPLA+IPZ
         ELSE
            ZDIST = ZPLA-IPZ
         ENDIF

*          Distance traveled along helix from perigee to ZPLA
         PDST = ZDIST/COSTH
*          ... and projexted in XY plane
         PDSTXY = ZDIST*TANTH

*          Rotation of track from perigee to ZPLA
         IF ( Q .LT. 0 ) THEN
            PHIPRI = ABS(PDSTXY/RTRK)
         ELSE
*             (this isn't real the rotation (should be 2pi),
*              but the right angle for the calculations below)
            PHIPRI = PI - ABS(PDSTXY/RTRK)
         ENDIF

*          Bring with in [0,2PI] (the track can have done many
*          turns before reaching ZPLA!)
         CALL ZUI2PI(PHIPRI,PHIPRI)

*          X^prime and Y^prime of intersection in system with origin
*          at the centre of the track helix, with the perigee in the xy-
*          plane on the x^prime axis, and with the y^prime axis oriented
*          so that p(Y^prime) > 0 at the perigee (p(X^prime) is ofcourse
*          0 there).
         XPRI = ABS(RTRK)*COS(PHIPRI)
         YPRI = ABS(RTRK)*SIN(PHIPRI)

*          X^bis and Y^bis of intersection: same as the prime coordinates,
*          except that the origin is moved to the detector origin (This
*          is a translation in x only, since poth the centre of the helix
*          and the perigee is on the x^prime axis !)
         XBIS = XPRI + Q*ABS(RTRK)+IPRP
         YBIS = YPRI

*          X and Y of intersection (ie. detector coordinates) : rotate
*          by the angle to the perigee.
*         IF ( IPRP .GT. 0.0 ) THEN
            PHIROT = PHI0
*         ELSE
*            PHIROT = PHI0 - PI/2.0
*         ENDIF
         X = XBIS * COS(PHIROT) - YBIS * SIN(PHIROT)
         Y = XBIS * SIN(PHIROT) + YBIS * COS(PHIROT)


*             Make TKR parameters
         TRACK(2) = SQRT(X**2+Y**2)
         IF(X.EQ.0. .AND. Y.EQ.0.) THEN
            TRACK(1) = 0.
         ELSE
            TRACK(1) = TRACK(2)*ATAN2(Y,X)
         ENDIF
         IF ( OTHSID ) THEN
            TRACK(3) = PI - THETA
         ELSE
            TRACK(3) = THETA
         ENDIF
         CALL ZUINPI(TRACK(3),TRACK(3))
         CALL ZUI2PI(PHI0 - Q*PI/2.0 + PHIPRI ,TRACK(4) )
         TRACK(5) = 1.0/PMOM

*             get momentum and crossing point
         CALL PXTPPX(ZPLA,TRACK,PEND,VEND)
         PEND(4) = P(4)
         PEND(5) = P(5)

         IF ( OTHSID ) THEN
            ZDIST = -ZPLA-ZSTART
         ELSE
            ZDIST =  ZPLA-ZSTART
         ENDIF

*          Distance traveled along helix from V to ZPLA
         DST = ABS(ZDIST/COSTH)
      ELSE

*             ZPLA = 0 -> perigee parameters requested.
*             Make TKR parameters at perigee.
*--         CALL ZTTCYL(P,V,ZPLA,FIRST,IPRP,IPZ,PHI0,
*--     .                  PEND,VEND,DST,TRACK,ERROR)
      ENDIF
      END
      SUBROUTINE ZTREFL(P,THT,OTHSID)

*      SUBROUTINE ZTREFL(P, THT, OTHSID)  :    REFLect a track in the x-y plane
*      --------------------------------------------------
*
*
*      If the track points to negative z, reflect it, otherwise don't.
*
* ARGUMENTS:
*
*  INPUT :
*
*    P : momentum of track.
*  OUTPUT :
*
*    THT : reflected theta.
*    OTHSID : TRUE if track was reflected.
*
* IMPLICIT INPUTS:
*
*    none.
*
* IMPLICIT OUTPUTS:
*
*    none.
*
* CALLED BY :
*
*    ZTTRAK, TRAcK a particle through the detector
*    ZTTCYL, Track intersect with CYLinder
*    ZTTPLA, Track intersect with PLAne
* CALLED ROUTINE(S):
*
*    none.
*
*
* SIDE EFFECTS:
*
*    none.
*
* AUTHORS:
*
*    M. Berggren, UFRJ
*
*
* CREATION DATE:
*
*    07-APR-1995
*
*
*%            C H A N G E   L O G
*%
*%       Date     | Name  | Description
*%----------------+-------+-----------------------------------------------------
*%[change_entry]
*

      IMPLICIT NONE
      REAL P(*), THT
      LOGICAL OTHSID

       REAL ABSP, TH,PHI,PI

       REAL VMOD

       PI = 4.0*ATAN(1.0)
       ABSP = VMOD(P,3)

*           Get Theta and phi of the track, then
*           modify theta to be between 0 and pi/2
          TH = ACOS ( P(3)/ABSP)
          THT = ABS(TH)
          PHI = ATAN2(P(2),P(1))
17        THT = ABS(THT)
          IF ( THT .GT. PI/2.0 ) THEN
             THT = PI - THT
             GOTO 17
          ENDIF
          IF ( TH .NE. THT ) THEN
             OTHSID = .TRUE.
          ELSE
             OTHSID = .FALSE.
          ENDIF

      END
      SUBROUTINE ZTTRIP(P,V,TRIPRP,TIPZ,PHIIP)

*
*
*      SUBROUTINE ZTTRIP(P, V, TRIPRP, TIPZ, PHIIP)  :
*                 calculate Impact Parameter
*      --------------------------------------------------
*
*
*     This routine calculates the impact-parameter in Rphi and z, with
*     geometrical sign, from the momentum of a track (p) and a point on its
*     trajectory (v). It also gives the phi angle of the point of closest
*     approach to the origin of the projection of the track in the x-y
*     (i.e. R phi) plane.
*
* ARGUMENTS:
*
*  INPUT :
*
*    P : 3-momentum of the track (in GeV)
*    V : The point where P is given (in cm).
*
*  OUTPUT :
*
*    TRIPRP : Impact parameter in Rphi (cm)
*    TIPZ : Impact parameter in z (cm)
*    PHIIP : Phi angle of point of closest approach (Radians)
*    All of these will be set to zero for tracks along the z-axis
*
* IMPLICIT INPUTS:
*
*    From /PXCONS/ : BMAG, the magnetic field; CONSA, Speed of light.
*
* IMPLICIT OUTPUTS:
*
*    none.
*
* CALLED BY :
*
*    ZTTCYL, Track intersect with CYLinder
*    ZTTPLA, Track intersect with PLAne
*    ZADECD, DECoDe data
* CALLED ROUTINE(S):
*
*         VDOT,VMOD,VUNIT (CERNLIB routines)
*
*
*
* SIDE EFFECTS:
*
*    none
*
*
* AUTHORS:
*
*    M. Berggren, LAL
*
*
* CREATION DATE:
*
*    27-jul-1992
*
*
*%            C H A N G E   L O G
*%
*%       Date     | Name  | Description
*%----------------+-------+-----------------------------------------------------
*%[change_entry]
*
*
       IMPLICIT NONE
       REAL P(*),V(*),TRIPRP,TIPZ,PHIIP

*+CDE,PXCONS.
      REAL BMAG, AMPI, PI, EBEAM, EBMDEF, EPS, CONSA
* General Constants.
      COMMON /PXCONS/ BMAG, AMPI, PI, EBEAM, EBMDEF
      PARAMETER( EPS=1.0E-20)
      PARAMETER( CONSA = .0029979246 )
*----
       REAL RV(3),ABSPT,COSPHI,SINPHI,IPSIGN,PHITR
       REAL*8 C,BFIELD,RHO,RRR,T2,TZIMP,SP,CAPR,Q ,ANGC
*       INTEGER LLL

       REAL VDOT,VMOD

       PI = 4.0*ATAN(1.0)

       IF ( P(1) .EQ. 0.0 .AND. P(2) .EQ. 0.0 ) THEN
*           Track right down the z-axis ....
          TRIPRP = 0.0
          TIPZ = 0.0
          PHIIP = 0.0
          RETURN
       ENDIF
       ABSPT = VMOD(P,2)
       RRR =0.0D0
*         charge
       Q = DBLE(P(7))
*           Particle direction
       CALL VUNIT(P,RV,3)

       IF ( Q .EQ. 0.0D0 ) THEN
          TRIPRP = 0.0
          TIPZ =0.0
*          CALL VSCALE(RV,-VDOT(V,RV,2),LPROJ,2)
*          CALL VADD(LPROJ,V,LCOMPL,2)
*          TRIPRP = VMOD(LCOMPL,2)
*          TIPZ =
       ELSE

*           Curvature of track
       C = DBLE(CONSA)
       BFIELD = DBLE(BMAG)
       RHO = DBLE(ABSPT)/(C*BFIELD)

*            Calculate impact parameter in  R-Phi

*              R distance to reference point
       RRR = DBLE(VMOD(V,2))

*              Phi angle at reference point
       SINPHI = RV(2)/VMOD(RV,2)
       COSPHI = RV(1)/VMOD(RV,2)
       PHITR = ATAN2(RV(2),RV(1))


*           Correct IP when B-field is present (geom. sign)
       T2 = DBLE(V(1)*SINPHI-V(2)*COSPHI)
       T2 = T2*Q
       IF ( DABS((RRR/RHO)**2 +
     .                2*T2/RHO) .LT. 1.0E-6) THEN
          TRIPRP = SNGL(-(0.5* (RRR**2/RHO) + T2 ))
       ELSE
          TRIPRP = SNGL(-(RHO - RHO * DSQRT ( 1D0 + (RRR/RHO)**2 +
     .                2D0*T2/RHO )))
       ENDIF

*           sign the ip
       IPSIGN = ((V(2) - P(1)*SNGL(Q*RHO)/VMOD(P,2))**2 +
     .                   (V(1) + P(2)*SNGL(Q*RHO)/VMOD(P,2))**2 -
     .                        SNGL(RHO)**2)
       TRIPRP = SIGN(TRIPRP,SNGL(Q)*IPSIGN)

*        Angle to point of closest approach:
*        Use sine theorem to find the amount of rotation between
*        the reference point and the perigee. This gives the
*        direction of the track at the perigee; the angle to the
*        perige is 90 deg less than that.
*        the triangle has sides a=q*rho, b=q*rho+triprp
*        c=rrr, angles A,B and C
*
*        Half the circumference of the cicumscribed circle
       SP = (RRR+Q*RHO+(Q*RHO+TRIPRP))/2.0D0
       IF (( SP+(-Q*RHO-TRIPRP) .EQ. 0.0 ) .OR.
     .       ( SP-Q*RHO .EQ. 0.0 ) ) THEN
*           Reference radius was to perigee...
          PHIIP = PHITR - PI/2.0
       ELSE
*            Radius of the cicumscribed circle
          CAPR = -RRR*Q*RHO*(-Q*RHO-TRIPRP)/
     .    (4D0*DSQRT(ABS(SP*(SP-Q*RHO)*(SP-RRR)*(SP+(-Q*RHO-TRIPRP)))))

          IF (RHO**2+(-Q*RHO-TRIPRP)**2 .GE. RRR**2) THEN
             ANGC = DASIN(Q*ABS(RRR)/(2.0D0*CAPR))
          ELSE
*              Track rotated by more than 90 degrees...
             ANGC = DBLE(PI) - DASIN(Q*ABS(RRR)/(2.0D0*CAPR))
          ENDIF
*            rotate track by C, then subtract pi/2
          PHIIP = SNGL(DBLE(PHITR) +
     .           ANGC - DBLE(PI/2.0))
*     .           angC - dsign(1.0d0,triprp)*dble(pi/2.0))
       ENDIF
       IF ( RRR .EQ. 0.0D1 ) THEN
          TIPZ = 0.0
       ELSEIF ( ABS(VDOT(RV,V,2)/((Q*RHO+TRIPRP) * VMOD(RV,2))) .LT.
     .             1.0E-6 ) THEN
          TZIMP = RHO * DBLE(VDOT(RV,V,2))/
     .                 (ABS(DBLE(Q)*RHO+TRIPRP) * DBLE(VDOT(RV,RV,2)))
          TIPZ = SNGL(DBLE(V(3)) - DBLE(RV(3))*TZIMP)
       ELSE
          IF ( ABS(VDOT(RV,V,2)/((Q*RHO+TRIPRP) * VMOD(RV,2))) .GT.
     .             1.0 ) THEN
*                PRINT *, ' EVENT : ',I
*             PRINT *, ' PROBLEM IN TRIP : ARGUMENT TO ASIN = ',
*     .       VDOT(RV,V,2)/((RHO+TRIPRP) * VMOD(RV,2))
*             PRINT *, ' RV = ', RV
*             PRINT *, '  V = ',  (V(LLL),LLL=1,3)
*             PRINT *, ' RHO AND TRIPRP = ',RHO,TRIPRP
*             PRINT *, ' ABSPT = ',ABSPT
             TIPZ = -1000.0
          ELSE
             TZIMP = RHO *
     .                DASIN( DBLE(VDOT(RV,V,2)) /
     .                  ( ABS(DBLE(Q)*RHO+TRIPRP)*DBLE(VMOD(RV,2))))/
     .                VMOD(RV,2)
             TIPZ = SNGL(DBLE(V(3)) - DBLE(RV(3))*TZIMP)
          ENDIF
       ENDIF
       ENDIF
      END
      SUBROUTINE ZUI2PI(PHIIN,PHIOUT)
*
*
*      SUBROUTINE ZUI2PI(PHIIN, PHIOUT)  : bring angle Into the
*                 interval 0 to 2PI
*      --------------------------------------------------
*
*
*         The angle PHIIN is set to the equivalent angle in
*         0 to 2pi.
*
* ARGUMENTS:
*
*  INPUT :
*
*    PHIIN : input angle.
*  OUTPUT :
*
*    PHIOUT : output angle = PHIIN, but given in the interval 0-2 pi
*
* IMPLICIT INPUTS:
*
*    none.
*
* IMPLICIT OUTPUTS:
*
*    none.
*
* CALLED BY :
*
*    any routine.
*
* CALLED ROUTINE(S):
*
*    none.
*
*
* SIDE EFFECTS:
*
*    none.
*
* AUTHORS:
*
*    K. Hultqvist, Stockholms Universitet
*
*
* CREATION DATE:
*
*    20-APR-2001
*
*
*%            C H A N G E   L O G
*%
*%       Date     | Name  | Description
*%----------------+-------+-----------------------------------------------------
*%[change_entry]
*
*
      IMPLICIT NONE
      REAL PHIIN,PHIOUT
      REAL TWOPI
      TWOPI = 8.0*ATAN(1.0)
      PHIOUT = PHIIN
1     CONTINUE
      IF ( PHIOUT .GE. TWOPI ) THEN
         PHIOUT = PHIOUT - TWOPI
         GOTO 1
      ELSEIF ( PHIOUT .LT. 0.0 ) THEN
         PHIOUT = PHIOUT + TWOPI
         GOTO 1
      ENDIF
      END
      SUBROUTINE ZUINPI(PHIIN,PHIOUT)
*
*
*      SUBROUTINE ZUINPI(PHIIN, PHIOUT)  : bring angle INto the interval 0 to PI
*      --------------------------------------------------
*
*
*         The polar angle PHIIN is set to the equivalent angle in
*         0 to pi (-PHIIN = PHIN, since that corresponds to a
*         rotation around the polar axis).
*
* ARGUMENTS:
*
*  INPUT :
*
*    PHIIN : input angle.
*  OUTPUT :
*
*    PHIOUT : output angle = PHIIN, but given in the interval 0-pi
*
* IMPLICIT INPUTS:
*
*    none.
*
* IMPLICIT OUTPUTS:
*
*    none.
*
* CALLED BY :
*
*    any routine.
*
* CALLED ROUTINE(S):
*
*    none.
*
*
* SIDE EFFECTS:
*
*    none.
*
* AUTHORS:
*
*    M. Berggren, UFRJ
*
*
* CREATION DATE:
*
*    11-APR-1995
*
*
*%            C H A N G E   L O G
*%
*%       Date     | Name  | Description
*%----------------+-------+-----------------------------------------------------
*%[change_entry]
*
      IMPLICIT NONE
      REAL PHIIN,PHIOUT
      REAL PI
      PI = 4.0*ATAN(1.0)
      PHIOUT = PHIIN
1     CONTINUE
      IF ( PHIOUT .GT. PI ) THEN
         PHIOUT = 2*PI - PHIOUT
         GOTO 1
      ELSEIF ( PHIOUT .LT. 0.0 ) THEN
         PHIOUT = -PHIOUT
         GOTO 1
      ENDIF
      END
      SUBROUTINE PXTPPX(Z,PAR,P,X)
*
*
*      SUBROUTINE PXTPPX(Z, PAR, P, X) :
*            transform Track parameters on Plan to PXpypz
*      --------------------------------------------------
*
*
*       This routine transforms the  TKR parameters given in PAR,
*       at Z to particle 3 momentum and the cartesian coordinates
*       of the point where PAR was given.
*
* ARGUMENTS:
*
*  INPUT :
*
*    Z : Z where track parameters are given. If Z = 0,
*             the parameters are assumed to be at the perigee.
*    PAR : Track parameters.
*  OUTPUT :
*
*    P : 3-momentum of the track.
*    X : Cartesian coordinates of the point where PAR was given.
*
* IMPLICIT INPUTS:
*
*    From /PXCONS/ : BMAG, magnetic field; CONSA, speed of light.
*
* IMPLICIT OUTPUTS:
*
*    none.
*
* CALLED BY :
*
*    any routine.
*
* CALLED ROUTINE(S):
*
*         ZUI2PI (To bring an angle In the range 0 to 2PI.)
*
*
* SIDE EFFECTS:
*
*    none.
*
* AUTHORS:
*
*    M. Berggren, UFRJ
*
*
* CREATION DATE:
*
*    12-APR-1995
*
*
*%            C H A N G E   L O G
*%
*%       Date     | Name  | Description
*%----------------+-------+-----------------------------------------------------
*%    01-04-20    | K.H.  | crash prot for track along z
*%[change_entry]
*

      IMPLICIT NONE

*+CDE,PXCONS.
      REAL BMAG, AMPI, PI, EBEAM, EBMDEF, EPS, CONSA
* General Constants.
      COMMON /PXCONS/ BMAG, AMPI, PI, EBEAM, EBMDEF
      PARAMETER( EPS=1.0E-20)
      PARAMETER( CONSA = .0029979246 )
*-----
      REAL Z,PAR(*),P(*),X(*)
      REAL PT,PHI,BZ,PMOM,RADIUS

      BZ = BMAG*CONSA


      IF ( Z .EQ. 0.0 ) THEN
*           Perige params
         RADIUS = Z
         PT = BZ/ABS(PAR(5))
         PHI = PAR(4)-PI/2.0
         CALL ZUI2PI(PHI,PHI)
         P(1) = ABS(PT)*COS(PAR(4))
         P(2) = ABS(PT)*SIN(PAR(4))
         P(3) = (1.0/TAN(PAR(3)))*ABS(PT)
         X(1) = PAR(1)*COS(PHI)
         X(2) = PAR(1)*SIN(PHI)
         X(3) = PAR(2)
      ELSE
*           track params
         IF(PAR(2) .NE. 0.) THEN
            PHI = PAR(1)/PAR(2)
         ELSE
            PHI = 0.
         ENDIF
         PMOM = 1.0/ABS(PAR(5))
         P(1) = SIN(PAR(3))*ABS(PMOM)*COS(PAR(4))
         P(2) = SIN(PAR(3))*ABS(PMOM)*SIN(PAR(4))
         P(3) = COS(PAR(3))*ABS(PMOM)
         X(1) = PAR(2)*COS(PHI)
         X(2) = PAR(2)*SIN(PHI)
         X(3) = Z
      ENDIF
      END
