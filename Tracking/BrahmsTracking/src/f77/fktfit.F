*  
* $Id: fktfit.F,v 1.3 2007-11-22 09:26:31 rasp Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2006/11/02 11:49:00  rasp
*
* Modified files adopted to FullLDCTracking processor
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
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
***********************************************************************
*                                                                     *
      SUBROUTINE FKTFIT(IDEB,IHIS,IOUTR,IRNKXY,IRNKZ,
     & NTEFIT,IPATTB,IPATTF,ISREFB,ISREFF,IERB,IERF,IERX,IDTER)
*                                                                     *
***********************************************************************
*                                                                     *
*   Steering for single track fit                                     *
*   Compute intersection of reference trajectory with measurement     *
*     surfaces, material surfaces and, if requested, extrapolation    *
*     surfaces                                                        *
*   Sort the intersections along the trajectory                       *
*   Perform a double Kalman filter (forward and backward)             *
*                                                                     *
*   Input  :  IDEB    :  level of debug for this event                *
*             IHIS    :  flag for histograms                          *
*             IOUT = MOD(IOUTR,10)
*             IOUT    :  if 0 skip outlayers logic                    *
*                     :  if 1 no ITC rejection                        *
*                     :  if 2 no TPC rejection                        *
*                     :  if 3 anything may be rejected                *
*                     :  if 5 either ITC or TPC must survive          *
*                     :  if 6 anything not decreasing rank may be     *
*                                                  rejected           *
*                     :  if 8 no FTD rejection                        *
*             ITAN = MOD(IOUTR,100)/10                                *
*             ITAN    : if 0 called from track search                 *
*             ITAN    : if > 0 called from track fit,                 *
*                       ITAN = number of iteration (1,2,...)          *
*             ICHI = IOUTR/100                                        *
*             ICHI    : if 0 incl. internal chi^2 and n.d.f. (normal) *
*             ICHI    : if 1 no internal chi^2 and n.d.f. (search)    *
*                                                                     *
*     Output :  NTEFIT  :  number of TE's included in backward fit    *
*             IPATTB  :  bit pattern of detectors (backward fit)      *
*             IPATTF  :  bit pattern of detectors (forward fit)       *
*             ISREFB  :  serial number of reference surface for       *
*                        parameters (backward fit)                    *
*             ISREFF  :  serial number of reference surface for       *
*                        parameters (forward fit)                     *
*             IERB    :  error flag for backward fit :                *
*                             0 if OK                                 *
*                             1 if insufficient information in TER's  *
*                             2 if weight matrix not positive         *
*                             3 if abnormal topology                  *
*             IERF    :  error flag for forward fit (same meaning)    *
*                                                                     *
***********************************************************************
#include "include/fksave.inc"
#include "include/fkclun.inc"
#include "include/fsfparams.inc"
#include "include/fkcom.inc"
#include "include/fkcom_plus.inc"
#include "include/fkcons.inc"
#include "include/fkddes.inc"
#include "include/fkddes1.inc"
#include "include/fkddes2.inc"
#include "include/fkexts.inc"
#include "include/fkfild.inc"


      real*4    sinbmx          ! maximum value for sin(beta) for fkxcyl
      data sinbmx / 0.999d0 /
      REAL      FDEDX,DFDEDX

      CHARACTER*6 FBTEXT(0:1), FLTEXT
      DATA FBTEXT /'forwd.','barrel'/
*
      DOUBLE PRECISION PI,TWOPI,PARINT(6),DER(8),
     /DQINT(5),WGINT(15),REMIN,REMAX,ZEMIN,ZEMAX,RINTER,ZINTER
      DOUBLE PRECISION SIGNP,STH,PABS,ETOT,DENERG
*
      PARAMETER (PI=3.14159265359D00,TWOPI=PI*2.0D00)
*
      DIMENSION IERX(NEXSMX),PROBS(NTEMX),IDTER(NTEMX)
*
      CHARACTER*11 TXTDQ(5)
      CHARACTER*4 TYPOLD
*
      real*8 amdedx
      data amdedx / 0.139d0 /
*
      dimension pchint(ntemx),ipchint(ntemx),iidet(ntemx)
      EXTERNAL FXPROB   ! MVL

      Real PhiMom,ThetaMom
      Real PhiSurface,ThetaSurface
      Real TanThetaSurface
      Real unitSurface(3)
      Real unitMom(3)
      Real xxel,xxrl



      IOUT=MOD(IOUTR,10)
      ITAN=MOD(IOUTR,100)/10
      ICHI=MOD(IOUTR/100,10)

C(KH)
      ideb_save=ideb
C(KH)

      ISCALOD=0
      ISCALTPC=0

      KHIST=LHIST
*
      IF(IDEB.GE.1) WRITE(IPRNT,'(/,''     ** entering FKTFIT'','//
     $     'I5,'' TE(s)   detect.'',20I4)') NTE,(IDET(I),I=1,NTE)
*
      IERB=0
      NSURF=0
      CALL VZERO(INUMB,NSRFMX)
      CALL VZERO(XRLS,NSRFMX)
      call vzero(xeloss,nsrfmx)
      CALL FKDPV0(DQTER,5*NTEMX)
*
      IF(PARREF(4).LT.PI/2.) THEN
        SIDPL=1.
      ELSE
        SIDPL=-1.
      ENDIF
*
*   initialize r and z limits
*
      REMIN=PARREF(1)
      REMAX=REMIN
      ZEMIN=PARREF(3)
      ZEMAX=ZEMIN
*
*   compute intersections with reference surface of all TER banks
*
      DO 10 ITE=1,NTE
*
        IF(IERTE(ITE).EQ.0) THEN
          IF(TYPE(ITE).EQ.'CYLI') THEN
             SPPAR=PARAM(1,ITE)
             CALL FKXCYL(PARREF,0,SPPAR,-500.,500.,sinbmx,0,
     &       PARINT,DER,ALRPH,IER)
            IF(IER.EQ.0) THEN
              RINTER=PARINT(1)
              ZINTER=PARINT(3)
            ENDIF
          ELSE
             SPPAR=PARAM(1,ITE)
             CALL FKXPLA(PARPLA,0,SPPAR,0.,500.,0,
     &       PARINT,DER,ALRPH,IER)
            IF(IER.EQ.0) THEN
              RINTER=DSQRT(PARINT(2)**2+PARINT(3)**2)
              ZINTER=PARINT(1)
            ENDIF
          ENDIF
          IF(IER.EQ.0) THEN
            NSURF=NSURF+1
            INUMB(NSURF)=ITE
            ALRPHI(NSURF)=ALRPH
            TYPSUR(NSURF)=TYPE(ITE)
            CALL FKDPCO(PARINT,PARAMS(1,NSURF),6)
*
*   update limits in r and z if needed
            IF(RINTER.LT.REMIN) REMIN=RINTER
            IF(RINTER.GT.REMAX) REMAX=RINTER
            IF(ZINTER.LT.ZEMIN) ZEMIN=ZINTER
            IF(ZINTER.GT.ZEMAX) ZEMAX=ZINTER
*
          ELSE
            IERTE(ITE)=4
            IF(IDEB.GE.1) WRITE(IPRNT,
     $           '(1X,''FKTFIT(TKF) TE'',I3,'' : error'','//
     $           'I3,'' in FKXCYL or FKXPLA'')') ITE,IER
*
          ENDIF
        ENDIF
   10 CONTINUE
*
*   if no TE accepted : set error flag and return
*
      IF(NSURF.EQ.0) THEN
C(KH)---begin
        IF (IDEB.GT.0) THEN
           WRITE(IPRNT,'(7X,''FKTFIT(TKF) : no TE accepted'')')
           WRITE(6,*) 'TE information for this failure:'
           DO ITE=1,NTE
              WRITE(6,*) 'TE ',ITE,': ',TYPE(ITE),'; ierte=',IERTE(ITE)
              WRITE(6,*) '  parameters: ',PARAM(1,ITE),PARAM(2,ITE),
     >             PARAM(3,ITE),PARAM(4,ITE),PARAM(5,ITE),PARAM(6,ITE)
           END DO
        END IF
C(KH)---end
        IERB=3
        IERF=3
        RETURN
      ENDIF
*
*   compute intersections with surfaces where extrapolation/interpolation
*   is requested
*
      ALRHPC=0.
      DO 20 IX=1,NEXS
        IF(TYPX(IX).EQ.'CYLI') THEN
          CALL FKXCYL(PARREF,0,RZSURF(IX),ZRMIN(IX),ZRMAX(IX),sinbmx,0,
     &    PARAMX(1,IX),DER,ALRPH,IERX(IX))
          IF(IERX(IX).EQ.0) THEN
            RINTER=PARAMX(1,IX)
            ZINTER=PARAMX(3,IX)
            IF(IX.EQ.NEXHPC) ALRHPC=ALRPH
          ENDIF
        ELSE
          CALL FKXPLA(PARPLA,0,SIDPL*RZSURF(IX),ZRMIN(IX),ZRMAX(IX),0,
     &    PARAMX(1,IX),DER,ALRPH,IERX(IX))
          IF(IERX(IX).EQ.0) THEN
            RINTER=DSQRT(PARAMX(2,IX)**2+PARAMX(3,IX)**2)
            ZINTER=PARAMX(1,IX)
*   forward tracks: ALRHPC=0! no extr. at EMF only if barrel track (SD)
            IF((ALRHPC.NE.0.).AND.(ALRPH.GT.ALRHPC)) IERX(IX)=4
*   forward tracks: limit to 180 deg turn angle
            IF(ALRPH.GT.ABS(PI/PARAMX(6,IX))) IERX(IX)=5
          ENDIF
        ENDIF
        IF(IERX(IX).EQ.0) THEN
          IF(IDEB.GE.2) WRITE(IPRNT,
     $         '(''       extr. surf.'',I4,2X,A4,''  param.'','//
     $         'T39,3F11.3,2F11.4,F11.6,''  proj. dist.'',F8.2)')
     $         IX,TYPX(IX),(PARAMX(J,IX),J=1,6),ALRPH
*
          NSURF=NSURF+1
          INUMB(NSURF)=-IX
          ALRPHI(NSURF)=ALRPH
          TYPSUR(NSURF)=TYPX(IX)
          CALL FKDPCO(PARAMX(1,IX),PARAMS(1,NSURF),6)
*
*   update limits in r and z if needed
          IF(RINTER.LT.REMIN) REMIN=RINTER
          IF(RINTER.GT.REMAX) REMAX=RINTER
          IF(ZINTER.LT.ZEMIN) ZEMIN=ZINTER
          IF(ZINTER.GT.ZEMAX) ZEMAX=ZINTER
        ENDIF
   20 CONTINUE
*
*   for the time being: distinguish betweeen barrel and forward tracks
*   using ALRHPC (SD) --> overruled by JW's FCA/FCB tests after FKKALM calls
*
      IF (ALRHPC.EQ.0.)  THEN
        IBAFO=0
      ELSE
        IBAFO=1
      ENDIF
*
*   store largest extrapolation length on each side
*
      ALRMAX=VMAX(ALRPHI,NSURF)
      ALRMIN=VMIN(ALRPHI,NSURF)
*
*   compute intersections with material surfaces
*
*   cylinders
*
c     write(*,*) ' '
c      write(*,*) 'Number of cyllindrical detectors = ',ncmat
      do 30 ic=1,ncmat
*
*   check the limits in r and z
*
        if (rcmat(ic).ge.remin .and. rcmat(ic).le.remax .and.
     $      zcmax(ic).ge.zemin .and. zcmin(ic).le.zemax) then
*
          call fkxcyl(parref,0,rcmat(ic),zcmin(ic),zcmax(ic),sinbmx,
     $         0,parint,der,alrph,ier)

          if (ier.eq.0 .and. alrph.ge.alrmin.and.alrph.le.alrmax) then
*
            if (ideb.ge.2) write(IPRNT,
     $           '(''       mat. cyl. '',I3,''  param.'','//
     $           'T39,3F11.3,2F11.4,F11.6,''  proj. dist.'',F8.2)')
     $           ic,parint,alrph

            nsurf = nsurf+1

            alrphi(nsurf) = alrph
            typsur(nsurf) = 'CYLI'
            xrls(nsurf)   = xrlc(ic)/abs(sin(parint(4))*dcos(parint(5)))
            call ucopy(parint,params(1,nsurf),2*6)
            xeloss(nsurf) =
     $           xelosc(ic)/abs(sin(parint(4))*cos(parint(5)))

          endif
        endif
   30 continue

*
*  ladders
*
c      write(*,*) ' '
c      write(*,*) 'Number of ladder detectors = ',nplmat

      do 32 ic=1,nplmat
          call fkxpln(parref,0,xplmat(ic),yplmat(ic),zplmat(ic),
     &         widplmat(ic),lenplmat(ic),phiplmat(ic),sinbmx,
     &         0,parint,der,alrph,cosn,ier)

           if (ier.eq.0 .and. alrph.ge.alrmin.and.alrph.le.alrmax) then

             if (ideb.ge.2) write(IPRNT,
     $           '(''       mat. cyl. '',I3,''  param.'','//
     $           'T39,3F11.3,2F11.4,F11.6,''  proj. dist.'',F8.2)')
     $           ic,parint,alrph

c             write(*,*) 'Final Ladder : ',parint(1)*cos(parint(2)),
c     &        parint(1)*sin(parint(2)),parint(3)

             nsurf = nsurf+1

             alrphi(nsurf) = alrph
             typsur(nsurf) = 'CYLI'
             xrls(nsurf)   = xrlpl(ic)/cosn
             call ucopy(parint,params(1,nsurf),2*6)
             xeloss(nsurf) = xelospl(ic)/cosn

           endif
 32   continue   

*
*   cones
*
c$$$      write(*,*) ' '
c$$$      write(*,*) 'Number of conical detectrors = ',nconmat
      do 31 ic=1,nconmat
*
*   check the limits in r and z
*
        rconmin = min(r1conmat(ic),r2conmat(ic))
        rconmax = max(r1conmat(ic),r2conmat(ic))
        zconmin = min(z1conmat(ic),z2conmat(ic))
        zconmax = max(z1conmat(ic),z2conmat(ic))
c$$$        write(*,*) ' '
c$$$        write(*,*) 'Conical detector ',ic,
c$$$     &             ' R1 = ',r1conmat(ic),' R2 = ',r2conmat(ic),
c$$$     &             ' Z1 = ',z1conmat(ic),' Z2 = ',z2conmat(ic)   
c$$$        write(*,*) ic,remin,remax,zemin,zemax
c$$$        write(*,*) rconmin,rconmax,zconmin,zconmax
c$$$        write(*,*) ' '
        if (rconmin.ge.remin .and. rconmax.le.remax .and.
     $      zconmin.ge.zemin .and. zconmax.le.zemax) then
*
           call fkxcon(parref,0,z1conmat(ic),z2conmat(ic),
     $                         r1conmat(ic),r2conmat(ic),
     $                sinbmx,0,parint,der,alrph,ier)
           if (ier.eq.0 .and. alrph.ge.alrmin.and.alrph.le.alrmax) then
*
              if (ideb.ge.2) write(IPRNT,
     $           '(''       mat. cyl. '',I3,''  param.'','//
     $           'T39,3F11.3,2F11.4,F11.6,''  proj. dist.'',F8.2)')
     $           ic,parint,alrph
*
              nsurf = nsurf+1
*
              alrphi(nsurf) = alrph
              typsur(nsurf) = 'CYLI'
              PhiSurface = parint(2)
              TanThetaSurface = (r2conmat(ic)-r1conmat(ic))/
     $             (z2conmat(ic)-z1conmat(ic))
              ThetaSurface = atan(TanThetaSurface) + 0.5*PI
              PhiMom = parint(5) + parint(2)
              ThetaMom = parint(4)
              unitSurface(1) = cos(PhiSurface)*sin(ThetaSurface)            
              unitSurface(2) = sin(PhiSurface)*sin(ThetaSurface)
              unitSurface(3) = cos(ThetaSurface)
              unitMom(1) = cos(PhiMom)*sin(ThetaMom)            
              unitMom(2) = sin(PhiMom)*sin(ThetaMom)
              unitMom(3) = cos(ThetaMom)
              
              cosalpha = 0.0
              do icomp=1,3
                 cosalpha = cosalpha + unitSurface(icomp)*unitMom(icomp)
              enddo

              xxrl = xrl1con(ic) +
     &             (parint(3)-z1conmat(ic))*(xrl2con(ic)-xrl1con(ic))/
     &             (z2conmat(ic)-z1conmat(ic))
              xxel = xel1con(ic) +
     &             (parint(3)-z1conmat(ic))*(xel2con(ic)-xel1con(ic))/
     &             (z2conmat(ic)-z1conmat(ic))
              
              xrls(nsurf)   = xxrl/abs(cosalpha)
              call ucopy(parint,params(1,nsurf),2*6)
              xeloss(nsurf) =
     $             xxel/abs(cosalpha)
           else
c$$$              WRITE(*,*) 'FKXCON : IER = ',IER 
           endif
        else
c$$$           WRITE(*,*) ' OUTSIDE REGION: RMIN = ',remin,
c$$$     $                               ' RMAX = ',remax,
c$$$     $                               ' ZMIN = ',zemin,
c$$$     $                               ' ZMAX = ',zemax
                    
        endif
 31   continue

*
*   planes
*
c$$$      write(*,*) '  '
c$$$      write(*,*) 'Number of planar detectors = ',npmat
c$$$      WRITE(*,*) REMIN,REMAX,ZEMIN,ZEMAX
c$$$      WRITE(*,*) PARPLA(1),PARPLA(2),PARPLA(3),PARPLA(4),PARPLA(5),
c$$$     &           PARPLA(6)
      DO 40 IP=1,NPMAT

        if ( ip .ge. npmamx ) then
          write(*,*) ' FTKFIT: Warning: too many extrapolation
     &      surfaces, ship to end! '
          goto 40
        end if
*
*   check the limits in r and z
*
        IF(SIDPL*ZPMAT(IP).GT.0.)THEN
         IF(RPMAX(IP).GE.REMIN .AND. RPMIN(IP).LE.REMAX .AND.
     &      ZPMAT(IP).GE.ZEMIN .AND. ZPMAT(IP).LE.ZEMAX) THEN
*
          
         CALL FKXPLA(PARPLA,0,ZPMAT(IP),RPMIN(IP),RPMAX(IP),0,
     &     PARINT,DER,ALRPH,IER)
          IF(IER.EQ.0 .AND. ALRPH.GE.ALRMIN.AND.ALRPH.LE.ALRMAX) THEN
*
c$$$             WRITE(*,75) IP,ZPMAT(IP),RPMIN(IP),RPMAX(IP),
c$$$     &            0.03/XRLP(IP),1e+3*XELOSP(IP)/0.03
c$$$ 75          FORMAT('ID = ',I2,' ; Z = ',F5.1,' ; R1 = ',F4.1,
c$$$     &            ' ; R2 = ',F4.1,' ; X0 = ',F5.2,' ; DEDX = ',F5.2)

             IF(IDEB.GE.2) WRITE(IPRNT,
     $            '(''       mat. plane'',I3,''  param.'','//
     $            'T39,3F11.3,2F11.4,F11.6,''  proj. dist.'',F8.2)')
     $            IP,PARINT,ALRPH
*
            NSURF=NSURF+1
            ALRPHI(NSURF)=ALRPH
            TYPSUR(NSURF)='PLAN'
            XRLS(NSURF)=XRLP(IP)/DABS(DCOS(PARINT(4)))
c$$$            WRITE(*,*) 'XINT=',PARINT(2),' YINT=',parint(3),
c$$$     &                 ' ZINT=',PARINT(1)
            CALL FKDPCO(PARINT,PARAMS(1,NSURF),6)
            XELOSS(NSURF)=XELOSP(IP)/ABS(COS(PARINT(4)))
          ENDIF
         ENDIF
        ENDIF
   40 CONTINUE
*
*  sort all intersections together according to the length in r-phi projection
*
      CALL SORTZV(ALRPHI,IORDS,NSURF,1,0,0)
*
      IF(IDEB.GE.2) WRITE(IPRNT,
     & '(''   surfaces ordered by dist. on xy proj.'',20I4,/,40X,20I4)')
     & (INUMB(IORDS(K)),K=1,NSURF)
*
*   store first and last measurement surface (TE) in IFSURF and ILSURF
*   topological check : the types of measurement surfaces should be,
*   either all cylinders, or all planes, or cylinders then planes
*   (never interspersed)
*
      IFSURF=0
      ILSURF=0
      TYPOLD='NONE'
*
*   determine the first and last surface with a measurement
      DO 50 IS=1,NSURF
        INUM=INUMB(IORDS(IS))
        IF(INUM.GT.0) THEN
*
          TYPOLD=TYPE(INUM)
          IF(IFSURF.EQ.0) IFSURF=IS
          ILSURF=IS
        ENDIF
   50 CONTINUE
*
*   abnormal situation : no intersection of ref. traj. with meas. surfaces
      IF(TYPOLD.EQ.'NONE') THEN
        WRITE(IPRNT,
     &  '('' (TKF) ref. traj. does not cross any meas. surface !'')')
        IERB=3
        IERF=3

        RETURN
      ENDIF
*
*   compute reference trajectory at all surfaces, including dEdx effects
*
*   surfaces crossed after the reference point
      CALL UCOPY(PARREF,PARINT,2*6)
      TYPOLD = 'CYLI'
      DO 52 IS=1,NSURF
        ISURF = IORDS(IS)
        IF(ALRPHI(ISURF).GE.0.) THEN
          IF(TYPSUR(ISURF).NE.TYPOLD) THEN
            CALL FKDPV0(DQINT,5)
            CALL FKDPV0(WGINT,15)
            CALL FKSWCH(TYPOLD,PARINT,DQINT,WGINT)
            TYPOLD = TYPSUR(ISURF)
          ENDIF
          if(typsur(isurf).eq.'CYLI') then
            iopt = 0
            idir = 0
            radf = real(params(1,isurf))
            call fkxcyl(parint,idir,radf,
     $           -500.0,500.0,sinbmx,iopt,parint,der,alrph,ier)
          else
            iopt = 0
            idir = 0
            zf   = real(params(1,isurf))
            call fkxpla(parint,idir,zf,0.0,500.0,iopt,
     &           parint,der,alrph,ier)
          endif
          call ucopy(parint,params(1,isurf),2*6)
          call ucopy(parint,paraft(1,isurf),2*6)
          if(inumb(isurf).lt.0) then
            call ucopy(parint,paramx(1,-inumb(isurf)),2*6)
          else if(inumb(isurf).eq.0) then
*
*  ---      dEdx without approximations mif 19/2/96
*
            SIGNP=SIGN(1.D0,PARINT(6))
            STH=SIN(PARINT(4))
            PABS=ABS(CONSB/(PARINT(6)*STH))
            ETOT=DSQRT(AMDEDX**2+PABS**2)
            call fkdedx(real(pabs),real(amdedx),fdedx,dfdedx)
            denerg    = fdedx*xeloss(isurf)
            etot      = etot - denerg
            pabs      = dsqrt(max(etot**2-amdedx**2,0.0001d0))
            parint(6) = consb/(pabs*sth)*signp
            call ucopy(parint,paraft(1,isurf),2*6)
          endif
        ENDIF
 52   CONTINUE
*
*     surfaces crossed before the reference point
*
      call ucopy(parref,parint,2*6)
      typold = 'CYLI'
      do 66 is=nsurf,1,-1
        isurf = iords(is)
        if(alrphi(isurf).le.0.) then
          if(typsur(isurf).ne.typold) then
            CALL FKDPV0(DQINT,5)
            CALL FKDPV0(WGINT,15)
            call fkswch(typold,parint,dqint,wgint)
            typold = typsur(isurf)
          endif
          if(typsur(isurf).eq.'CYLI') then
            iopt = 0
            idir = 0
            radf = real(params(1,isurf))
            call fkxcyl(parint,idir,radf,
     $           -500.0,500.0,sinbmx,iopt,parint,der,alrph,ier)
          else
            iopt = 0
            idir = 0
            zf   = real(params(1,isurf))
            call fkxpla(parint,idir,zf,0.0,500.0,iopt,
     &       parint,der,alrph,ier)
          endif
          call ucopy(parint,params(1,isurf),2*6)
          call ucopy(parint,paraft(1,isurf),2*6)
          if(inumb(isurf).lt.0) then
            call ucopy(parint,paramx(1,-inumb(isurf)),2*6)
          else if (inumb(isurf).eq.0) then
*
* --- dEdx without approximations mif 19/2/96
*
            SIGNP=SIGN(1.D0,PARINT(6))
            STH=SIN(PARINT(4))
            PABS=ABS(CONSB/(PARINT(6)*STH))
            ETOT=DSQRT(AMDEDX**2+PABS**2)
            call fkdedx(real(pabs),real(amdedx),fdedx,dfdedx)
            denerg    = fdedx*xeloss(isurf)
            ETOT = ETOT + DENERG
            PABS=DSQRT(ETOT**2-AMDEDX**2)
            PARINT(6)= CONSB/(PABS*STH)*SIGNP

            call ucopy(parint,params(1,isurf),2*6)
          endif
        endif
 66   continue
*
*   store the array of (measured quantities - extrapolated ones)
*   put zero for unmeasured quantities
*
      do 67 is=1,nsurf
        isurf = iords(is)
        ite   = inumb(isurf)
        if (ite.gt.0) then
          do 68 i=1,5
            if(jbit(measur(ite),i).eq.1) then
              dqter(i,ite) = param(i+1,ite)-params(i+1,isurf)
              if ((type(ite).eq.'CYLI' .and. i.eq.1) .or. i.eq.4)
     &             dqter(i,ite) = dmod(dqter(i,ite)+twopi+pi,twopi)-pi
              if (ideb.ge.1) write(txtdq(i),'(f11.6)') dqter(i,ite)
            else IF(IDEB.GE.1) THEN
              txtdq(I)='           '
            endif
 68       continue
*
          IF(IDEB.GE.1) WRITE(IPRNT,'(''       int. TE'',I3,2X,'//
     $         'A4,''  param.'',T39,3F11.3,2F11.4,'//
     $         'F11.6,''  proj. dist.'','//
     $         'F8.2,/,''       deviat. from ref.'',T50,5A11)')
     $       ITE,TYPE(ITE),(PARAMS(I,ISURF),I=1,6),ALRPHI(ISURF),TXTDQ
*
        end if
 67   continue
*
      CALL VZERO(PROBS,NTEMX)
*
*   backward KALMAN filter -----------------------------------------------------
*
      IPASS = 0
   60 CONTINUE
      IPASS = IPASS + 1
*
      CALL FKDPV0(DQBCK1,5*NTEMX)
      CALL FKDPV0(WGBCK1,15*NTEMX)
      CALL VZERO (NDBCK1,NTEMX)
      CALL VZERO (CHBCK1,NTEMX)
      CALL FKDPV0(DQBCK2,5*NTEMX)
      CALL FKDPV0(WGBCK2,15*NTEMX)
      CALL VZERO (NDBCK2,NTEMX)
      CALL VZERO (CHBCK2,NTEMX)
      CALL VZERO (IRFBCK,NTEMX)
      CALL VZERO (IRZBCK,NTEMX)
      CALL FKDPV0(DQBCKE,5*NEXSMX)
      CALL FKDPV0(WGBCKE,15*NEXSMX)
*
      IF(IDEB.GT.1) WRITE(IPRNT,
     & '(/,7X,''Begin backward Kalman Filter'')')
*
      CALL FKKALM(-1,ILSURF,IDEB,IHIS,
     &            NTEFIT,IPATTB,ISREFB,IERB)
*
*   forward KALMAN filter -----------------------------------------------------
*
      CALL FKDPV0(DQFRW1,5*NTEMX)
      CALL FKDPV0(WGFRW1,15*NTEMX)
      CALL VZERO (NDFRW1,NTEMX)
      CALL VZERO (CHFRW1,NTEMX)
      CALL FKDPV0(DQFRW2,5*NTEMX)
      CALL FKDPV0(WGFRW2,15*NTEMX)
      CALL VZERO (NDFRW2,NTEMX)
      CALL VZERO (CHFRW2,NTEMX)
      CALL VZERO (IRFFRW,NTEMX)
      CALL VZERO (IRZFRW,NTEMX)
      CALL FKDPV0(DQFRWE,5*NEXSMX)
      CALL FKDPV0(WGFRWE,15*NEXSMX)
*
      IF(IDEB.GT.1) WRITE(IPRNT,
     & '(/,7X,''Begin forward Kalman Filter'')')
*
      CALL FKKALM(1,IFSURF,IDEB,IHIS,
     &            NDUMMY,IPATTF,ISREFF,IERF)
*
      IF(IERB.NE.0.OR.IERF.NE.0) RETURN
*
*** New distinction BARREL vs. FORWARD (fwd. defined by hitting FCA) ***
      IF (IPASS.EQ.1)   THEN
        IJWBAR=JBIT(IPATTB,20)
        IF(IJWBAR.EQ.0)THEN
          IBAFO=1
        ELSE
          IBAFO=0
        ENDIF
      ENDIF

*   IBAFO (0 = forward, 1 = barrel) is the only flag used from now on ***
*
**********************************************************************
*   Following code removes internal TE chisquares and
*   TE d.o.f. : mif: do for search, but not for fit.

       IF (ICHI.EQ.1) THEN

         CHITE=0.
         NDGTE=0
         IF(NTEFIT.GT.1) THEN
           DO 63 ITE=1,NTE
             IF(IERTE(ITE).EQ.0.AND.IFLAG(ITE).EQ.0) THEN
               NDGTE=NDGTE+NDTER(ITE)
               CHITE=CHITE+CHTER(ITE)
             ENDIF
   63      CONTINUE
           NDBACK=NDBACK-NDGTE
           CHBACK=CHBACK-CHITE
           NDFORW=NDFORW-NDGTE
           CHFORW=CHFORW-CHITE
         END IF
       END IF
*
*   Skip outlier logic if switched off R.E./R.F.
*
      if (iout.eq.0) return
      IF (IDEB.GE.1) WRITE(IPRNT,'(//1X,3A,I4)')
     &   'FKTFIT: ', FBTEXT(IBAFO),
     &   ' track passed bw.&fw. Kalman filters - PASS =', IPASS
*
*   Skip outlier logic if track fit's total probability is acceptable
*   (or if n.d.f. is equal to zero)
*
      IF (NDBACK.GT.0) THEN
        PROBX = PROB(CHBACK,NDBACK)
c MVL        PROBX = FXPROB(CHBACK,NDBACK)

c MVL      IF(CHBACK.GT.0.0.and.CHBACK.lT.300.0)THEN
c MVL        PROBX=max(PROB(CHBACK,NDBACK),1.e-25)
c MVL      ELSE
c MVL        write(*,*) ' PROBX --- reason2',CHBACK,NDBACK
c MVL        PROBX=1.e-25
c MVL      ENDIF
c MVL

        IF (PROBX.GE.PRTOT(IBAFO)) THEN
          IF (IDEB.GE.1) WRITE(IPRNT,'(/9X,A,I8,G13.6)')
     &       'Good track fit - NDF, PROB =', NDBACK, PROBX
          GOTO 99
        ENDIF
      ELSE
        IF (IDEB.GE.1) WRITE(IPRNT,'(/9X,A,I8,A)')
     &     'No defined fit - NDF =', NDBACK, ' - return'
        GOTO 99
      ENDIF
      if (itan.ne.0) then
        ideb_save = ideb
      end if
      IF (IDEB.GE.1) then
        WRITE(IPRNT,'(/9X,A,I8,G13.6/9X,A)')
     &       'Bad track fit  - NDF, PROB =', NDBACK, PROBX,
     &       'Entering the outlier logic for TEs:'
        write(IPRNT,'(12x,10Z8)') (idter(i),i=1,nte)
      end if
*
**********************************************************************
*   Enter outlier logic if this track would be fired by the total chisq
*   probability cut: either "special" or "standard" logic
**********************************************************************
*
*   If less than 3 TEs: special outlier logic (favouring the TPC)
*
      IF     (NTEFIT.EQ.1) THEN
*       if there is only one TE, it is a TPC, thus:
        IF (IDEB.GE.1) WRITE(IPRNT,'(/9X,A)')
     &     'Only one TE (it is a TPC) - keep it & return'
        GOTO 99
      ELSEIF (NTEFIT.EQ.2) THEN
*       look for a possibly existing TE of the TPC:
        JTPCTE=0
        JOTHER=0
        DO 65 ITE=1,NTE
          IF(IERTE(ITE).EQ.0.AND.IFLAG(ITE).EQ.0) THEN
            IF(IDET(ITE).EQ.detid(i_TPC)) THEN
              JTPCTE=ITE
            ELSE
              JOTHER=ITE
            ENDIF
          ENDIF
  65    CONTINUE
        IF(JTPCTE.NE.0) THEN
*         there is a TPC among the two TEs -> mark the other as outlier,
*         submit once more for the track fit with the other TE removed:
          IFLAG(JOTHER) = 1
          IF (IDEB.GE.1) WRITE(IPRNT,'(/9X,A,I6,A)')
     &       'Two TEs (IDET = detid(i_TPC),', IDET(JOTHER),
     &       ') - keep TPC & refit'
          GOTO 60
        ELSE
*--------------------- choose either "strong" logic ------------------
*         there is no TPC present -> bad track, reject it completely:
          IERB = 4
          IERF = 4
          IF (IDEB.GE.1) WRITE(IPRNT,'(/9X,A)')
     &       'Two TEs but without a TPC - reject this track'
*--------------------- or "loose" logic ------------------------------
*         there is no TPC present -> bad track, but keep it for now:
*          IF (IDEB.GE.1) WRITE(IPRNT,'(/9X,A)')
*     &       'Two TEs but without a TPC - keep track & return'
*---------------------------------------------------------------------
          GOTO 99
        ENDIF
      ENDIF
*
*
*   If at least 3 TEs: standard outlier logic (combine forward and backward
*   filters' predictions to calculate total chisquares with one TE removed)
*
      IF (IDEB.GE.1) THEN
        WRITE(IPRNT,'(/9X,A,I4,A,i2,a)')
     $       'No. of the TEs =', NTEFIT,
     &       ' >= 3 - standard outlier logic (iout = ',iout,' ) :'
        WRITE(IPRNT,'(/9X,A/9X,A/9X,A/)')
     &       '  removed TE    interpolation: fit without TE    outlier',
     &       '--------------  -----------------------------    candi-',
     &       'INUM IDET RANK  NDF        CHISQ         PROB    date ?'
      ENDIF
*
*   Loop over all active TEs
*
      INUMX =  0
      NFLAG =  0
      NDMAX = -1
      IODIN=0
      IFWDIN=0
      DO 70 INUM = 1,NTE
        pchint(inum) = -1.
        IF(IFLAG(INUM).EQ.0) THEN
*
*   if rank of information without this TE is incomplete: enforce positive
*   definitness by slightly increasing the weight matrices ("Billoir's trick")
*
          IRNKRF=MIN0(3,IRFFRW(INUM)+IRFBCK(INUM))
          IRNKZT=MIN0(2,IRZFRW(INUM)+IRZBCK(INUM))
          IF(IRNKRF.LT.1) THEN
            WGBCK1(1,INUM)=WGBCK1(1,INUM)+1.
            WGFRW1(1,INUM)=WGFRW1(1,INUM)+1.
          ENDIF
          IF(IRNKRF.LT.2) THEN
            WGBCK1(10,INUM)=WGBCK1(10,INUM)+1.
            WGFRW1(10,INUM)=WGFRW1(10,INUM)+1.
          ENDIF
          IF(IRNKRF.LT.3) THEN
            WGBCK1(15,INUM)=WGBCK1(15,INUM)+1.
            WGFRW1(15,INUM)=WGFRW1(15,INUM)+1.
          ENDIF
*
          IF(IRNKZT.LT.1) THEN
            WGBCK1(3,INUM)=WGBCK1(3,INUM)+1.
            WGFRW1(3,INUM)=WGFRW1(3,INUM)+1.
          ENDIF
          IF(IRNKZT.LT.2) THEN
            WGBCK1(6,INUM)=WGBCK1(6,INUM)+1.
            WGFRW1(6,INUM)=WGFRW1(6,INUM)+1.
          ENDIF
*
*   compute interpolation (weighted mean of filters without this TE)
*
          CALL FKWGMN(IRNKRF+IRNKZT,
     &      DQBCK1(1,INUM),WGBCK1(1,INUM),NDBCK1(INUM),CHBCK1(INUM),
     &      IRFBCK(INUM)+IRZBCK(INUM),
     &      DQFRW1(1,INUM),WGFRW1(1,INUM),NDFRW1(INUM),CHFRW1(INUM),
     &      IRFFRW(INUM)+IRZFRW(INUM),
     &      DQINT,WGINT,NDINT,CHINT,IERR)
          IF(IERR.NE.0) THEN
            WRITE(IPRNT,*) ' FKTFIT(TKF) : error in FKWGMN'
            IF(IDEB.GE.1) WRITE(IPRNT,'(7X,''FKWGMN failed '')')
*            PROBS(INUM)=-1.
            GO TO 70
          ENDIF
**********************************************************************
*   subtract internal TE chisquares and TE n.d.f. (SD)
*
          IF (ICHI.EQ.1) THEN
            NDINT=NDINT-NDGTE+NDTER(INUM)
            CHINT=CHINT-CHITE+CHTER(INUM)
          END IF
**********************************************************************
          if (ndint.gt.0)  then
            pchint(inum) = prob(chint,ndint)
c MVL           IF(CHINT.le.0.0) write(*,*) ' Pchint',CHINT,ndint
c MVL            pchint(inum) = fxprob(chint,ndint)

c MVL      IF(CHINT.GT.0.0.and.CHINT.lT.300.0)THEN
c MVL        pchint(inum) = max(PROB(CHINT,NDINT),1.e-25)
c MVL      ELSE
c MVL        write(*,*) ' Pchint --- reason3',CHINT
c MVL        pchint(inum) = 1.e-25
c MVL      ENDIF
c MVL
c            print*,chint,ndint,pchint(inum) ,CHBCK1(INUM),CHFRW1(INUM)
          else
            pchint(inum) = 0.
          endif
*
*   compute smoothed chisquare (difference total - interpolation) and
*   filtered chisquare - used for histograms only
*
          DCHI2=CHBACK-CHINT
          NDDEG=NDBACK-NDINT
          IF (NDDEG.GT.0.AND.DCHI2.GT.0.) THEN
            PROBS(INUM)=PROB(DCHI2,NDDEG)
          ELSE
            PROBS(INUM)=1.
          ENDIF
          DCHI2F=CHBCK2(INUM)-CHBCK1(INUM)-CHTER(INUM)
          ndint0 = ndint
*
* mif:   redefine NDINT as a score relative to the importance of the detector:
*
          IF (IDET(INUM).EQ.detid(i_TPC)) THEN
C -----------------------------TPC
                                      NDINT=1
          ELSE IF (IDET(INUM).EQ.detid(i_VTX)) THEN
C -----------------------------VTX
                                      NDINT=2
          ELSE IF (IDET(INUM).EQ.detid(i_ITC)) THEN
C -----------------------------ITC
                                      NDINT=3
          ELSE IF (IDET(INUM).EQ.detid(i_FTD)) THEN
C -----------------------------FTD
                                      NDINT=6
                                      IFWDIN=1
          END IF

*
*    mif --check that rank of matrix is not decreased by throwing out TE for iout=6
*
        IF (IOUT.EQ.6) THEN
           IF (IRNKRF.LT.IRNKXY .OR. IRNKZT.LT.IRNKZ) THEN
             if (ideb.ge.1) write(IPRNT,'(9x,i4,2i3,i8,2e13.5,'//
     &          '4x,''no-forbidden by IOUT='',i2)')
     &          inum,idet(inum),ndint,ndint0,chint,pchint(inum),iout
             goto 70
           END IF
        END IF
*
*    mif --check that ITC is not thrown out for iout=1
*
        IF (IOUT.EQ.1) THEN
           IF (IDET(INUM).EQ.detid(i_ITC)) THEN
             if (ideb.ge.1) write(IPRNT,'(9x,i4,i6,i8,2e13.5,'//
     &          '4x,''no-forbidden by IOUT='',i2)')
     &          inum, idet(inum), ndint, chint, pchint(inum), iout
             GOTO 70
           END IF
        END IF
*
*    mif --check that TPC is not thrown out for iout=2
*
        IF (IOUT.EQ.2) THEN
           IF (IDET(INUM).EQ.detid(i_TPC)) THEN
             if (ideb.ge.1) write(IPRNT,'(9x,i4,i5,i9,2e13.5,'//
     &          '4x,''no-forbidden by IOUT='',i2)')
     &          inum, idet(inum), ndint, chint, pchint(inum), iout
             GOTO 70
           END IF
        END IF
*
* ME --check that FTD is not thrown out for iout=8
*
         IF (IOUT.EQ.8) THEN
           IF (IDET(INUM).EQ.detid(i_FTD)) THEN
             if (ideb.ge.1) write(IPRNT,'(9x,i4,i5,i9,2e13.5,'//
     &          '4x,''no-forbidden by IOUT='',i2)')
     &          inum, idet(inum), ndint, chint, pchint(inum), iout
             GOTO 70
           END IF
         END IF
*
*   only if acceptable total probability after interpolation: store
*   serial number of TE with highest n.d.f. (in case of equal n.d.f.s
*   take the one with higher probability)
*
          if (pchint(inum).ge.prtot(ibafo)) then
            nflag = nflag + 1
            if     (ndint.gt.ndmax) then
              inumx = inum
              ndmax = ndint
              prmax = pchint(inum)
            elseif (ndint.eq.ndmax) then
              if (pchint(inum).gt.prmax) then
                inumx = inum
                prmax = pchint(inum)
              endif
            endif
            fltext = '- yes '
          else
            fltext = '- no  '
          endif
          if (ideb.ge.1) write(IPRNT,'(9x,i4,2i5,i4,2e13.5,4x,a)')
     $         inum, idet(inum), ndint, ndint0, chint, pchint(inum),
     $         fltext
        ENDIF
*
   70 CONTINUE
*
*   test whether the standard outlier logic has been succesful - if yes,
*   flag the "worst TE" and refit the track without it
*
*     logic:  if NO TE is flaged (NFLAG=0)
*             and fk3trk is NOT called from track search (itan>0)
*             remove TE with highest probability and fit
*             again (maximum 2 times --> ipass < 3)
*             --> one can remove max. 2 TEs from track befor rejecting
*
      i_reject = 0
      if (nflag.eq.0) then
        if (itan.eq.0) then
          i_reject = 1
        else
CC          if (ipass.ge.3) then
          if (ipass.ge.30) then
            i_reject = 1
          end if
        end if
      end if

      if (i_reject.eq.1) then
*
*                      reject track
*--------------------- choose either "strong" logic ------------------
*       bad interpol. prob. for all TEs after substraction of 2 TEs
*       --> bad track, reject it completely:
*
        ierb = 4
        ierf = 4
        if (ideb.ge.1) write(IPRNT,'(/9x,a)')
     &     'Outlier test - no success - reject this track'
*--------------------- or "loose" logic ------------------------------
*       bad interpol. prob. for all TEs -> bad track, but keep it for now:
*        IF (IDEB.GE.1) WRITE(IPRNT,'(/9X,A)')
*---------------------------------------------------------------------
*     &     'Outlier test - no success - keep track & return'
        goto 99
      else
        if (nflag.eq.0.and.itan.gt.0) then
*
*     bad interpol. prob. for all TEs -> bad track, find TE with
*     worst prob., remove it and fit again, but don't reject
*     TPC TE, if it is not the outerst detector
*
          call sortzv(pchint,ipchint,nte,1,1,0)
          inumx = ipchint(1)
          !
          ! Don't reject TPC, if there is an outer detector (OD,FCA,FCB)
          !
          call sortzv(idet,iidet,nte,0,1,0)
          if (idet(inumx).eq.detid(i_TPC).and.
     $        idet(iidet(1)).gt.idet(inumx)) then
            inumx = ipchint(2)
            if (ideb.ge.1) write(IPRNT,'(/9x,a)')
     $           'TPC is not outerst detector, do not reject TE'
          end if
        end if
*
*     outlier test was succesful
*
*
*      mark the "worst TE" as an outlier,
*      submit once more for the track fit with that TE removed:
*
        iflag(inumx) = 1
        if (ideb.ge.1) write(IPRNT,'(/9x,a,2i6,a)')
     &     'removal of INUM/IDET =', INUMX, IDET(INUMX),
     &     ' - keep other TEs & refit'
*
*     reject track, if no ITC rejection (IOUT=1) but worst TE is ITC
*
        if (iout.eq.1 .and. idet(inumx).eq.detid(i_ITC)) then
          ierb = 4
          ierf = 4
          if (ideb.ge.1) write(IPRNT,'(/9x,a)')
     $         'no ITC rejection (IOUT=1) but worst TE is ITC --> '//
     $         'track rejected'
          return
        end if
*
*     reject track, if no TPC rejection (IOUT=2) but worst TE is TPC
*
        if (iout.eq.2 .and. idet(inumx).eq.detid(i_TPC)) then
          ierb = 4
          ierf = 4
          if (ideb.ge.1) write(IPRNT,'(/9x,a)')
     $         'no TPC rejection (IOUT=2) but worst TE is TPC --> '//
     $         'track rejected'
          return
        end if
*
* -ME-reject track, if no FTD rejection (IOUT=8) but worst TE is FTD
*
        if (iout.eq.8 .and. idet(inumx).eq.detid(i_FTD)) then
          ierb = 4
          ierf = 4
          if (ideb.ge.1) write(IPRNT,'(/9x,a)')
     $         'no FTD rejection (IOUT=8) but worst TE is FTD --> '//
     $         'track rejected'
          return
        end if
*
        goto 60
      endif
*
   99 CONTINUE
      if (itan.ne.0) ideb = ideb_save
      RETURN
      END

