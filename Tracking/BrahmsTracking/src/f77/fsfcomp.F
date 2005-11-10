*  
* $Id: fsfcomp.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfcomp(paramf,te1,det1,te2,det2,
     1                   dcx,dcy,dcRPhiR,dcz,
     2                   dctheta,dcphi,dcthest,dcphist,
     3                   goodcomp)
*******************************************************************************
*
*     check if final helix parameters coincide with dat2
*     in detail: check each measured quantity (except z/R)
*                against the extrapolations
*                prediction and declare it as a goodcomp if all checks are
*                successful
*
*     INPUT:  paramf  - extrapolated parameters
*             te1     - index of TE1 to extrapolate to TE2
*             det1    - detector index
*             te2     - index of second TE
*             det2    - detector index
*             dcx,dcy,dcRPhiR,dcz,
*             dctheta,dcphi,dcthest,dcphist - cut values !
*     OUTPUT: goodcomp - .TRUE. for good pair
*                                                                 Markus Elsing
*******************************************************************************
      implicit none
#include "include/fsfb.inc"
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfpi.inc"
      logical goodcomp
      integer te1,det1,te2,det2,datdescr1,datdescr2,jbit
      real    dcx(detmax,detmax),dcy(detmax,detmax)
      real    dcRPhiR(detmax,detmax),dcz(detmax,detmax)
      real    dctheta(detmax,detmax),dcphi(detmax,detmax)
      real    dcthest(detmax,detmax),dcphist(detmax,detmax)
      real    paramf(6),deltaphi
      real    delta,fsfpdiff,thest,thestp,phist,phistp,phip

      goodcomp = .FALSE.
      datdescr1=datdescr(te1,det1)
      datdescr2=datdescr(te2,det2)
c
c --- test on RPhi/R and z ----------------------------------------------------
c
      if(1.eq.jbit(datdescr2,1))then   ! rphi version !
c --- transform to RPhi/R
        phip=atan2(paramf(3),paramf(2))        ! Phi
        if (phip.lt.0.) phip = phip + 2.*pi
c --- RPhi/R
        delta=abs(fsfpdiff(phip,RPhi(te2,det2)/R(te2,det2)))
        if(delta.ge.dcRPhiR(det1,det2).and.
     1     dRPhi(te2,det2).lt.1e5)then
          if(fsfdbg.gt.1.or.(labelsok.and.fsfdbg.gt.0))then
            print*,'FSFCOMP: RPhi/R is off:',phip,'-',
     1             RPHi(te2,det2)/R(te2,det2),' sigma:',
     1             dcRPhiR(det1,det2)
          endif
          return
        endif
c --- z
        if(abs(paramf(1)-z(te2,det2)).ge.dcz(det1,det2).and.
     1     dz(te2,det2).lt.1e5) then
          if(fsfdbg.gt.1.or.(labelsok.and.fsfdbg.gt.0))then
            print*,'FSFCOMP: z is off: ',paramf(1),'-',
     1             z(te2,det2),' sigma:',dcz(det1,det2)
          endif
          return
        endif
c
c --- test on x and y ---------------------------------------------------------
c
      else                             ! xyz version
c --- x
        if(abs(paramf(2)-x(te2,det2)).ge.dcx(det1,det2).and.
     1     dx(te2,det2).lt.1e5) then
          if(fsfdbg.gt.1.or.(labelsok.and.fsfdbg.gt.0))then
            print*,'FSFCOMP: x is off:',paramf(2),'-',x(te2,det2),
     1             ' sigma:',dcx(det1,det2)
          endif
          return
        endif
c --- y
        if(abs(paramf(3)-y(te2,det2)).ge.dcy(det1,det2).and.
     1     dy(te2,det2).lt.1e5) then
           if(fsfdbg.gt.1.or.(labelsok.and.fsfdbg.gt.0))then
             print*,'FSFCOMP: y is off:',paramf(3),'-',y(te2,det2),
     1              ' sigma:',dcy(det1,det2)
           endif
           return
        endif
c --- theta*
        if (dy(te2,det2).lt.1e5.and.dx(te2,det2).lt.1e5.and.
     1       z(te2,det2).ne.0. .and.paramf(1).ne.0.) then
          thest  = ATAN(sqrt(x(te2,det2)**2+y(te2,det2)**2)
     1             /abs(z(te2,det2)))
          thestp = ATAN(sqrt(paramf(2)**2+paramf(3)**2)
     1             /abs(paramf(1)))
          thest  = min(thest ,pi-thest )
          thestp = min(thestp,pi-thestp)
          if(abs(thestp-thest).ge.dcthest(det1,det2)) then
            if(fsfdbg.gt.1.or.(labelsok.and.fsfdbg.gt.0))then
              print*,'FSFCOMP: theta* is off:',thestp,'-',thest,
     1               ' sigma:',dcthest(det1,det2)
            endif
            return
          endif
c --- phi*
          phist  = ATAN2(y(te2,det2),x(te2,det2))
          phistp = ATAN2(paramf(3),paramf(2))
          if (phist .lt.0.) phist  = phist  + 2.*pi
          if (phistp.lt.0.) phistp = phistp + 2.*pi
          deltaphi = abs(fsfpdiff(phist,phistp))
          if(deltaphi.ge.dcphist(det1,det2)) then
            if(fsfdbg.gt.1.or.(labelsok.and.fsfdbg.gt.0))then
              print*,'FSFCOMP: phi* is off:',phistp,'-',phist,
     1               ' sigma:',dcphist(det1,det2)
              print*,'         deltaphi',deltaphi,'=',
     1              phistp,'-',phist
            endif
            return
          endif
        endif
      endif
c
c --- test on theta, phi -------------------------------------------------
c
c --- theta
      if(abs(paramf(4)-theta(te2,det2)).ge.dctheta(det1,det2).and.
     1   dtheta(te2,det2).lt.1e5) then
         if(fsfdbg.gt.1.or.(labelsok.and.fsfdbg.gt.0))then
           print*,'FSFCOMP: theta is off:',paramf(4),'-',
     1            theta(te2,det2),' sigma:',dctheta(det1,det2)
           if(abs(paramf(4)-theta(te2,det2)).ge.0.5)then
              print*,'WARNING in FSFCOMP: dtheta > 0.5 !!!'
           endif
         endif
         return
      endif
c --- phi
      deltaphi=abs(fsfpdiff(paramf(5),phi(te2,det2)))
      if(deltaphi.ge.dcphi(det1,det2).and.
     1   dphi(te2,det2).lt.1e5) then
         if(fsfdbg.gt.1.or.(labelsok.and.fsfdbg.gt.0))then
           print*,'FSFCOMP: phi is off:',paramf(5),'-',
     1            phi(te2,det2),' sigma:',dcphi(det1,det2)
           print*,'         deltaphi',deltaphi,'=',
     1            paramf(5),'-',phi(te2,det2)
         endif
         return
      endif
c
c --- ok, take it !
c
      goodcomp=.TRUE.
      end


