*  
* $Id: fsfpconst.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfpconst(itepd,pd,ierr)
**********************************************************************
*
* Momentum beam constraint for string search.
*
* INPUT:  itepd    - index in teid array
*         pd       - primary detector
* OUTPUT: ierr  <> 0 if beam constraint rejected
*                                                        Markus ELsing
**********************************************************************
      implicit none
#include "include/exmach.inc"
#include "include/brcntl.inc"
#include "include/fsfb.inc"
#include "include/fsfpi.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfflags.inc"
#include "include/fsfdebug.inc"
      real    xx,yy
      real    overrt_old,overrt_new,difphi,oldtheta,sdifphi
      integer ierr,pd,itepd,jbit

c --- we need some phi here
      if(0.eq.jbit(datdescr(itepd,pd),9))then !phi
        if(fsfdbg.gt.1)then
          print*,'FSFPCONST: TE without phi rejected !'
        endif
        ierr = 1
        return
      else
        ierr = 0
      endif

c --- save curvature and its error
      overrt_old=overrt(itepd,pd)

c --- correct curvature to fit beam spot
      xx=x(itepd,pd)
      yy=y(itepd,pd)
      overrt_new=-(2* ( (yy-bspoty)*cos(phi(itepd,pd))
     2                 -(xx-bspotx)*sin(phi(itepd,pd))) )
     3                /((xx-bspotx)**2+(yy-bspoty)**2)

c --- some debug output
      if (fsfdbg.gt.0)then
        print *,' FSFPCONST: old 1/p : ',overrt_old,
     1                       ' label : ',telabl(itepd,pd)
        print *,'            new 1/p : ',overrt_new
        print *,'            sim 1/p : ',simoverrt(itepd,pd)
      endif

c --- Use this contraint only if the change is below 5sigma
      if(doverrt(itepd,pd).ne.0)then
        if(abs(overrt_new-overrt_old)/doverrt(itepd,pd).gt.5)then
          if(fsfdbg.gt.1)then
            print *,'FSFPCONST: TE rejected by 5 sigma cut !'
          endif
          ierr = 6
          return
        endif
      endif

      overrt(itepd,pd) = overrt_new

* --- ME return old theta if overrt to small
      if(abs(overrt_new).lt.1.0e-06)return

c --- now get theta from z(beamspot) to z(TE)
      oldtheta = theta(itepd,pd)
      sdifphi = overrt(itepd,pd)*((bspoty-yy)*sin(phi(itepd,pd))
     1                           +(bspotx-xx)*cos(phi(itepd,pd)))
      if (sdifphi.lt.1.and.sdifphi.gt.-1) then
        difphi=asin(sdifphi)
      else
        return            ! keep old number for the time
      endif
      theta(itepd,pd)=atan(difphi/overrt(itepd,pd)/(bspotz-z(itepd,pd)))
      if (theta(itepd,pd).lt.0) theta(itepd,pd) = pi+theta(itepd,pd)
c --- catch flips of theta, why ?
      if ((z(itepd,pd).lt.0.and.theta(itepd,pd).lt.pi/2.).or.
     1    (z(itepd,pd).gt.0.and.theta(itepd,pd).gt.pi/2.)) then
        if (fsfdbg.gt.0) print *,' FSFPCONST: flip theta !'
        theta(itepd,pd) = pi-theta(itepd,pd)
      endif
c --- some debug output
      if (fsfdbg.gt.0)then
        print *,' FSFPCONST: old theta : ',oldtheta,
     1                         ' label : ',telabl(itepd,pd)
        print *,'            new theta : ',theta(itepd,pd)
        print *,'            sim theta : ',simtheta(itepd,pd)
        if (abs(theta(itepd,pd)-simtheta(itepd,pd)).gt.pi/2.
     1                     .and.simtheta(itepd,pd).ne.0.) then
          print *,' WARNING in FSFPCONST: crazy theta ?'
        endif
      endif

      end


