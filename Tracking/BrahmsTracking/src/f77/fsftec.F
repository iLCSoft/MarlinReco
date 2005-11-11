*  
* $Id: fsftec.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsftec(firstdet,pconst,method)
*******************************************************************************
*
*     Copy TE information to local arrays and do some
*     preprocessing for the extrapolation starting from the 'firstdet'.
*
*      INPUT:  firstdet - IDcode of detector to use as pivot
*              pconst   - if TRUE redefine TE using p from beam spot constraint
*              method   - work on detector, if method(det)<>0
*                                                                  Daniel Wicke
*******************************************************************************
      implicit none
#include "include/brcntl.inc"
#include "include/fsfb.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfxerrs.inc"
#include "include/fsfpi.inc"
      integer jbit,FSFSHWR,method(detmax)
      integer idet,ite,ndat,firstdet,i,errindex(5),index,ierr
      logical pconst
c -----------------------------------------------------------------------------
c     errindex shall hold the position of the error of the different
c     coordinates:
c
c     errindex(1) -> t1
c     errindex(2) -> t2
c     errindex(3) -> theta
c     errindex(4) -> phi
c     errindex(5) -> 1/p or 1/p_t or E
c -----------------------------------------------------------------------------

      do idet=1,detmax
       if (method(idet).ne.0) then
        do ite=1,teindexused(idet)

           datdescr(ite,idet)=ispace(teindex(ite,idet)+4)

* --- get if the TE is after a secondary interaction
           showerte(ite,idet)=FSFSHWR(teid(ite,idet))

           ndat=tesize(ite,idet)

           if(0.eq.jbit(datdescr(ite,idet),1))then !x y z
             x(ite,idet)=space(teindex(ite,idet)+10)
             y(ite,idet)=space(teindex(ite,idet)+11)
             R(ite,idet)   =sqrt(x(ite,idet)**2+y(ite,idet)**2)
             Rphi(ite,idet)=atan2(y(ite,idet),x(ite,idet))
c --- catch phi problems in TEs - rest of the code would crash !
             if(Rphi(ite,idet).lt.0)
     1          Rphi(ite,idet)=Rphi(ite,idet)+2.*pi
             if(Rphi(ite,idet).gt.2*pi)
     1          Rphi(ite,idet)=Rphi(ite,idet)-2.*pi
             Rphi(ite,idet)=R(ite,idet)*Rphi(ite,idet)
           else                                    ! Rphi
             R(ite,idet)    = space(teindex(ite,idet)+10)
             Rphi(ite,idet) = space(teindex(ite,idet)+11)
             x(ite,idet)=R(ite,idet)*cos(RPHI(ite,idet)/R(ite,idet))
             y(ite,idet)=R(ite,idet)*sin(RPHI(ite,idet)/R(ite,idet))
           endif
           z(ite,idet)    =space(teindex(ite,idet)+12)
           theta(ite,idet)=space(teindex(ite,idet)+13)
           phi(ite,idet)  =space(teindex(ite,idet)+14)
c --- catch phi problems in TEs - rest of the code would crash !
           if (phi(ite,idet).lt.0.)   phi(ite,idet)=phi(ite,idet)+2*pi
           if (phi(ite,idet).gt.2*pi) phi(ite,idet)=phi(ite,idet)-2*pi

           if(1.eq.jbit(datdescr(ite,idet),10))then  !dat1(15) is 1/p
             if(1.eq.jbit(datdescr(ite,idet),8) !Check for theta to be valid
     1          .and.theta(ite,idet).ne.0)then
                overrt(ite,idet)=space(teindex(ite,idet)+15)
     1                       *bfieldz/abs(sin(theta(ite,idet)))  !1/R_T
             else
                overrt(ite,idet)=0
             endif
           else if(1.eq.jbit(datdescr(ite,idet),11))then  !dat1(15) is 1/p_T
             overrt(ite,idet)=space(teindex(ite,idet)+15)
     1                       *bfieldz !1/R_T
           else
             overrt(ite,idet)=0
           endif

           index=1
           do i=1,4
             if(1.eq.jbit(datdescr(ite,idet),5+i))then
                errindex(i)=16+(index*(index+1))/2
                index=index+1
             else
                errindex(i)=-1000000 !This should crash it, if used.
             endif
           enddo
           errindex(5)=16+(index*(index+1))/2 !Fill this one always

           if(0.eq.jbit(datdescr(ite,idet),1))then !x y z
              if(jbit(datdescr(ite,idet),6).eq.1)then
                 dx(ite,idet)=
     1             sqrt(abs(space(teindex(ite,idet)+errindex(1))))
              else
                 dx(ite,idet)=1e10
              endif
              if(jbit(datdescr(ite,idet),7).eq.1)then
                 dy(ite,idet)=
     1             sqrt(abs(space(teindex(ite,idet)+errindex(2))))
              else
                 dy(ite,idet)=1e10
              endif
              drphi(ite,idet)=1e10
              dz   (ite,idet)=1e10
           else                                    ! R RPHi z
              if(jbit(datdescr(ite,idet),6).eq.1)then
                 dRPhi(ite,idet)=
     1             sqrt(abs(space(teindex(ite,idet)+errindex(1))))
              else
                 dRPhi(ite,idet)=1e10
              endif
              if(jbit(datdescr(ite,idet),7).eq.1)then
                 dz(ite,idet)=
     1             sqrt(abs(space(teindex(ite,idet)+errindex(2))))
              else
                 dz(ite,idet)=1e10
              endif
              dx(ite,idet)=1e10
              dy(ite,idet)=1e10
           endif
           if(jbit(datdescr(ite,idet),8).eq.1)then
             dtheta(ite,idet)=
     1         sqrt(abs(space(teindex(ite,idet)+errindex(3))))
           else
             dtheta(ite,idet)=1e10
           endif
           if(jbit(datdescr(ite,idet),9).eq.1)then
             dphi(ite,idet)=
     1         sqrt(abs(space(teindex(ite,idet)+errindex(4))))
           else
             dphi(ite,idet)=1e10
           endif
           if(jbit(datdescr(ite,idet),10).eq.1)then
              doverrt(ite,idet)=bfieldz/abs(sin(Theta(ite,idet)))
     1            *sqrt(abs(space(teindex(ite,idet)+errindex(5))))
           else if(jbit(datdescr(ite,idet),11).eq.1)then
              doverrt(ite,idet)=bfieldz
     1            *sqrt(abs(space(teindex(ite,idet)+errindex(5))))
           else
              doverrt(ite,idet)=1e10
           endif
        enddo
       endif
      enddo

      if(pconst)then
        idet=firstdet

c ---   Calculate new track radius for primary detector TEs, assuming
c       that the track starts at the beam spot
        do ite=1,teindexused(idet)

c ---     do beam constraint in separate routine, uses now full trackfit
          call fsfpconst(ite,idet,ierr)

          if(ierr.eq.0)then
            if(fsfdbg.gt.1)then
              print*,'FSFTEC: momentum beamspot constraint accepted.'
            endif
            doverrt(ite,idet)=0
            call sbit1(datdescr(ite,idet),10) !We have some radius.
            call sbit0(datdescr(ite,idet),11) !We have no 1/Pt.
          else
            if(fsfdbg.gt.1)then
              print*,'FSFTEC: momentum beamspot constraint rejected.'
            endif
            overrt(ite,idet)=0
            doverrt(ite,idet)=1e10
            call sbit0(datdescr(ite,idet),10) !We have no radius.
            call sbit0(datdescr(ite,idet),11) !We have no 1/Pt.
          endif
        enddo
      endif
      end

