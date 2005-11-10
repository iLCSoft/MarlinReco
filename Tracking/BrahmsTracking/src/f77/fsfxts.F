*  
* $Id: fsfxts.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfxts(parami,paramf,ite,idet,ierr)
**********************************************************************
*     Routine to extrapolate PARAMI to detector
*     the TE ITE in detector IDET, using FKXPLA or FKXCYL.
*
*     parami and paramf convention are those of FKXPLA
*     INPUT
*       parami - (z,x,y,theta,phi,1/R)
*       ite    - TE index in the FSFARRAYs.
*       idet   - detector index in the FSFARRAYs.
*     OUTPUT
*       paramf - (z,x,y,theta,phi,1/R)
*       ierr   - .ne.0 if error occured.
*                In this case paramf will be set equal to parami.
*
*                                                         Daniel Wicke
**********************************************************************
      implicit none
c +cde, EXXPRM. ! - no implicit none
#include "include/exxgeo.inc"
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfpi.inc"
      real parami(6),paramf(6)
      integer ite,idet,ierr,i
c      real rphiz(6),rphizf(6),alrphi,tmp,sqrt,atan2
      real alrphi,sqrt,atan2
      double precision rphiz(6),rphizf(6),tmp
      integer exxerr,jbit

      call ucopy (parami,paramf,6)

c --- set common vars from EXXGEO for helix-extrapolation
      thxmin=-100000.0
      thxmax= 100000.0
      phymin=-100000.0
      phymax= 100000.0

      ierr=0
      if(parami(6).eq.0)then
        print*,'ERROR in FSFXTS: no radius given. Skipped extrapolation'
        call ucopy(parami,paramf,6)
        ierr=2
c --- catch theta or phi out of range
      else if (parami(4).lt.0.or.parami(4).gt.pi.or.
     1         parami(5).lt.0.or.parami(5).gt.2*pi) then
        print*,'WARNING in FSFXTS: theta/phi off at input'
        print*,'parameters :',(parami(i),i=1,6)
        ierr = 2
      else if(jbit(ispace(teindex(ite,idet)+4),1).eq.1)then !RPhiz
        rphiz(1)=sqrt(parami(2)**2+parami(3)**2)
        if(rphiz(1).eq.0.) then
          rphiz(2)=parami(5)  ! use Phi=phi
          rphiz(5)=0.         ! beta = 0.
        else
          rphiz(2)=atan2(parami(3),parami(2))
          if (rphiz(2).lt.0) rphiz(2)=rphiz(2)+2.*pi
          call fsfbeta (dble(parami(5)),rphiz(2),rphiz(5))
        endif
        rphiz(3)=parami(1)
        rphiz(4)=parami(4)
        rphiz(6)=parami(6)

        if(rphiz(4).eq.0.or.rphiz(6).eq.0)then
           exxerr=-1
        else
           call fkxcyl(rphiz,0,space(teindex(ite,idet)+10),
     1                 -1e10,1e10,1.0,0,rphizf,tmp,alrphi,exxerr)
        endif
        if(exxerr.ne.0)then
           if (fsfdbg.gt.1)then
             print*,'FSFXTS: FKXCYL failed. ierr =',exxerr
           endif
           call ucopy(parami,paramf,6)
           ierr=1
        else
           paramf(2)=rphizf(1)*cos(rphizf(2))
           paramf(3)=rphizf(1)*sin(rphizf(2))
           paramf(1)=rphizf(3)
           paramf(4)=rphizf(4)
           paramf(5)=rphizf(5)+rphizf(2)
           if (paramf(5).gt.2*pi) paramf(5) = paramf(5)-2*pi
           if (paramf(5).lt.0.)   paramf(5) = paramf(5)+2*pi
           paramf(6)=rphizf(6)
       endif
      else                                       !xyz
        if(parami(4).eq.0.or.parami(6).eq.0)then
           exxerr=-1
        else
           do i=1,6
             rphiz(i)=parami(i)
             rphizf(i)=paramf(i)
           end do
           call fkxpla(rphiz,0,space(teindex(ite,idet)+12),
     1                 0.0,100000.0,0,rphizf,tmp,alrphi,exxerr)
        endif
        if(exxerr.ne.0)then
           if (fsfdbg.gt.1)then
             print*,'FSFXTS: FKXPLA failed. ierr =',exxerr
           endif
           call ucopy(parami,paramf,6)
           ierr=1
        else
c --- ME BUGFIX
           do i=1,6
             paramf(i)=rphizf(i)
           enddo
        endif
      endif
      end

