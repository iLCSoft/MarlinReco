*  
* $Id: fsfhelix.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfhelix(fsfpair,te1,det1,te2,det2)
*******************************************************************************
*
*      Do helix extrapolation and try to pair TE1 with TE2.
*
*      INPUT:  te1     - index of TE1 to extrapolate to TE2
*              det1    - detector index
*              te2     - index of second TE
*              det2    - detector index
*      OUTPUT: fsfpair - .TRUE. for good pair
*                                                                 Markus Elsing
*******************************************************************************
      implicit none
#include "include/fsfb.inc"
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfxerrs.inc"
#include "include/fsfpi.inc"
      logical fsfpair
      integer te1,det1,te2,det2,datdescr1,datdescr2,jbit
      real    parami(6),paramf(6),xyz(3)
      integer ierr,ij

      logical new_extrapolation
      integer te1_old,det1_old,det2_old
      real    old_value

      logical first
      data    first /.true./
      save first

* --- save cached values for next call (for safety)
      save te1_old,det1_old,det2_old,old_value,parami,paramf

      fsfpair=.FALSE.
      datdescr1=datdescr(te1,det1)
      datdescr2=datdescr(te2,det2)

* --- initialise extrapolation cache
      if (first) then
         first=.false.
         te1_old=0
      end if

* --- do we need to reextrapolated, otherwise use old result

      new_extrapolation = .TRUE.

* --- same detector+TE to start and target is same detector ?
      if (te1.eq.te1_old.and.det1.eq.det1_old.and.det2.eq.det2_old.and.
* --- check also phi for example to be sure not to miss the next event
     1    abs(parami(5)-phi(te1,det1)).lt.0.0001) then

* --- in polar coordinates or not ?
        if (jbit(datdescr2,1).eq.1.) then
* --- cut on 1um in R
          if (abs(R(te2,det2)-old_value).lt.0.0001) then
            if(fsfdbg.gt.1)then
              print *,'FSFHELIX: cashed, no extapolation !'
            endif
            new_extrapolation = .FALSE.
          endif
        else
* --- cut on 1um in Z
          if (abs(Z(te2,det2)-old_value).lt.0.0001) then
            if(fsfdbg.gt.1)then
              print *,'FSFHELIX: cashed, no extapolation !'
            endif
            new_extrapolation = .FALSE.
          endif
        endif
      endif

      if (new_extrapolation) then

        if(fsfdbg.gt.1)then
          print *,'FSFHELIX: -> new extapolation to be done !'
        endif

c --- check if needed values are measured, if not: return
        if(0.eq.jbit(datdescr1,3) .or.
     1     0.eq.jbit(datdescr1,4) .or.
     2     0.eq.jbit(datdescr1,5) .or.
     3     0.eq.jbit(datdescr1,8) .or.
     4     0.eq.jbit(datdescr1,9) .or.
     5    (0.eq.jbit(datdescr1,10).and.0.eq.jbit(datdescr1,11)).or.
     6    (0.eq.jbit(datdescr2,3).and.0.eq.jbit(datdescr2,5)))then
          if(fsfdbg.gt.1)then
            print *,'FSFHELIX: Not all needed coordinates measured.'
          endif
          return
        endif

c --- convert coordinates of first TE to cartesian if neccessary.
        call fsfp2c(1.eq.jbit(datdescr1,1),te1,det1,xyz)

c --- set startparameters for helix-extrapolation

        parami(1)=xyz(3)        !z
        parami(2)=xyz(1)        !x
        parami(3)=xyz(2)        !y
        parami(4)=theta(te1,det1)
        parami(5)=phi(te1,det1)
        parami(6)=overrt(te1,det1)

c --- extrapolate
        if(parami(4).eq.0.or.parami(6).eq.0)then
          ierr=-1
        else
          CALL fsfxts (parami,paramf,te2,det2,ierr)
        endif

        if(ierr.ne.0)then
          if(fsfdbg.gt.1)then
            print*,'FSFHELIX: helix-extrapolation failed'
          endif
          if(.NOT.(ierr.eq.1.or.ierr.eq.3.or.ierr.eq.-1))then
            print*,'WARNING in FSFHELIX: FSFXTS returned unknown error.'
          endif
          return
        endif

* --- save important values for checking next time the cashing
        te1_old   = te1
        det1_old  = det1
        det2_old  = det2
        if (jbit(datdescr2,1).eq.1.) then
          old_value = R(te2,det2)
        else
          old_value = Z(te2,det2)
        endif

      endif

#if defined(HBOOK)
      if (fsfhist.gt.0) then
        ij=10000*det1+50*det2
c --- labels ok ?
        if(telabl(te1,det1).ne.0.and.
     1    telabl(te1,det1).eq.abs(telabl(te2,det2)))then ! det1 needs labl>0 !
          ij = ij + 10
        else
          ij = ij + 20
        endif
        CALL fsfhst (paramf,te1,det1,te2,det2,ij)
      endif
#endif

c --- now call test routine to compare extrapolation and TEs
        CALL fsfcomp (paramf,te1,det1,te2,det2,
     1                xcx,xcy,xcRPhiR,xcz,
     2                xctheta,xcphi,xcthest,xcphist,
     3                fsfpair)
      END

