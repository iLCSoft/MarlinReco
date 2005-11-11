*  
* $Id: fsftte.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsftte(idxv,detcode,detcode2,bad)
********************************************************************
*
* Test, if TE combination makes sense, reject funny ones by
* setting bad=.TRUE.
*
* INPUT:  idxv     - index array of TE combination for local arrays
*         detcode  - list of detector combinations to be rejected.
*                    Terminated by 0.
*         detcode2 - list of detector combinations where more TEs
*                    than detctors are required. Terminated by 0.
* OUTPUT: bad      - .TRUE. if killed
*                                                      Markus Elsing
********************************************************************
      implicit none
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
      integer detcode(*),detcode2(*),idxv(detmax)
      logical bad
      integer i,j,numdet,ibitted,numbit

      ibitted = 0
      bad     = .FALSE.
      numdet  = 0

c --- loop over idxv and create bitmask
      do i=1,detmax
c --- is detector used ?
        if (idxv(i).ne.0)then
          if (detid(i).gt.0)then
            call sbit1(ibitted,detid(i)+1)
            numdet = numdet + 1
          else ! this is the only left over solution, its i_STRNG ?
            ibitted = ispace(teindex(idxv(i),i)+2)
          endif
        endif
      enddo

c --- test if bad detector combination
      j = 1
      do while (detcode(j).ne.0)
        if(ibitted.eq.detcode(j))then
          if(fsfdbg.gt.0)then
            print*,'FSFTTE: bad detector combination rejected !'
          endif
          bad = .TRUE.
          return
        endif
        j = j + 1
      enddo

      if (numbit(ibitted).eq.numdet) then
c --- test detector combinations where more TEs than detectors
        j = 1
        do while (detcode2(j).ne.0)
          if(ibitted.eq.detcode2(j))then
c --- protect against multi point VFT
            if (idxv(i_ftd).eq.0) then
              bad = .TRUE.
            else if (idxv(i_vtx).eq.0) then
              bad = .TRUE.
            endif
          endif
          if (bad) then
            if(fsfdbg.gt.0)then
              print*,'FSFTTE: bad single combination rejected !'
            endif
            return
          endif
          j = j + 1
        enddo
      endif
      end

