*  
* $Id: fsftts.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsftts(its,detcode,detcode2)
********************************************************************
*
* Test, if TS makes sense, reject funny ones by setting ITS=0.
*
* INPUT:  its      - index of TS to be tested
*         detcode  - list of detector combinations to be rejected.
*                    Terminated by 0.
*         detcode2 - list of detector combinations where more TEs
*                    than detctors are required. Terminated by 0.
* OUTPUT: its      - zero if killed
*                                                      Markus Elsing
********************************************************************
      implicit none
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfbits.inc"
      integer its,j,detcode(*),detcode2(*)
      integer iftd,ivtx,iteindex,numbit,idetus
      logical fsfftdm,fsfvtxm

      if(its.eq.0) return

      idetus = ispace(its+2)
c --- set down pseudoTE bits (FMIVONLY)
      call sbit0(idetus,detid_pvtx  +1)
      call sbit0(idetus,detid_zfudge+1)
CC detid_shower is TPC !!!
CC      call sbit0(idetus,detid_shower+1)

      if (numbit(idetus).gt.ispace(its+5)) then
        print*,'FSFTTS: corruption in TS detector word !'
        its=0
        return
      endif

c --- test if bad detector combination
      j = 1
      do while (detcode(j).ne.0)
        if(idetus.eq.detcode(j))then
          if(fsfdbg.gt.0)then
            print*,'FSFTTS: bad detector combination rejected !'
          endif
          its=0
          return
        endif
        j = j + 1
      enddo

      if (numbit(idetus).eq.ispace(its+5)) then
c --- test detector combinations where more TEs than detctors
        j = 1
        do while (detcode2(j).ne.0)
          if(idetus.eq.detcode2(j))then
c --- protect against multi point FTD
            call fsfgette(its,i_FTD,iftd,iteindex)
            if(iftd.eq.0)then
              its = 0
            else
c --- submodule should do the job ?
              if (.not.fsfftdm(ispace(teindex(iteindex,i_FTD)+1)))
     1          its = 0
            endif
c --- protect against multi point VTX
            call fsfgette(its,i_VTX,ivtx,iteindex)
            if(ivtx.eq.0)then
              its = 0
            else
c --- submodule should do the job ?
              if (.not.fsfvtxm(ispace(teindex(iteindex,i_VTX)+1)))
     1          its = 0
            endif
          endif
          if (its.eq.0) then
            if(fsfdbg.gt.0)then
              print*,'FSFTTS: bad single combination rejected !'
            endif
            return
          endif
          j = j + 1
        enddo
      endif
      end

