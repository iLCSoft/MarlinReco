*  
* $Id: fsfnsub.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfnsub(flags)
**********************************************************************
*
* Kill substrings in TS list to make things more easy afterwards.
*
*     INPUTL
*      flags - bit 0 set: remove doubles    (use parameter FLAG_DBL)
*              bit 1 set: remove substrings (use parameter FLAG_SUB)
*
*                                                        Markus Elsing
**********************************************************************
      implicit none
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
      integer i,j,k,l,nkill,ndble,nsub,flags,jbit
      logical tefound,different
      integer nstart,nend

c --- loop over ts and remove sub-strings
c
      nkill=0
      ndble=0
      nsub =0

      if(fsfdbg.gt.0)then
         if(jbit(flags,1).eq.1) print*,'FSFNSUB: remove double TS'
         if(jbit(flags,2).eq.1) print*,'FSFNSUB: remove sub-TS'
      endif

c --- final number of TS active
      if(fsfdbg.gt.0)then
        nstart = 0
        do i=1,tsindexused
          if(tsindex(i).ne.0) nstart = nstart+1
        enddo
      endif

      do i=1,tsindexused
c --- already killed ?
        if(tsindex(i).ne.0)then
c --- second loop
          do j=1,tsindexused
c --- already killed ?
            if(tsindex(j).ne.0.and.i.ne.j)then
c --- TS i has more TEs than TS j
              if(ispace(tsindex(i)+5).ge.ispace(tsindex(j)+5))then
                different=.false.
c --- loop over TEs
                do k=1,ispace(tsindex(j)+5)
                  tefound=.false.
                  do l=1,ispace(tsindex(i)+5)
c --- same TE found ?
                    if (ispace(tsindex(j)+tssize+k).eq.
     1                  ispace(tsindex(i)+tssize+l))   tefound=.true.
                  enddo
c --- if TE not found, then we have different TS
                  if(.not.tefound) different=.true.
                enddo
c --- see FSFSTS, this kills TS
                if(.not.different)then
                  if(jbit(flags,1).eq.1.and.
     1               ispace(tsindex(i)+5).eq.ispace(tsindex(j)+5))then
                    if(fsfdbg.gt.0)then
                      print'(A,4(A,i5))',' FSFNSUB: double TS found,',
     1                     ' TS ',i,' module:',ispace(tsindex(i)+1),
     2                     ' TS ',j,' module:',ispace(tsindex(j)+1)
                    endif
                    if(fsfdbg.gt.1)then
                      do k=1,ispace(tsindex(i)+5)
                        print'(A,i8)','   TE-list 1: ',
     1                       ispace(tsindex(i)+tssize+k)
                      enddo
                      do k=1,ispace(tsindex(j)+5)
                        print'(A,i8)','   TE-list 2: ',
     1                       ispace(tsindex(j)+tssize+k)
                      enddo
                    endif
                    tsindex(j) = 0
                    nkill      = nkill+1
                    ndble      = ndble+1
                  else if(jbit(flags,2).eq.1.and.
     1               ispace(tsindex(i)+5).gt.ispace(tsindex(j)+5))then
                    if(fsfdbg.gt.0)then
                      print'(A,4(A,i5))',' FSFNSUB: sub-TS found,',
     1                     ' TS ',i,' module:',ispace(tsindex(i)+1),
     2                     ' sub-TS ',j,' module:',ispace(tsindex(j)+1)
                    endif
                    if(fsfdbg.gt.1)then
                      do k=1,ispace(tsindex(i)+5)
                        print'(A,i8)','   TE-list:   ',
     1                       ispace(tsindex(i)+tssize+k)
                      enddo
                      do k=1,ispace(tsindex(j)+5)
                        print'(A,i8)','   subTE-list:',
     1                       ispace(tsindex(j)+tssize+k)
                      enddo
                    endif
                    tsindex(j) = 0
                    nkill      = nkill+1
                    nsub       = nsub+1
                  endif
                endif
              endif
            endif
          enddo
        endif
      enddo

c --- final number of TS active
      if(fsfdbg.gt.0)then
         nend = 0
         do i=1,tsindexused
           if(tsindex(i).ne.0) nend = nend+1
         enddo
c
         print*,'FSFNSUB: number of TS at start : ',nstart
         print*,'         number of TS dropped  : ',nkill
         print*,'         double TS found       : ',ndble
         print*,'         sub TS found          : ',nsub
         print*,'         number of TS at end   : ',nend
      endif
      end

