*  
* $Id: fsftsdump.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsftsdump
*******************************************************************************
*
*     Print short list of all TS in memory
*
*
*     WRITTEN BY                                     Kristian Harder
*
*******************************************************************************
      implicit none
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
      integer i,j,isave

      write(6,*)
      write(6,*) 'FSFTSDUMP: current list of TSs:'
      isave=0
      do i=1,tsindexused
        if (tsindex(i).gt.0) then
          isave=isave+1
          write(6,*) 'FSFTSDUMP: TS ',i,': ',ispace(tsindex(i)+5),' TEs'
          write(6,*) '        detector code ',ispace(tsindex(i)+2)
          do j=1,ispace(tsindex(i)+5)
            write(6,*) '        TE ',ispace(tsindex(i)+tssize+j)
          end do
        end if
      end do
      write(6,*) 'total number of TS: ',isave

      return
      end

