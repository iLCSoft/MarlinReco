*  
* $Id: fsfrefit.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfrefit
*******************************************************************************
*
*      Final refit without error scaling to cut down statistics.
*
*******************************************************************************
      IMPLICIT NONE

#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/brcntl.inc"
#include "include/fsfflags.inc"
#include "include/fkdebug.inc"

      INTEGER its,nrej,nfit,ilabel

      if (fsfdbg.gt.0) then
        print *,'FSFREFIT: --------------------------------',
     +          '----------------'
        print *,'FSFREFIT: final refit of track candidates,',
     +          ' no errorscaling'
      endif

c --- some statistics
      nfit = 0
      nrej = 0

c --- loop over TSlist and refit final tracks
      do its = 1,tsindexused
        if (tsindex(its).ne.0) then
          if (fsfdbg.gt.0) then
            call fsftlla(ispace(tsindex(its)+5),
     1                  ispace(tsindex(its)+tssize+1),ilabel)
            print *,'FSFREFIT: refit TS #',its,' => label ',ilabel
          endif
          nfit = nfit + 1
c --- refit without error scaling and with outlayer to remove bad hits
          call fsfits (tsindex(its),FLAG_NOSCALING+FLAG_OUTLAYER
c          write(6,*) 'fsfrefit: warning: using error scaling!'
c          call fsfits (tsindex(its),FLAG_OUTLAYER
     1                 +FLAG_EXTRAPOL,ispace(tsindex(its)+1)-3000)
          if (tsindex(its).eq.0) then
            nrej = nrej + 1
            if (fsfdbg.gt.0) then
              print *,'FSFREFIT: => TS rejected !'
            endif
          endif
        endif
      enddo

      if (fsfdbg.gt.0) then
        print *,'FSFREFIT: number of tracks refitted:',nfit
        print *,'FSFREFIT: number of tracks rejected:',nrej
        print *,'FSFREFIT: --------------------------------',
     +          '----------------'
      endif

      end


