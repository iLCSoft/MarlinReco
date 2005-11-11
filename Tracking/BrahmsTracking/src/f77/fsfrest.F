*  
* $Id: fsfrest.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      SUBROUTINE FSFREST (pd)
*******************************************************************************
*
*     Create pseudo TEs for all strings to search for left over TEs.
*
*     INPUT: pd - primary detector name
*
*                                                                Markus Elsing
*******************************************************************************
      IMPLICIT NONE
#include "include/fsfparams.inc"
#include "include/fsfflags.inc"
#include "include/brcntl.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfpi.inc"
#include "include/fsfsim.inc"
      integer its,ndat,ilabel,ishower,ierr,pd
      logical fsfgetspace
c --- loop over all TS
      do its = 1,tsindexused
c --- still active ?
        if (tsindex(its).ne.0) then
c --- test space array !
          if(.not.fsfgetspace('FSFREST ',tssize+maxtev))RETURN
c --- OK, create TE from it
          call fsftste (tsindex(its),spaceused,
     1                  tssize+ispace(tsindex(its)+5),ndat,ishower)
c --- valid TS if ndat > 0
          if (ndat.eq.0) then
            print*,'WARNING in FSFREST: failed to create dummy TE'
          else
c --- now add in all TEs in this string
            call ucopy (space(tsindex(its)+tssize+1),
     1                  space(spaceused+ndat+1),ispace(tsindex(its)+5))
            ndat = ndat+ispace(tsindex(its)+5)
c --- update TE pointers and indices, TEID = tsindex(its) !!!
            call fsfwte ('FSFREST',pd,tsindex(its),ndat,ierr)
c --- keep shower TE information
            showerte(teindexused(pd),pd)=ishower
c --- get label
            if (delsim) CALL FSFTLLA (ispace(tsindex(its)+5),
     1                        ispace(tsindex(its)+tssize+1),ilabel)
            telabl(teindexused(pd),pd)=ilabel
            if(fsfdbg.gt.0)then
              print*,'FSFREST: pseudo TE created from TS of type ',
     1               ispace(tsindex(its)+1),' label ',ilabel
            endif
          endif
        endif
      enddo

      end


