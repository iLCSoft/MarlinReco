*  
* $Id: fsfwte.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfwte(name,detindex,idv,ndat,ierr)
*******************************************************************************
*
*      put information about TE into arrays and update space array.
*
*      INPUT:  name(10) - name of calling routine
*              detindex - detector index in DELSEARCH convention
*              idv      - tanangra identifer
*              ndat     - length of array
*      OUTPUT: ierr     - 1 for error
*                                                                 Markus Elsing
*******************************************************************************
      implicit none
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
      CHARACTER*(10) name
      INTEGER        detindex,idv,ndat,ierr
*
      if(teindexused(detindex).lt.temax)then
        fsfnte(detindex)                       =fsfnte(detindex) +1
        teindexused(detindex)                  =teindexused(detindex)+1
        teindex(teindexused(detindex),detindex)=spaceused
        teid   (teindexused(detindex),detindex)=idv
        tesize (teindexused(detindex),detindex)=ndat
        spaceused                              =spaceused+ndat
        if (spaceused.gt.spacesize) then
          STOP 'FATAL error in FSFWTE: space array corrupted'
        endif
        ierr = 0
      else
        print*,'ERROR in ',name,': ',
     1         '**** TEMAX to small to read this TE ****'
        print*,'            Enlarge parameter temax recompile'
        print*,'            and rerun to cure this problem!!'
        print*,'detindex =',detindex
        CALL TKVETO(1)
        ierr = 1
      endif
      end

