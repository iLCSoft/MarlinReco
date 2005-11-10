*  
* $Id: fsfgetspace.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      logical function fsfgetspace(name,iwords)
*******************************************************************************
*
*      tests, if space is available in space array.
*
*      INPUT:  name(10) - name of calling routine
*              iwords   - number of words
*                                                                 Markus Elsing
*******************************************************************************
      implicit none
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
      CHARACTER*(10) name
      INTEGER        iwords
*
      if (spaceused+iwords.gt.spacesize) then
        print*,'ERROR in ',name,': ',
     1         '**** Not enough space in array ****'
        print*,'            Data not used for processing'
        print*,'            Enlarge parameter spacesize, recompile,'
        print*,'            and rerun to cure this problem!!'
        fsfgetspace = .FALSE.
        CALL TKVETO(1)
      else
        fsfgetspace = .TRUE.
* ---   update maxspaceused for statistics
        if(spaceused+iwords.gt.maxspaceused)then
           maxspaceused=spaceused+iwords
        endif
      endif
      end

