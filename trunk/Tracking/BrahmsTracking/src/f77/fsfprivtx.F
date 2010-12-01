*  
* $Id: fsfprivtx.F,v 1.1 2005-08-02 09:32:40 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfprivtx
*******************************************************************************
*
*      Fills FMPVTX sequence with primary vertex fit information.
*
*      WRITTEN BY                                  Markus Elsing
*      CHANGES BY                                  Kristian Harder, August 99
*******************************************************************************
*
      implicit none
#include "include/brcntl.inc"
#include "include/fsfclun.inc"
#include "include/fmpvtx.inc"
#include "include/ipsmear.inc"
      INTEGER I

c --- fill up local variables array with copy of beam spot
      pvtx(1)    = 0.          ! x position of vertex
      pvtx(2)    = 0.          ! y position of vertex
      pvtx(3)    = 0.          ! z position of vertex
      dpvtx(1)   = SMRIPX**2   ! error in x (from BRAHMS -KH-)
      dpvtx(2)   = 0.0
      dpvtx(3)   = SMRIPY**2   ! error in y
      dpvtx(4)   = 0.0
      dpvtx(5)   = 0.0
      dpvtx(6)   = SMRIPZ**2   ! error in z
      if(fsfdbg.gt.0)then
        write(iprnt,900) (pvtx(i),i=1,3)
 900    format (1x,'FSFPRIVTX: assume beam spot at : ',3(f6.3,' '))
        write(iprnt,901) (sqrt(dpvtx(i)),i=1,3)
 901    format (1x,'FSFPRIVTX: beam spot errors : ',3(f10.7,' '))
        write(iprnt,902) (sqrt(dpvtx(i)),i=4,6)
 902    format (1x,'FSFPRIVTX:                    ',3(f10.7,' '))
      endif
      end

