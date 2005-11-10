*  
* $Id: fsfp2c.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfp2c(really,ite,idet,xyz)
*******************************************************************************
*
*     Transform polar coordinates to cartesian coordinates.
*     (if REALLY.EQ..FALSE. just copy.)
*
*      INPUT:  really - flag to toggle transformation
*              ite    - TE index in arrays x,y,z
*              idet   - detector index
*      OUTPUT: xyz    - xyz position
*                                                                  Daniel Wicke
*******************************************************************************
      implicit none
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
      logical really
      real xyz(3)
      integer ite,idet

      if(.NOT.really)then
         xyz(1)=x(ite,idet)
         xyz(2)=y(ite,idet)
         xyz(3)=z(ite,idet)
      else
         xyz(1)=R(ite,idet)*cos(RPHI(ite,idet)/R(ite,idet))
         xyz(2)=R(ite,idet)*sin(RPHI(ite,idet)/R(ite,idet))
         xyz(3)=z(ite,idet)
      endif

      end

