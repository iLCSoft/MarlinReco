*  
* $Id: fsfgette.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfgette(its,idet,ite,iteindex)
**********************************************************************
*
* Find TE from detector idet in string its. The index of the TE in the
* TS list is returned.
*
* INPUT:  its  - tsindex
*         idet - i_TPC=TPC,i_FTD=FTD,i_ITC=ITC,...
* OUTPUT: ite      - index of the TE found in list of this TS
*         iteindex - index of TE in local TE arrarys
*                                                        Markus Elsing
**********************************************************************
      implicit none
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
      integer ite,idet,its,iteindex,jbit
c
      ite      = 0
      iteindex = 0
c --- first test bits !
      if (jbit(ispace(its+2),detid(idet)+1).eq.0) RETURN
c --- call routine to scan te list
      call fsfinddet(ispace(its+tssize+1),ispace(its+5),idet,
     1               ite,iteindex)
      end

