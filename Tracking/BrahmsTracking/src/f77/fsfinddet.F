*  
* $Id: fsfinddet.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfinddet(telist,tes,idet,ite,iteindex)
**********************************************************************
*
* Find TE from detector idet in telist. The index of the TE in the
* TS list is returned.
*
* INPUT:  telist - array of TE identifier
*         tes    - number of TEs in list
*         idet   - i_TPC=TPC,i_FTD=FTD,i_ITC=ITC...
* OUTPUT: ite      - index of the TE found in list of this TS
*         iteindex - index of TE in local TE arrarys
*                                                        Markus Elsing
**********************************************************************
      implicit none
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
      integer k,l,ite,idet,tes,iteindex,telist(*)
      external evtese
      logical  evtese
c
      ite      = 0
      iteindex = 0
c
c --- scan internal arrays to find full information
      do k = 1,tes
        do l = 1,teindexused(idet)
          if(telist(k).eq.teid(l,idet)) then
            ite      = k
            iteindex = l
            return
          endif
        enddo
      enddo
      end

