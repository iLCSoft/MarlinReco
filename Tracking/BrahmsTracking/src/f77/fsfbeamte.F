*  
* $Id: fsfbeamte.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfbeamte(datte,n)
*******************************************************************************
*
*     Create beam spot x,y TE at Z beam spot.
*
*     Input: DATTE(,) - TER data to be fitted
*            n        - number of te in datte
*
*                                                                 Markus Elsing
*******************************************************************************
      IMPLICIT NONE
#include "include/exmach.inc"
#include "include/ipsmear.inc"
#include "include/fsfparams.inc"
#include "include/fsfbits.inc"

      integer n
      real    DATTE(MXTE,NTEMX)
      integer itenew(mxte)
      real    tenew(mxte)
      equivalence (itenew,tenew)

      call vzero(itenew,MXTE)

      itenew(1)=detid_pvtx ! interaction region (tanagra index)
      itenew(4)=124 ! measured xyz, xy with errors
      itenew(5)=20  ! pointer to massid (dummy)
      itenew(6)=3   ! charge unknown
      tenew(10)=bspotx
      tenew(11)=bspoty
      tenew(12)=bspotz
      tenew(17)=smripx**2
      tenew(19)=smripy**2

c --- increase number of TEs
      n=n+1
      call ucopy (tenew,datte(1,n),mxte)
      end


