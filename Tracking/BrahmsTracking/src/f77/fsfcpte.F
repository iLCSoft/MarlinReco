*  
* $Id: fsfcpte.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfcpte(ite,idet,helixpar)
*******************************************************************************
*
*     Copy TE parameter to helixpar array.
*
*     INPUT:  ite      - nindex of TE
*             idet     - detector index
*     OUTPUT: helixpar - parameters
*                                                                  Daniel Wicke
*******************************************************************************
      implicit none
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/brcntl.inc"
#include "include/fsfflags.inc"
#include "include/fsfreject.inc"
#include "include/fsfpi.inc"
      integer ite,idet
      real helixpar(6)
      integer jbit

      if(1.eq.jbit(datdescr(ite,idet),1))then !RRPhiz
        if(R(ite,idet).eq.0)then
          helixpar(2)=0
          helixpar(3)=0
        else
          helixpar(2)=R(ite,idet)*cos(RPhi(ite,idet)/R(ite,idet))
          helixpar(3)=R(ite,idet)*sin(RPhi(ite,idet)/R(ite,idet))
        endif
      else
        helixpar(2)=x(ite,idet)
        helixpar(3)=y(ite,idet)
      endif
      helixpar(1)=z(ite,idet)
      helixpar(4)=theta(ite,idet)
      helixpar(5)=phi(ite,idet)
      helixpar(6)=overrt(ite,idet)

      end


