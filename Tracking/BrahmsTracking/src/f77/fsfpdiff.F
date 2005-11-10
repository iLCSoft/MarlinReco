*  
* $Id: fsfpdiff.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      function fsfpdiff(phi1,phi2)
*************************************************************************
*                                                            Daniel Wicke
*************************************************************************
      implicit none
#include "include/fsfpi.inc"
      real fsfpdiff,phi1,phi2

      fsfpdiff=phi1-phi2
      !* we need the correct sign
      if(fsfpdiff.gt.Pi)then
         fsfpdiff=fsfpdiff-2*Pi
      endif
      if(fsfpdiff.lt.-Pi)then
         fsfpdiff=fsfpdiff+2*Pi
      endif
      return
      end

