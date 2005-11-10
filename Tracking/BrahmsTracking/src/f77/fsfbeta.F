*  
* $Id: fsfbeta.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfbeta (phi2,RPhiR,beta)
*******************************************************************************
*
*     Calculate beta = phi - RPhi/R
*
*      INPUT:  phi2   - direction at R
*              RPhiR  - Phi from position
*      OUTPUT: beta   - phi-Phi
*                                                                 Markus Elsing
*******************************************************************************
      implicit none
#include "include/fsfpi.inc"
      double precision phi,RPhiR,beta,phi2

      phi=phi2
      if (phi2.lt.0.)   phi   = phi  +2.*pi
      if (RphiR.lt.0.) RphiR = RphiR+2.*pi

c --- ME bugfix !!!
      beta = dmod(phi-RphiR+dble(2.*pi+pi),dble(2.*pi))-dble(pi)

CC      beta = phi - RphiR
CC      if (beta.gt.pi)  beta = 2.*pi-beta
CC      if (beta.lt.-pi) beta = 2.*pi+beta
      end

