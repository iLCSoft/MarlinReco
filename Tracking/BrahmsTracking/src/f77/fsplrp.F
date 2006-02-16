*  
* $Id: fsplrp.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsplrp(a2,ite,idet,xp,yp,hxyp,index)
*******************************************************************************
*
*      Polar inversion.
*
*      INPUT:  a2    - ?
*              ite   - TE in list to be transformed
*              idet  - detector number to be worked on
*      OUTPUT: xp    - x transformed
*              yp    - y transformed
*              index - index in x and y array
*              hxyz  - error matrix (3 elements)
*                                                                  Daniel Wicke
*******************************************************************************
      implicit none
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
      integer ite,idet,index
      real PhiP,C,S,a2,xp(*),yp(*),hxyp(3,detmax),R4I,xerr
      real DJ11,DJ12,DJ22,R2I,yerr

c --- Polar Inversion with center at DELPHI origin:

      PhiP=RPhi(ite,idet)/R(ite,idet)
      C=cos(PhiP)
      S=sin(PhiP)
      xp(index)=a2*c/R(ite,idet)
      yp(index)=a2*s/R(ite,idet)

CC      xerr=abs(-S*dRPhi(ite,idet))
CC      yerr=C*dRPhi(ite,idet)
CC      xerr=abs(xerr)
      xerr=1.0
      yerr=1.0

      R2I           = 1/R(ite,idet)**2
      R4I           = 1/R(ite,idet)**4
CC      HXYP(1,index) =  S*S*R4I*xerr
CC      HXYP(2,index) = -S*C*R4I*xerr
CC      HXYP(3,index) =  C*C*R4I*xerr
CCc --- suppress error corelation
CC      HXYP(2,index) =  0.

      DJ11       =(S*S-C*C)*R2I
      DJ12       =-2.*C*S*R2I
      DJ22       = -DJ11
      HXYP(1,index) = DJ11*DJ11*xerr+DJ12*DJ12*yerr
CC     +               +2.*DJ11*DJ12*space(teindex(ite,idet)+18)
      HXYP(2,index) = DJ11*DJ12*xerr+DJ12*DJ22*yerr
CC     +               +(DJ11*DJ22+DJ12*DJ12)*space(teindex(ite,idet)+18)
      HXYP(3,index) = DJ22*DJ22*yerr+DJ12*DJ12*xerr
CC     +               +2.*DJ22*DJ12*space(teindex(ite,idet)+18)

      end

