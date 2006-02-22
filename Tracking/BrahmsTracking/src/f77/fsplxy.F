*  
* $Id: fsplxy.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      SUBROUTINE FSPLXY(a2,ite,idet,xp,yp,hxyp,index)
*******************************************************************************
*
*      Polar inversion, xy version.
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
      real R2,a2,xp(*),yp(*),hxyp(3,detmax),dj11,dj12,dj22,R4I
      real xerr,yerr

c --- Polar Inversion with center at DELPHI origin:

      R2=x(ite,idet)**2+y(ite,idet)**2
      xp(index)=a2*x(ite,idet)/R2
      yp(index)=a2*y(ite,idet)/R2
CC      xerr=dx(ite,idet) !was: 1.0
CC      yerr=dy(ite,idet) !was: 1.0
      xerr=1.0
      yerr=1.0

      R4I        = 1./(R2*R2)
      DJ11       =(y(ITE,IDET)*y(ITE,IDET)-x(ite,idet)*x(ite,idet))*R4I
      DJ12       =-2.*x(ite,idet)*y(ITE,IDET)*R4I
      DJ22       = -DJ11
      HXYP(1,index) = DJ11*DJ11*xerr+DJ12*DJ12*yerr
CC     +               +2.*DJ11*DJ12*space(teindex(ite,idet)+18)
      HXYP(2,index) = DJ11*DJ12*xerr+DJ12*DJ22*yerr
CC     +               +(DJ11*DJ22+DJ12*DJ12)*space(teindex(ite,idet)+18)
      HXYP(3,index) = DJ22*DJ22*yerr+DJ12*DJ12*xerr
CC     +               +2.*DJ22*DJ12*space(teindex(ite,idet)+18)
CCc --- suppress error corelation
CC      HXYP(2,index) =  0.

      end

