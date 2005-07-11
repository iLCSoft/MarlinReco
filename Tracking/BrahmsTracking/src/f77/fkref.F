*  
* $Id: fkref.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE FKREF(NTE,DTE,ref,irf)
      implicit none
      save
#include "include/coildims.inc"
#include "include/fkparm.inc"
c#include "include/gconst.inc"
#include "include/alcons.inc"

      integer nte
      integer jbit
      integer i,j
      real    DTE(mxte,NTE),ref(*)
      integer irf(*)
      real x1, y1, z1, x, y, z
      real xp(ntemx),yp(ntemx),hxy(3,ntemx)
      real r2, thm, ct, pt, dphi, axy(2),chixy,haxy(3)
      integer imts
      parameter (imts = 1+8)
*
CSJAWRT
c      write(*,*) "hello I'm in fkref @$@#$@#$@#$@#$#@$@#$@#$@#$"
      thm = 0
      do i=1,nte
         if (jbit(dte(4,i),1).eq.0) then
            x     = dte(10,i)
            y     = dte(11,i)
            z     = dte(12,i)
         else
            x     = cos(dte(11,i)/dte(10,i))*dte(10,i)
            y     = sin(dte(11,i)/dte(10,i))*dte(10,i)
            z     = dte(12,i)
         endif
*
         if (i.eq.1) then
            x1    = x
            y1    = y
            z1    = z
            j     = 0            !(KH)
         else
     >        if ((z.ne.z1).and.((x.ne.x1).or.(y.ne.y1))) then  !(KH)
C(KH)       j     = i-1
            j=j+1                !(KH)
            x     = x-x1
            y     = y-y1
            z     = z-z1
            r2    = x**2+y**2
            xp(j) = x/r2
            yp(j) = y/r2
            hxy(1,j) = (0.01/r2)**2
            hxy(2,j) = 0.
            hxy(3,j) = (0.01/r2)**2
            thm      = thm+sqrt(r2)/z
         endif
      enddo
C      call fsflft(xp,yp,hxy,axy,haxy,2,nte-1,chixy)
      call fsflft(xp,yp,hxy,axy,haxy,2,j,chixy) !(KH)
      CT     = 1.+AXY(1)**2
      PT     = 0.5*SQRT(CT)*BFIELD*FCT/AXY(2)
      thm = atan(thm/real(j))
      if (thm.lt.0) thm = thm+pi
      irf(1)=0
      irf(2)=0
      irf(3)=imts
      ref(17) = x1
      ref(18) = y1
      ref(19) = z1
      ref(20) = thm
      ref(21) = ATAN(AXY(1))
      ref(22) = sin(thm)/pt
      call fsfand(atan2(yp(1),xp(1)),ref(21),dphi,i)
      if (abs(dphi).gt.piby2) then
         ref(21) = ref(21) +PI
         ref(22) = -ref(22)
      endif
      if (ref(21).lt.0.) ref(21) = ref(21) + twopi
      if (ref(21).gt.twopi) ref(21) = ref(21) - twopi
*
  999 return
      end


