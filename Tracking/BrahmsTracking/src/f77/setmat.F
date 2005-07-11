*  
* $Id: setmat.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      subroutine setmat
C****************************************************************************
c* Version number: 0.0      Date of last modification  30.11.98             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
c*------------------------------------------------------------------------- *
C* Main Authors:                                                            *
c* Name          Institute                        Contact email             *
C* K. Moenig     DESY-Zeuthen                     moenig@ifh.de             *
c*                                                                          *
c****************************************************************************
C  set up the material database for the track fit
C ---------------------------------------------------------------------------
      implicit none
      save
#include "include/coildims.inc"
c#include "turnon.inc"
c#include "pixpipe.inc"
c#include "ccdpipe.inc"
c#include "vxdpix.inc"
c#include "vxdccd.inc"
c#include "itcdims.inc"
#include "include/tpcdims.inc"
c#include "ftddims.inc"
c#include "sitdims.inc"
c#include "fchdims.inc"
#include "include/fkparm.inc"
#include "include/fkddes.inc"
#include "include/fkexts.inc"
#include "include/fkfild.inc"
#include "include/fkcons.inc"
c#include "flag.inc"
*
      integer  nexfch
      common /extr/nexfch
*
      real xrbery, xrsili, xralu, xrargon, xrstrw, xrkapton, xrcopper
      parameter (xrbery=35.28, xrsili=9.36, xralu=8.9, xrargon=10971.)
      parameter (xrkapton=28.6,xrstrw=363.,xrcopper=1.43)
      real dedxbery, dedxsili, dedxalu,dedxargon, dedxstrw, dedxkapton
      real dedxcopper
      parameter (dedxbery=1.85*1.59e-3, dedxsili=2.33*1.66e-3)
      parameter (dedxalu=2.70*1.62e-3, dedxargon=0.0018*1.52e-3)
      parameter (dedxkapton=1.42*1.7e-3,dedxstrw=0.16*1.7e-3)
      parameter (dedxcopper=8.96*1.6e-3)
      integer nsttpc,nstitc
      parameter (nsttpc=10,nstitc=4)
      real xstep

      real xs
      integer i
*
*     magnetic field
*
      consb = bfield*fct

      ncmat = 0
      npmat = 0
      nexs  = 0
ccc* beam pipe
ccc      if (lvxdp) then
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = bpp0(1)
ccc         zcmin(ncmat)  = -bpp0(3)
ccc         zcmax(ncmat)  = bpp0(3)
ccc         xrlc(ncmat)   = (bpp0(2)-bpp0(1))/xrbery
ccc         xelosc(ncmat) = (bpp0(2)-bpp0(1))*dedxbery
ccc*
ccc*          radius of 1st conical part is inner radius
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = bpp1(2)
ccc         zcmin(ncmat)  = bpp0(3)
ccc         zcmax(ncmat)  = bpp0(3)+2.*bpp1(1)
ccc* fudge factor is the longer path through the cone for 8 deg
ccc         xrlc(ncmat)   = (bpp1(3)-bpp1(2))/xrbery*
ccc     $        (sin(0.14)/sin(0.14-atan(0.5*(bpp1(4)-bpp1(2))/bpp1(1))))
ccc         xelosc(ncmat) = (bpp1(3)-bpp1(2))*dedxbery*
ccc     $        (sin(0.14)/sin(0.14-atan(0.5*(bpp1(4)-bpp1(2))/bpp1(1))))
ccc*
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = rcmat(ncmat-1)
ccc         zcmin(ncmat)  = -zcmax(ncmat-1)
ccc         zcmax(ncmat)  = -zcmin(ncmat-1)
ccc         xrlc(ncmat)   = xrlc(ncmat-1)
ccc         xelosc(ncmat) = xelosc(ncmat-1)
ccc*
ccc      endif
ccc*
ccc      if (lvxdc) then
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = BPCINR1
ccc         zcmin(ncmat)  = -BPCINZ1
ccc         zcmax(ncmat)  = BPCINZ1
ccc         xrlc(ncmat)   = BPCT1/xrbery
ccc         xelosc(ncmat) = BPCT1*dedxbery
ccc*
ccc*          radius of 1st conical part is middle of cone
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = 0.5*(BPCINR1+BPCINR2)
ccc         zcmin(ncmat)  = bpcinz1*rcmat(ncmat)/bpcinr1
cccc         zcmax(ncmat)  = bpcinz3*rcmat(ncmat)/bpcinr2
ccc         zcmax(ncmat)  = 100.
ccc* fudge factor is the longer path through the cone for 10 deg
ccc         xrlc(ncmat)   = bpct2/xrbery*
ccc     $        (sin(0.175)/sin(0.175-atan((BPCINR2-BPCINR1)/
ccc     $        (BPCINZ2-BPCINZ1))))
ccc         xelosc(ncmat) = bpct2*dedxbery*
ccc     $        (sin(0.175)/sin(0.175-atan((BPCINR2-BPCINR1)/
ccc     $        (BPCINZ2-BPCINZ1))))
ccc*
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = rcmat(ncmat-1)
ccc         zcmin(ncmat)  = -zcmax(ncmat-1)
ccc         zcmax(ncmat)  = -zcmin(ncmat-1)
ccc         xrlc(ncmat)   = xrlc(ncmat-1)
ccc         xelosc(ncmat) = xelosc(ncmat-1)
ccc*
ccc
ccc      endif
ccc*

      nexs          = nexs+1
      typx(nexs)    = 'CYLI'
C FIXME:SJA: this line is commened out and the value replace with 0 to ensure that the material
C     description of brahm and marlin "agree"
c      rzsurf(nexs)  = rcmat(1)-0.01
      rzsurf(nexs)  = 0.
      zrmin(nexs)   = -1000.
      zrmax(nexs)   = 1000.


ccc*  microvertex
ccc      if (lvxdp) then
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = vtxd1(1)+0.01
ccc         zcmin(ncmat)  = -vtxd1(3)
ccc         zcmax(ncmat)  = vtxd1(3)
ccc         xrlc(ncmat)   = (vtxd1(2)-vtxd1(1))/xrsili
ccc         xelosc(ncmat) = (vtxd1(2)-vtxd1(1))*dedxsili
ccc*
ccc         nexs          = nexs+1
ccc         typx(nexs)    = 'CYLI'
ccc         rzsurf(nexs)  = rcmat(ncmat)-0.01
ccc         zrmin(nexs)   = zcmin(ncmat)
ccc         zrmax(nexs)   = zcmax(ncmat)
ccc*
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = vtxd2(1)+0.01
ccc         zcmin(ncmat)  = -vtxd2(3)
ccc         zcmax(ncmat)  = vtxd2(3)
ccc         xrlc(ncmat)   = (vtxd2(2)-vtxd2(1))/xrsili
ccc         xelosc(ncmat) = (vtxd2(2)-vtxd2(1))*dedxsili
ccc*
ccc         nexs          = nexs+1
ccc         typx(nexs)    = 'CYLI'
ccc         rzsurf(nexs)  = rcmat(ncmat)-0.01
ccc         zrmin(nexs)   = zcmin(ncmat)
ccc         zrmax(nexs)   = zcmax(ncmat)
ccc*
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = vtxd3(1)+0.01
ccc         zcmin(ncmat)  = -vtxd3(3)
ccc         zcmax(ncmat)  = vtxd3(3)
ccc         xrlc(ncmat)   = (vtxd3(2)-vtxd3(1))/xrsili
ccc         xelosc(ncmat) = (vtxd3(2)-vtxd3(1))*dedxsili
ccc*
ccc         nexs          = nexs+1
ccc         typx(nexs)    = 'CYLI'
ccc         rzsurf(nexs)  = rcmat(ncmat)-0.01
ccc         zrmin(nexs)   = zcmin(ncmat)
ccc         zrmax(nexs)   = zcmax(ncmat)
ccc*
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = vxtr(1)+0.01
ccc         zcmin(ncmat)  = -vxtr(3)
ccc         zcmax(ncmat)  = vxtr(3)
ccc         xrlc(ncmat)   = (vxtr(2)-vxtr(1))/xrsili
ccc         xelosc(ncmat) = (vxtr(2)-vxtr(1))*dedxsili
ccc*
ccc         nexs          = nexs+1
ccc         typx(nexs)    = 'CYLI'
ccc         rzsurf(nexs)  = rcmat(ncmat)-0.01
ccc         zrmin(nexs)   = zcmin(ncmat)
ccc         zrmax(nexs)   = zcmax(ncmat)
ccc*  endcaps (they are closer to the horizontal, so save them as cylinders)
ccc         do xs=-1.,1.,2.
ccc            ncmat         = ncmat+1
ccc            rcmat(ncmat)  = 0.5*(vxf1(2)+vxf1(4))
ccc            zcmin(ncmat)  = min(xs*(vxzf1+vxf1(1))*rcmat(ncmat)/vxf1(4),
ccc     $                          xs*(vxzf1-vxf1(1))*rcmat(ncmat)/vxf1(2))
ccc            zcmax(ncmat)  = max(xs*(vxzf1+vxf1(1))*rcmat(ncmat)/vxf1(4),
ccc     $                          xs*(vxzf1-vxf1(1))*rcmat(ncmat)/vxf1(2))
ccc* fudge factor is the longer path through the cone for 25 deg
ccc            xrlc(ncmat)   = bpct2/xrbery*sin(0.44)/sin(0.44+0.52)
ccc            xelosc(ncmat) = bpct2*dedxbery*sin(0.44)/sin(0.44+0.52)
ccc
ccc            xrlc(ncmat)   = vxthkf/xrsili
ccc            xelosc(ncmat) = vxthkf*dedxsili
ccc*
ccc            nexs          = nexs+1
ccc            typx(nexs)    = 'CYLI'
ccc            rzsurf(nexs)  = rcmat(ncmat)
ccc            zrmin(nexs)   = zcmin(ncmat)
ccc            zrmax(nexs)   = zcmax(ncmat)
ccc         enddo
ccc*
ccc         do xs=-1.,1.,2.
ccc            ncmat         = ncmat+1
ccc            rcmat(ncmat)  = 0.5*(vxf2(2)+vxf2(4))
ccc            zcmin(ncmat)  = min(xs*(vxzf2+vxf2(1))*rcmat(ncmat)/vxf2(4),
ccc     $                          xs*(vxzf2-vxf2(1))*rcmat(ncmat)/vxf2(2))
ccc            zcmax(ncmat)  = max(xs*(vxzf2+vxf2(1))*rcmat(ncmat)/vxf2(4),
ccc     $                          xs*(vxzf2-vxf2(1))*rcmat(ncmat)/vxf2(2))
ccc            xrlc(ncmat)   = vxthkf/xrsili
ccc            xelosc(ncmat) = vxthkf*dedxsili
ccc*
ccc            nexs          = nexs+1
ccc            typx(nexs)    = 'CYLI'
ccc            rzsurf(nexs)  = rcmat(ncmat)
ccc            zrmin(nexs)   = zcmin(ncmat)
ccc            zrmax(nexs)   = zcmax(ncmat)
ccc         enddo
ccc*
ccc      endif
ccc*
ccc      if (lvxdc) then
ccc         if (ccdtype.eq.1) then
ccc            ncmat         = ncmat+1
ccc            rcmat(ncmat)  = vxdrl1
ccc            zcmin(ncmat)  = -vxdzl1
ccc            zcmax(ncmat)  = vxdzl1
ccc            xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
ccc            xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
ccc*
ccc            nexs          = nexs+1
ccc            typx(nexs)    = 'CYLI'
ccc            rzsurf(nexs)  = rcmat(ncmat)
ccc            zrmin(nexs)   = zcmin(ncmat)
ccc            zrmax(nexs)   = zcmax(ncmat)
ccc         endif
ccc*
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = vxdrl2
ccc         zcmin(ncmat)  = -vxdzl2
ccc         zcmax(ncmat)  = vxdzl2
ccc         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
ccc         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
ccc*
ccc         nexs          = nexs+1
ccc         typx(nexs)    = 'CYLI'
ccc         rzsurf(nexs)  = rcmat(ncmat)
ccc         zrmin(nexs)   = zcmin(ncmat)
ccc         zrmax(nexs)   = zcmax(ncmat)
ccc*
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = vxdrl3
ccc         zcmin(ncmat)  = -vxdzl3
ccc         zcmax(ncmat)  = vxdzl3
ccc         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
ccc         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
ccc*
ccc         nexs          = nexs+1
ccc         typx(nexs)    = 'CYLI'
ccc         rzsurf(nexs)  = rcmat(ncmat)
ccc         zrmin(nexs)   = zcmin(ncmat)
ccc         zrmax(nexs)   = zcmax(ncmat)
ccc*
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = vxdrl4
ccc         zcmin(ncmat)  = -vxdzl4
ccc         zcmax(ncmat)  = vxdzl4
ccc         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
ccc         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
ccc*
ccc         nexs          = nexs+1
ccc         typx(nexs)    = 'CYLI'
ccc         rzsurf(nexs)  = rcmat(ncmat)
ccc         zrmin(nexs)   = zcmin(ncmat)
ccc         zrmax(nexs)   = zcmax(ncmat)
ccc*
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = vxdrl5
ccc         zcmin(ncmat)  = -vxdzl5
ccc         zcmax(ncmat)  = vxdzl5
ccc         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
ccc         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
ccc*
ccc         nexs          = nexs+1
ccc         typx(nexs)    = 'CYLI'
ccc         rzsurf(nexs)  = rcmat(ncmat)
ccc         zrmin(nexs)   = zcmin(ncmat)
ccc         zrmax(nexs)   = zcmax(ncmat)
ccc*     electronics
ccc         do xs=-1.,1.,2.
ccc            ncmat         = ncmat+1
ccc            rcmat(ncmat)  = VXDELR11
ccc            zcmin(ncmat)  = min(xs*VXDELZ11,xs*VXDELZ12)
ccc            zcmax(ncmat)  = max(xs*VXDELZ11,xs*VXDELZ12)
ccc            xrlc(ncmat)   = VXDELE/xrsili
ccc            xelosc(ncmat) = VXDELE*dedxsili
ccc*
ccc            ncmat         = ncmat+1
ccc            rcmat(ncmat)  = VXDELR21
ccc            zcmin(ncmat)  = min(xs*VXDELZ21,xs*VXDELZ22)
ccc            zcmax(ncmat)  = max(xs*VXDELZ21,xs*VXDELZ22)
ccc            xrlc(ncmat)   = VXDELE/xrsili
ccc            xelosc(ncmat) = VXDELE*dedxsili
ccc*
ccc            ncmat         = ncmat+1
ccc            rcmat(ncmat)  = VXDELR31
ccc            zcmin(ncmat)  = min(xs*VXDELZ31,xs*VXDELZ32)
ccc            zcmax(ncmat)  = max(xs*VXDELZ31,xs*VXDELZ32)
ccc            xrlc(ncmat)   = VXDELE/xrsili
ccc            xelosc(ncmat) = VXDELE*dedxsili
ccc*
ccc            ncmat         = ncmat+1
ccc            rcmat(ncmat)  = VXDELR41
ccc            zcmin(ncmat)  = min(xs*VXDELZ41,xs*VXDELZ42)
ccc            zcmax(ncmat)  = max(xs*VXDELZ41,xs*VXDELZ42)
ccc            xrlc(ncmat)   = VXDELE/xrsili
ccc            xelosc(ncmat) = VXDELE*dedxsili
ccc*
ccc            ncmat         = ncmat+1
ccc            rcmat(ncmat)  = VXDELR51
ccc            zcmin(ncmat)  = min(xs*VXDELZ51,xs*VXDELZ52)
ccc            zcmax(ncmat)  = max(xs*VXDELZ51,xs*VXDELZ52)
ccc            xrlc(ncmat)   = VXDELE/xrsili
ccc            xelosc(ncmat) = VXDELE*dedxsili
ccc         enddo
ccc*     support tube
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = vxdsupr
ccc         zcmin(ncmat)  = -vxdsupz
ccc         zcmax(ncmat)  = vxdsupz
ccc         xrlc(ncmat)   = vxdsupt/xrbery
ccc         xelosc(ncmat) = vxdsupt*dedxbery
ccc*     endplate
ccc*
ccc         do xs=-1.,1.,2.
ccc            npmat         = npmat+1
ccc            zpmat(npmat)  = xs*(VXDEPZ2+VXDEPZ1)/2.
ccc            rpmin(npmat)  = VXDEPR1
ccc            rpmax(npmat)  = VXDEPR2
ccc            xrlp(npmat)   = VXDEPT/xralu
ccc            xelosp(npmat) = VXDEPT*dedxalu
ccc         enddo
ccc*
ccc*     cryostat barrel
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = vxdcr2
ccc         zcmin(ncmat)  = -vxdcrz
ccc         zcmax(ncmat)  = vxdcrz
ccc         xrlc(ncmat)   = vxdcrt2/xralu
ccc         xelosc(ncmat) = vxdcrt2*dedxalu
ccc*     cryostate endplate
ccc*
ccc         do xs=-1.,1.,2.
ccc            npmat         = npmat+1
ccc            zpmat(npmat)  = xs*(VXDCRZ + VXDCRT1)
ccc            rpmin(npmat)  = VXDCR1
ccc            rpmax(npmat)  = VXDCR2+VXDCRT1+VXDCRT2
ccc            xrlp(npmat)   = VXDCRT2/xralu
ccc            xelosp(npmat) = VXDCRT2*dedxalu
ccc         enddo
ccc*
ccc      endif
ccc* ftd
ccc      if (lftd) then
ccc         do i=1,nftd
ccc            do xs=-1.,1.,2.
ccc               npmat         = npmat+1
ccc               zpmat(npmat)  = xs*ftdz(i)
ccc               rpmin(npmat)  = ftdrin(i)
ccc               rpmax(npmat)  = ftdrout(i)-ftdorth
ccc               if (i.le.nftdi) then
ccc                  xrlp(npmat)   = ftdthk(i)/xrsili *8.72/2.33
ccc                  xelosp(npmat) = ftdthk(i)*dedxsili *8.72/2.33
ccc               else
ccc                  xrlp(npmat)   = ftdthk(i)/xrsili
ccc                  xelosp(npmat) = ftdthk(i)*dedxsili
ccc               endif
ccc            enddo
ccc*
ccc            nexs          = nexs+1
ccc            typx(nexs)    = 'PLAN'
ccc            rzsurf(nexs)  = zpmat(npmat)
ccc            zrmin(nexs)   = rpmin(npmat)
ccc            zrmax(nexs)   = rpmax(npmat)
ccc         enddo
ccc*
ccc* outer support rings
ccc*
ccc         do i=nftdi+1,nftd
ccc            do xs=-1.,1.,2.
ccc               npmat         = npmat+1
ccc               zpmat(npmat)  = xs*ftdz(i)
ccc               rpmin(npmat)  = ftdrout(i)-ftdorth
ccc               rpmax(npmat)  = ftdrout(i)
ccc               xrlp(npmat)   = ftdorln/xrkapton
ccc               xelosp(npmat) = ftdorln*dedxkapton
ccc            enddo
ccc*
ccc         enddo
ccc*
ccc* outer cylinder and cables
ccc*
ccc         do xs=-1.,1.,2.
ccc            ncmat         = ncmat+1
ccc            rcmat(ncmat)  = ftdrout(nftd)+ftdorth+0.5*ftdspth
ccc            zcmin(ncmat)  = min(xs*ftdz(nftdi+2),xs*ftdz(nftd))
ccc            zcmax(ncmat)  = max(xs*ftdz(nftdi+2),xs*ftdz(nftd))
ccc            xrlc(ncmat)   = ftdspth/xrkapton+ftdcath/xrcopper
ccc            xelosc(ncmat) = ftdspth*dedxkapton+ftdcath*dedxcopper
ccc         enddo
ccc      endif
ccc* itc
ccc      if (litc) then
ccc*     inner tube
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = RITCINN
ccc         zcmin(ncmat)  = -ITCHLFZ
ccc         zcmax(ncmat)  = ITCHLFZ
ccc         xrlc(ncmat)   = ITCTHBI/xralu
ccc         xelosc(ncmat) = ITCTHBI*dedxalu
ccc*
ccc         nexs          = nexs+1
ccc         typx(nexs)    = 'CYLI'
ccc         rzsurf(nexs)  = rcmat(ncmat)+1
ccc         zrmin(nexs)   = zcmin(ncmat)
ccc         zrmax(nexs)   = zcmax(ncmat)
ccc*     gas
ccc         xstep = (ritcout-ritcinn)/real(nstitc)
ccc         do i=1,nstitc
ccc           ncmat         = ncmat+1
ccc           rcmat(ncmat)  = ritcinn+(real(i)-0.5)*xstep
ccc           zcmin(ncmat)  = -ITCHLFZ
ccc           zcmax(ncmat)  = ITCHLFZ
ccc           xrlc(ncmat)   = xstep/xrargon
ccc           xelosc(ncmat) = xstep*dedxargon
ccc        enddo
ccc*
ccc
ccc*     outer tube
ccc         ncmat         = ncmat+1
ccc         rcmat(ncmat)  = RITCOUT
ccc         zcmin(ncmat)  = -ITCHLFZ
ccc         zcmax(ncmat)  = ITCHLFZ
ccc         xrlc(ncmat)   = ITCTHBO/xralu
ccc         xelosc(ncmat) = ITCTHBO*dedxalu
ccc*     endplate
ccc         do xs=-1.,1.,2.
ccc            npmat         = npmat+1
ccc            zpmat(npmat)  = xs*(ITCHLFZ-ITCTHKE/2.)
ccc            rpmin(npmat)  = RITCINN
ccc            rpmax(npmat)  = RITCOUT
ccc            xrlp(npmat)   = 0.2
ccc            xelosp(npmat) = 0.2*xralu*dedxalu
ccc         enddo
ccc
ccc      endif
ccc* SIT
ccc      if (lsit) then
ccc         do i=1,nsit
ccc            ncmat         = ncmat+1
ccc            rcmat(ncmat)  = sitrad(i)
ccc            zcmin(ncmat)  = -sitlen(i)
ccc            zcmax(ncmat)  = sitlen(i)
ccc            xrlc(ncmat)   = sitthk/xrsili *8.72/2.33
ccc            xelosc(ncmat) = sitthk*dedxsili *8.72/2.33
ccc*
ccc            nexs          = nexs+1
ccc            typx(nexs)    = 'CYLI'
ccc            rzsurf(nexs)  = rcmat(ncmat)
ccc            zrmin(nexs)   = zcmin(ncmat)
ccc            zrmax(nexs)   = zcmax(ncmat)
ccc         enddo
ccc      endif
ccc

* TPC
c      if (ltpc) then
*     inner tube
         ncmat         = ncmat+1
c(KH)    rcmat(ncmat)  = RTPCINN
         rcmat(ncmat)  = RTPCINN+TPCTHBI/2.
         zcmin(ncmat)  = -TPCHLFZ
         zcmax(ncmat)  = TPCHLFZ
         xrlc(ncmat)   = TPCTHBI/xralu
         xelosc(ncmat) = TPCTHBI*dedxalu
*
         nexs          = nexs+1
         typx(nexs)    = 'CYLI'
         rzsurf(nexs)  = rcmat(ncmat)+1
         zrmin(nexs)   = zcmin(ncmat)
         zrmax(nexs)   = zcmax(ncmat)
*     gas
c(KH)    xstep = (rtpcout-rtpcinn)/real(nsttpc)
         xstep = (rtpcout-tpcthbo-tpcacri-tpcthbi)/real(nsttpc)
         do i=1,nsttpc
           ncmat         = ncmat+1
c(KH)      rcmat(ncmat)  = rtpcinn+(real(i)-0.5)*xstep
           rcmat(ncmat)  = rtpcinn+tpcthbi+(real(i)-0.5)*xstep
           zcmin(ncmat)  = -TPCHLFZ
           zcmax(ncmat)  = TPCHLFZ
           xrlc(ncmat)   = xstep/xrargon
           xelosc(ncmat) = xstep*dedxargon
        enddo
*

*     outer tube
         ncmat         = ncmat+1
         rcmat(ncmat)  = RTPCOUT
         zcmin(ncmat)  = -TPCHLFZ
         zcmax(ncmat)  = TPCHLFZ
         xrlc(ncmat)   = TPCTHBO/xralu
         xelosc(ncmat) = TPCTHBO*dedxalu
*     endplate
         do xs=-1.,1.,2.
            npmat         = npmat+1
            zpmat(npmat)  = xs*(TPCHLFZ-TPCTHKE/2.)
            rpmin(npmat)  = RTPCINN
            rpmax(npmat)  = RTPCOUT
            xrlp(npmat)   = 0.35
            xelosp(npmat) = 0.35*xralu*dedxalu
         enddo

c      endif

ccc* FCH
ccc      if (lfch) then
ccc         do xs=-1.,1.,2.
ccc            npmat         = npmat+1
ccc            zpmat(npmat)  = xs*(ZFCH-FCHTHK/2.)
ccc            rpmin(npmat)  = RTPCINN
ccc            rpmax(npmat)  = FCHROUT
ccc            xrlp(npmat)   = FCHTHK/xrstrw
ccc            xelosp(npmat) = FCHTHK*dedxstrw
ccc         enddo
ccc*
ccc         nexs          = nexs+1
ccc         nexfch        = nexs
ccc         typx(nexs)    = 'PLAN'
ccc         rzsurf(nexs)  = zpmat(npmat)
ccc         zrmin(nexs)   = rpmin(npmat)
ccc         zrmax(nexs)   = rpmax(npmat)
ccc      endif
ccc      IF ( printout .eq. 1 ) then
ccc         print*
ccc         print '('' material cylinders for track fit'')'
ccc         print '('' Nr    radius     z-min     z-max    X0     Eloss'')'
ccc         do i=1,ncmat
ccc            print '(i3,3f10.2,2f8.5)',
ccc     $           i,rcmat(i),zcmin(i),zcmax(i),xrlc(i),xelosc(i)
ccc         enddo
ccc         print*
ccc         print '('' material planes for track fit'')'
ccc         print '('' Nr       z       r-min     r-max    X0     Eloss'')'
ccc         do i=1,npmat
ccc            print '(i3,3f10.2,2f8.5)',
ccc     $           i,zpmat(i),rpmin(i),rpmax(i),xrlp(i),xelosp(i)
ccc         enddo
ccc      endif
*                                    chisquare probability cut values, |
*                                    indices = (forward, barrel):      |
*                                    total chisquare of the track fit  |
      PRTOT(0) = 0.0001
      PRTOT(1) = 0.0001
*                                    smoothed chisquare (outlier test) |
      PRSMO(0) = 0.0001
      PRSMO(1) = 0.0001
*
  999 return
      end


