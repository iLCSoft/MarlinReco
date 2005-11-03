*  
* $Id: setmat.F,v 1.3 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2005/08/02 10:21:38  aplin
* placed initialisation of BFIELD here and corrected debug output
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
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
#include "include/tkffread.inc"
#include "include/coildims.inc"
#include "include/allpipe.inc"
#include "include/maskdims.inc"
c#include "turnon.inc"
#include "include/pixpipe.inc"
#include "include/vxdpix.inc"
#include "include/ccdpipe.inc"
#include "include/vxdccd.inc"
c#include "itcdims.inc"
#include "include/tpcdims.inc"
c#include "ftddims.inc"
#include "include/sitdims.inc"
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

CSJA FIXME: include here as a temporary solution all the geometric information for the VTX and beampipe.

      MKTHMAX=83.1

C all dimensions are in cm

C Thickness of different parts
      BPT0 = 0.05
      BPT1 = 0.1
      BPT2 = 0.1
      BPT3 = 0.1
      BPT4 = 0.1
      BPT5 = 0.1
C Outer length:
      BPPZ=800.

C Beampipe dimensions for PIXEL vertex detector
C (Based on code supplied by M. Battaglia):

C Central part:
c Inner Radius:
      BPP0(1)=1.4
c Outer Radius
      BPP0(2)=1.4 + BPT0
c half length
      BPP0(3)=12.5

C Conical parts:
c Half length
      BPP1(1)= 31.25
c first inner radius
      BPP1(2)=BPP0(1)
c first outer radius
      BPP1(3)=BPP1(2) + BPT1
c Second inner radius
      BPP1(4)=3.5
c Second outer radius
      BPP1(5)=BPP1(4) + BPT1

C----------------------------------------------------------------------
C Beampipe dimensions for CCD vertex detector:
C Radiation length of beryllium is taken to be 35.28 cm
c Radiation length of steel is taken to be 1.76 cm
c for CCD option beampipe has the structure, working from the ip outwards:
c Be cylinder, Be cone, Be cylinder, steel cylinder
c thicknesses :
C !!!! NEEDS CHECKING !!!!

*        BPCT1=1.4E-3 *35.28
*        BPCT2=2.8E-3 *35.28
*        BPCT3=2.8E-3 *35.28
*        BPCT4=5.7E-3 *1.76
      BPCT1 = BPT0
      BPCT2 = BPT1
      BPCT3 = BPT2
      BPCT4 = BPT3

C BP inner radii:

c      IF ( CCDTYPE .EQ. 1 ) THEN
        BPCINR1=1.4
        BPCINR2=2.2
c      ELSE IF ( CCDTYPE .EQ. 2 ) THEN
c        BPCINR1 = 2.0
c        BPCINR2 = 2.0
c      END IF

      BPCINR3=2.0
      BPCINR3 = 2.2

c      IF (ILAT.EQ.1) THEN
c        BPCINZ4 = LTZ1-3.5
c        BPCINR4 = BPCINZ4 * TAN(MKTHMAX/1000.) - BPT4
c        BPCINR4 = BPCINZ4 * TAN(MKTHMAX/1000.)
c      END IF

C Z Positions of CHANGES in radius
      BPCINZ1=3.1
      BPCINZ2=15.0
      BPCINZ3=15.2

CCC       BPCINZ1 = 3.7
      BPCINZ1 = 6.1
CCC      BPCINZ2 = 4.7
C*kmo      BPCINZ2 = 7.2
      BPCINZ2 = 17.2
      BPCINZ3 = 17.2

c      BPCINR2 = BPCINZ2* TAN(MKTHMAX/1000.)
      BPCINZ3 = BPCINR3 / TAN(MKTHMAX/1000.)

C     from here on CCD and APS beam pipe are the same.


csja call vxdgeom here to set up the vtx geometry so that material definitions for vtx can be set
      call VXDGEOM

CCC FIXME: SET BFIELD HERE
      BFIELD = 40.0

      consb = bfield*fct

      ncmat = 0
      npmat = 0
      nexs  = 0
* beam pipe
csja      if (lvxdp) then
csja         ncmat         = ncmat+1
csja         rcmat(ncmat)  = bpp0(1)
csja         zcmin(ncmat)  = -bpp0(3)
csja         zcmax(ncmat)  = bpp0(3)
csja         xrlc(ncmat)   = (bpp0(2)-bpp0(1))/xrbery
csja         xelosc(ncmat) = (bpp0(2)-bpp0(1))*dedxbery
csja*
csja*          radius of 1st conical part is inner radius
csja         ncmat         = ncmat+1
csja         rcmat(ncmat)  = bpp1(2)
csja         zcmin(ncmat)  = bpp0(3)
csja         zcmax(ncmat)  = bpp0(3)+2.*bpp1(1)
csja* fudge factor is the longer path through the cone for 8 deg
csja         xrlc(ncmat)   = (bpp1(3)-bpp1(2))/xrbery*
csja     $        (sin(0.14)/sin(0.14-atan(0.5*(bpp1(4)-bpp1(2))/bpp1(1))))
csja         xelosc(ncmat) = (bpp1(3)-bpp1(2))*dedxbery*
csja     $        (sin(0.14)/sin(0.14-atan(0.5*(bpp1(4)-bpp1(2))/bpp1(1))))
csja*
csja         ncmat         = ncmat+1
csja         rcmat(ncmat)  = rcmat(ncmat-1)
csja         zcmin(ncmat)  = -zcmax(ncmat-1)
csja         zcmax(ncmat)  = -zcmin(ncmat-1)
csja         xrlc(ncmat)   = xrlc(ncmat-1)
csja         xelosc(ncmat) = xelosc(ncmat-1)
csja*
csja      endif
csja*

C FIXME: the CCD is set to default for now
csja      if (lvxdc) then
         ncmat         = ncmat+1
         rcmat(ncmat)  = BPCINR1
         zcmin(ncmat)  = -BPCINZ1
         zcmax(ncmat)  = BPCINZ1
         xrlc(ncmat)   = BPCT1/xrbery
         xelosc(ncmat) = BPCT1*dedxbery
*
*          radius of 1st conical part is middle of cone
         ncmat         = ncmat+1
         rcmat(ncmat)  = 0.5*(BPCINR1+BPCINR2)
         zcmin(ncmat)  = bpcinz1*rcmat(ncmat)/bpcinr1
c         zcmax(ncmat)  = bpcinz3*rcmat(ncmat)/bpcinr2
         zcmax(ncmat)  = 100.
* fudge factor is the longer path through the cone for 10 deg
         xrlc(ncmat)   = bpct2/xrbery*
     $        (sin(0.175)/sin(0.175-atan((BPCINR2-BPCINR1)/
     $        (BPCINZ2-BPCINZ1))))
         xelosc(ncmat) = bpct2*dedxbery*
     $        (sin(0.175)/sin(0.175-atan((BPCINR2-BPCINR1)/
     $        (BPCINZ2-BPCINZ1))))
*
         ncmat         = ncmat+1
         rcmat(ncmat)  = rcmat(ncmat-1)
         zcmin(ncmat)  = -zcmax(ncmat-1)
         zcmax(ncmat)  = -zcmin(ncmat-1)
         xrlc(ncmat)   = xrlc(ncmat-1)
         xelosc(ncmat) = xelosc(ncmat-1)
*

csja      endif
*

      nexs          = nexs+1
      typx(nexs)    = 'CYLI'
C FIXME:SJA: this line is commened out and the value replace with 0 to ensure that the material
C     description of brahm and marlin "agree"
      rzsurf(nexs)  = rcmat(1)-0.01
c      rzsurf(nexs)  = 0.
      zrmin(nexs)   = -1000.
      zrmax(nexs)   = 1000.


ccc*  microvertex
csja      if (lvxdp) then
csja         ncmat         = ncmat+1
csja         rcmat(ncmat)  = vtxd1(1)+0.01
csja         zcmin(ncmat)  = -vtxd1(3)
csja         zcmax(ncmat)  = vtxd1(3)
csja         xrlc(ncmat)   = (vtxd1(2)-vtxd1(1))/xrsili
csja         xelosc(ncmat) = (vtxd1(2)-vtxd1(1))*dedxsili
csja*
csja         nexs          = nexs+1
csja         typx(nexs)    = 'CYLI'
csja         rzsurf(nexs)  = rcmat(ncmat)-0.01
csja         zrmin(nexs)   = zcmin(ncmat)
csja         zrmax(nexs)   = zcmax(ncmat)
csja*
csja         ncmat         = ncmat+1
csja         rcmat(ncmat)  = vtxd2(1)+0.01
csja         zcmin(ncmat)  = -vtxd2(3)
csja         zcmax(ncmat)  = vtxd2(3)
csja         xrlc(ncmat)   = (vtxd2(2)-vtxd2(1))/xrsili
csja         xelosc(ncmat) = (vtxd2(2)-vtxd2(1))*dedxsili
csja*
csja         nexs          = nexs+1
csja         typx(nexs)    = 'CYLI'
csja         rzsurf(nexs)  = rcmat(ncmat)-0.01
csja         zrmin(nexs)   = zcmin(ncmat)
csja         zrmax(nexs)   = zcmax(ncmat)
csja*
csja         ncmat         = ncmat+1
csja         rcmat(ncmat)  = vtxd3(1)+0.01
csja         zcmin(ncmat)  = -vtxd3(3)
csja         zcmax(ncmat)  = vtxd3(3)
csja         xrlc(ncmat)   = (vtxd3(2)-vtxd3(1))/xrsili
csja         xelosc(ncmat) = (vtxd3(2)-vtxd3(1))*dedxsili
csja*
csja         nexs          = nexs+1
csja         typx(nexs)    = 'CYLI'
csja         rzsurf(nexs)  = rcmat(ncmat)-0.01
csja         zrmin(nexs)   = zcmin(ncmat)
csja         zrmax(nexs)   = zcmax(ncmat)
csja*
csja         ncmat         = ncmat+1
csja         rcmat(ncmat)  = vxtr(1)+0.01
csja         zcmin(ncmat)  = -vxtr(3)
csja         zcmax(ncmat)  = vxtr(3)
csja         xrlc(ncmat)   = (vxtr(2)-vxtr(1))/xrsili
csja         xelosc(ncmat) = (vxtr(2)-vxtr(1))*dedxsili
csja*
csja         nexs          = nexs+1
csja         typx(nexs)    = 'CYLI'
csja         rzsurf(nexs)  = rcmat(ncmat)-0.01
csja         zrmin(nexs)   = zcmin(ncmat)
csja         zrmax(nexs)   = zcmax(ncmat)
csja*  endcaps (they are closer to the horizontal, so save them as cylinders)
csja         do xs=-1.,1.,2.
csja            ncmat         = ncmat+1
csja            rcmat(ncmat)  = 0.5*(vxf1(2)+vxf1(4))
csja            zcmin(ncmat)  = min(xs*(vxzf1+vxf1(1))*rcmat(ncmat)/vxf1(4),
csja     $                          xs*(vxzf1-vxf1(1))*rcmat(ncmat)/vxf1(2))
csja            zcmax(ncmat)  = max(xs*(vxzf1+vxf1(1))*rcmat(ncmat)/vxf1(4),
csja     $                          xs*(vxzf1-vxf1(1))*rcmat(ncmat)/vxf1(2))
csja* fudge factor is the longer path through the cone for 25 deg
csja            xrlc(ncmat)   = bpct2/xrbery*sin(0.44)/sin(0.44+0.52)
csja            xelosc(ncmat) = bpct2*dedxbery*sin(0.44)/sin(0.44+0.52)
csja
csja            xrlc(ncmat)   = vxthkf/xrsili
csja            xelosc(ncmat) = vxthkf*dedxsili
csja*
csja            nexs          = nexs+1
csja            typx(nexs)    = 'CYLI'
csja            rzsurf(nexs)  = rcmat(ncmat)
csja            zrmin(nexs)   = zcmin(ncmat)
csja            zrmax(nexs)   = zcmax(ncmat)
csja         enddo
csja*
csja         do xs=-1.,1.,2.
csja            ncmat         = ncmat+1
csja            rcmat(ncmat)  = 0.5*(vxf2(2)+vxf2(4))
csja            zcmin(ncmat)  = min(xs*(vxzf2+vxf2(1))*rcmat(ncmat)/vxf2(4),
csja     $                          xs*(vxzf2-vxf2(1))*rcmat(ncmat)/vxf2(2))
csja            zcmax(ncmat)  = max(xs*(vxzf2+vxf2(1))*rcmat(ncmat)/vxf2(4),
csja     $                          xs*(vxzf2-vxf2(1))*rcmat(ncmat)/vxf2(2))
csja            xrlc(ncmat)   = vxthkf/xrsili
csja            xelosc(ncmat) = vxthkf*dedxsili
csja*
csja            nexs          = nexs+1
csja            typx(nexs)    = 'CYLI'
csja            rzsurf(nexs)  = rcmat(ncmat)
csja            zrmin(nexs)   = zcmin(ncmat)
csja            zrmax(nexs)   = zcmax(ncmat)
csja         enddo
csja*
csja      endif
csja*


c FIXME: CCD set to default for now
csja      if (lvxdc) then
csja         if (ccdtype.eq.1) then
            ncmat         = ncmat+1
            rcmat(ncmat)  = vxdrl1
            zcmin(ncmat)  = -vxdzl1
            zcmax(ncmat)  = vxdzl1
            xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
            xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
*
            nexs          = nexs+1
            typx(nexs)    = 'CYLI'
            rzsurf(nexs)  = rcmat(ncmat)
            zrmin(nexs)   = zcmin(ncmat)
            zrmax(nexs)   = zcmax(ncmat)
csja         endif
*
         ncmat         = ncmat+1
         rcmat(ncmat)  = vxdrl2
         zcmin(ncmat)  = -vxdzl2
         zcmax(ncmat)  = vxdzl2
         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
*
         nexs          = nexs+1
         typx(nexs)    = 'CYLI'
         rzsurf(nexs)  = rcmat(ncmat)
         zrmin(nexs)   = zcmin(ncmat)
         zrmax(nexs)   = zcmax(ncmat)
*
         ncmat         = ncmat+1
         rcmat(ncmat)  = vxdrl3
         zcmin(ncmat)  = -vxdzl3
         zcmax(ncmat)  = vxdzl3
         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
*
         nexs          = nexs+1
         typx(nexs)    = 'CYLI'
         rzsurf(nexs)  = rcmat(ncmat)
         zrmin(nexs)   = zcmin(ncmat)
         zrmax(nexs)   = zcmax(ncmat)
*
         ncmat         = ncmat+1
         rcmat(ncmat)  = vxdrl4
         zcmin(ncmat)  = -vxdzl4
         zcmax(ncmat)  = vxdzl4
         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
*
         nexs          = nexs+1
         typx(nexs)    = 'CYLI'
         rzsurf(nexs)  = rcmat(ncmat)
         zrmin(nexs)   = zcmin(ncmat)
         zrmax(nexs)   = zcmax(ncmat)
*
         ncmat         = ncmat+1
         rcmat(ncmat)  = vxdrl5
         zcmin(ncmat)  = -vxdzl5
         zcmax(ncmat)  = vxdzl5
         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
*
         nexs          = nexs+1
         typx(nexs)    = 'CYLI'
         rzsurf(nexs)  = rcmat(ncmat)
         zrmin(nexs)   = zcmin(ncmat)
         zrmax(nexs)   = zcmax(ncmat)
*     electronics
         do xs=-1.,1.,2.
            ncmat         = ncmat+1
            rcmat(ncmat)  = VXDELR11
            zcmin(ncmat)  = min(xs*VXDELZ11,xs*VXDELZ12)
            zcmax(ncmat)  = max(xs*VXDELZ11,xs*VXDELZ12)
            xrlc(ncmat)   = VXDELE/xrsili
            xelosc(ncmat) = VXDELE*dedxsili
*
            ncmat         = ncmat+1
            rcmat(ncmat)  = VXDELR21
            zcmin(ncmat)  = min(xs*VXDELZ21,xs*VXDELZ22)
            zcmax(ncmat)  = max(xs*VXDELZ21,xs*VXDELZ22)
            xrlc(ncmat)   = VXDELE/xrsili
            xelosc(ncmat) = VXDELE*dedxsili
*
            ncmat         = ncmat+1
            rcmat(ncmat)  = VXDELR31
            zcmin(ncmat)  = min(xs*VXDELZ31,xs*VXDELZ32)
            zcmax(ncmat)  = max(xs*VXDELZ31,xs*VXDELZ32)
            xrlc(ncmat)   = VXDELE/xrsili
            xelosc(ncmat) = VXDELE*dedxsili
*
            ncmat         = ncmat+1
            rcmat(ncmat)  = VXDELR41
            zcmin(ncmat)  = min(xs*VXDELZ41,xs*VXDELZ42)
            zcmax(ncmat)  = max(xs*VXDELZ41,xs*VXDELZ42)
            xrlc(ncmat)   = VXDELE/xrsili
            xelosc(ncmat) = VXDELE*dedxsili
*
            ncmat         = ncmat+1
            rcmat(ncmat)  = VXDELR51
            zcmin(ncmat)  = min(xs*VXDELZ51,xs*VXDELZ52)
            zcmax(ncmat)  = max(xs*VXDELZ51,xs*VXDELZ52)
            xrlc(ncmat)   = VXDELE/xrsili
            xelosc(ncmat) = VXDELE*dedxsili
         enddo
*     support tube
         ncmat         = ncmat+1
         rcmat(ncmat)  = vxdsupr
         zcmin(ncmat)  = -vxdsupz
         zcmax(ncmat)  = vxdsupz
         xrlc(ncmat)   = vxdsupt/xrbery
         xelosc(ncmat) = vxdsupt*dedxbery
*     endplate
*
         do xs=-1.,1.,2.
            npmat         = npmat+1
            zpmat(npmat)  = xs*(VXDEPZ2+VXDEPZ1)/2.
            rpmin(npmat)  = VXDEPR1
            rpmax(npmat)  = VXDEPR2
            xrlp(npmat)   = VXDEPT/xralu
            xelosp(npmat) = VXDEPT*dedxalu
         enddo
*
*     cryostat barrel
         ncmat         = ncmat+1
         rcmat(ncmat)  = vxdcr2
         zcmin(ncmat)  = -vxdcrz
         zcmax(ncmat)  = vxdcrz
         xrlc(ncmat)   = vxdcrt2/xralu
         xelosc(ncmat) = vxdcrt2*dedxalu
*     cryostate endplate
*
         do xs=-1.,1.,2.
            npmat         = npmat+1
            zpmat(npmat)  = xs*(VXDCRZ + VXDCRT1)
            rpmin(npmat)  = VXDCR1
            rpmax(npmat)  = VXDCR2+VXDCRT1+VXDCRT2
            xrlp(npmat)   = VXDCRT2/xralu
            xelosp(npmat) = VXDCRT2*dedxalu
         enddo

*
csja      endif


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
* SIT

csja FIXME: the SIT geom will be called from here for now.
         call sitgeom

csja      if (lsit) then
         do i=1,nsit
            ncmat         = ncmat+1
            rcmat(ncmat)  = sitrad(i)
            zcmin(ncmat)  = -sitlen(i)
            zcmax(ncmat)  = sitlen(i)
            xrlc(ncmat)   = sitthk/xrsili *8.72/2.33
            xelosc(ncmat) = sitthk*dedxsili *8.72/2.33
*
            nexs          = nexs+1
            typx(nexs)    = 'CYLI'
            rzsurf(nexs)  = rcmat(ncmat)
            zrmin(nexs)   = zcmin(ncmat)
            zrmax(nexs)   = zcmax(ncmat)
         enddo
csja      endif


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

         IF (TKSTDBG.GT.0) then
            print*
            print '('' material cylinders for track fit'')'
            print '('' Nr    radius     z-min     z-max    X0     Eloss'')'
            do i=1,ncmat
               print '(i3,3f10.2,2f8.5)',
     $              i,rcmat(i),zcmin(i),zcmax(i),xrlc(i),xelosc(i)
            enddo
            print*
            print '('' material planes for track fit'')'
            print '('' Nr       z       r-min     r-max    X0     Eloss'')'
            do i=1,npmat
               print '(i3,3f10.2,2f8.5)',
     $              i,zpmat(i),rpmin(i),rpmax(i),xrlp(i),xelosp(i)
            enddo
         endif
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


