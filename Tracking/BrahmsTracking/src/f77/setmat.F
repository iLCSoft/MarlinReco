*  
* $Id: setmat.F,v 1.7 2007-09-05 09:47:29 rasp Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.5  2006/11/02 11:49:00  rasp
*
* Modified files adopted to FullLDCTracking processor
*
* Revision 1.4  2006/06/28 15:29:04  aplin
* The B-Field is now variable for LEPTracking via the gear xml file. The B-Field is specified in the TPCParameters as follows: <parameter name="tpcBField" type="double"> 4.0  </parameter>
*
* The value is passed internaly to the F77 code via the same function which passes the TPC geometry i.e. gettpcgeom(float* innerrad, float* outerrad, int* npadrows, float* maxdrift, float* tpcpixz, float* ionpoten, float* tpcrpres, float* tpczres, float* tpcbfield). It is set in setmat.F. tpcgeom.F had to be modified as it also uses gettpcgeom, although it does not make use of the B-Field.
*
* Revision 1.3  2005/11/03 15:16:14  aplin
* Added the Trackstring creation and the biulding of full Track candiates (TK's) which have passed the Delphi Ambiguity resolver fxambi. The material description of the vtx detector, as for the TPC, is hard coded in setmat. Presently the VTX and SIT resolutions are hard coded in LEPTrackingProcessor. The debug output has been reduced and can be controlled via TKSTDBG etc. in tkinit.F. delsolve contains the delphi ambuguity resolver written in C and is contained in the directory named C. The Tk's are written back into the C++ side in tktrev. The corresponding Tk bank structure analogous to the TE bank structure has been added in tktkbank whilst the access wrapper functions are contained in LEPTracking.
*
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

      REAL INNERRAD, OUTERRAD, MAXDRIFT, TPCPIXZ, IONPOTEN 
      REAL TPCRPRES, TPCZRES, TPCBFIELD
      
      INTEGER NPADROWS

      INTEGER TRY
      INTEGER gettpcgeom
      external gettpcgeom

      LOGICAL LFIRST
      
      DATA LFIRST /.TRUE./

      SAVE LFIRST


C CALL wrapper function to get tpc parameters from GEAR including the BFIELD
      TRY = gettpcgeom(INNERRAD, OUTERRAD, NPADROWS, MAXDRIFT, TPCPIXZ, 
     +     IONPOTEN, TPCRPRES, TPCZRES, TPCBFIELD)

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
      BFIELD =  TPCBFIELD*10.0
c      BFIELD = 40.0
      write(*,*) "BFIELD = ", BFIELD/10.0, " Tesla"

      consb = bfield*fct

c      ncmat = 0
c      npmat = 0
c      nexs  = 0
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
c$$$         ncmat         = ncmat+1
c$$$         rcmat(ncmat)  = BPCINR1
c$$$         zcmin(ncmat)  = -BPCINZ1
c$$$         zcmax(ncmat)  = BPCINZ1
c$$$         xrlc(ncmat)   = BPCT1/xrbery
c$$$         xelosc(ncmat) = BPCT1*dedxbery
*
*          radius of 1st conical part is middle of cone
c$$$         ncmat         = ncmat+1
c$$$         rcmat(ncmat)  = 0.5*(BPCINR1+BPCINR2)
c$$$         zcmin(ncmat)  = bpcinz1*rcmat(ncmat)/bpcinr1
c$$$c         zcmax(ncmat)  = bpcinz3*rcmat(ncmat)/bpcinr2
c$$$         zcmax(ncmat)  = 100.
c$$$* fudge factor is the longer path through the cone for 10 deg
c$$$         xrlc(ncmat)   = bpct2/xrbery*
c$$$     $        (sin(0.175)/sin(0.175-atan((BPCINR2-BPCINR1)/
c$$$     $        (BPCINZ2-BPCINZ1))))
c$$$         xelosc(ncmat) = bpct2*dedxbery*
c$$$     $        (sin(0.175)/sin(0.175-atan((BPCINR2-BPCINR1)/
c$$$     $        (BPCINZ2-BPCINZ1))))
c$$$*
c$$$         ncmat         = ncmat+1
c$$$         rcmat(ncmat)  = rcmat(ncmat-1)
c$$$         zcmin(ncmat)  = -zcmax(ncmat-1)
c$$$         zcmax(ncmat)  = -zcmin(ncmat-1)
c$$$         xrlc(ncmat)   = xrlc(ncmat-1)
c$$$         xelosc(ncmat) = xelosc(ncmat-1)
*

csja      endif
*

c$$$      nexs          = nexs+1
c$$$      typx(nexs)    = 'CYLI'
C FIXME:SJA: this line is commened out and the value replace with 0 to ensure that the material
C     description of brahm and marlin "agree"
c      rzsurf(nexs)  = rcmat(1)-0.01
c      rzsurf(nexs)  = 0.
c      zrmin(nexs)   = -1000.
c      zrmax(nexs)   = 1000.


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
c$$$            ncmat         = ncmat+1
c$$$            rcmat(ncmat)  = vxdrl1
c$$$            zcmin(ncmat)  = -vxdzl1
c$$$            zcmax(ncmat)  = vxdzl1
c$$$            xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
c$$$            xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
c$$$*
c$$$            nexs          = nexs+1
c$$$            typx(nexs)    = 'CYLI'
c$$$            rzsurf(nexs)  = rcmat(ncmat)
c$$$            zrmin(nexs)   = zcmin(ncmat)
c$$$            zrmax(nexs)   = zcmax(ncmat)
c$$$csja         endif
c$$$*
c$$$         ncmat         = ncmat+1
c$$$         rcmat(ncmat)  = vxdrl2
c$$$         zcmin(ncmat)  = -vxdzl2
c$$$         zcmax(ncmat)  = vxdzl2
c$$$         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
c$$$         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
c$$$*
c$$$         nexs          = nexs+1
c$$$         typx(nexs)    = 'CYLI'
c$$$         rzsurf(nexs)  = rcmat(ncmat)
c$$$         zrmin(nexs)   = zcmin(ncmat)
c$$$         zrmax(nexs)   = zcmax(ncmat)
c$$$*
c$$$         ncmat         = ncmat+1
c$$$         rcmat(ncmat)  = vxdrl3
c$$$         zcmin(ncmat)  = -vxdzl3
c$$$         zcmax(ncmat)  = vxdzl3
c$$$         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
c$$$         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
c$$$*
c$$$         nexs          = nexs+1
c$$$         typx(nexs)    = 'CYLI'
c$$$         rzsurf(nexs)  = rcmat(ncmat)
c$$$         zrmin(nexs)   = zcmin(ncmat)
c$$$         zrmax(nexs)   = zcmax(ncmat)
c$$$*
c$$$         ncmat         = ncmat+1
c$$$         rcmat(ncmat)  = vxdrl4
c$$$         zcmin(ncmat)  = -vxdzl4
c$$$         zcmax(ncmat)  = vxdzl4
c$$$         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
c$$$         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
c$$$*
c$$$         nexs          = nexs+1
c$$$         typx(nexs)    = 'CYLI'
c$$$         rzsurf(nexs)  = rcmat(ncmat)
c$$$         zrmin(nexs)   = zcmin(ncmat)
c$$$         zrmax(nexs)   = zcmax(ncmat)
c$$$*
c$$$         ncmat         = ncmat+1
c$$$         rcmat(ncmat)  = vxdrl5
c$$$         zcmin(ncmat)  = -vxdzl5
c$$$         zcmax(ncmat)  = vxdzl5
c$$$         xrlc(ncmat)   = vxdclt/xrsili+vxdsup/xrbery
c$$$         xelosc(ncmat) = vxdclt*dedxsili+vxdsup*dedxbery
c$$$*
c$$$         nexs          = nexs+1
c$$$         typx(nexs)    = 'CYLI'
c$$$         rzsurf(nexs)  = rcmat(ncmat)
c$$$         zrmin(nexs)   = zcmin(ncmat)
c$$$         zrmax(nexs)   = zcmax(ncmat)
*     electronics
c$$$         do xs=-1.,1.,2.
c$$$            ncmat         = ncmat+1
c$$$            rcmat(ncmat)  = VXDELR11
c$$$            zcmin(ncmat)  = min(xs*VXDELZ11,xs*VXDELZ12)
c$$$            zcmax(ncmat)  = max(xs*VXDELZ11,xs*VXDELZ12)
c$$$            xrlc(ncmat)   = VXDELE/xrsili
c$$$            xelosc(ncmat) = VXDELE*dedxsili
c$$$*
c$$$            ncmat         = ncmat+1
c$$$            rcmat(ncmat)  = VXDELR21
c$$$            zcmin(ncmat)  = min(xs*VXDELZ21,xs*VXDELZ22)
c$$$            zcmax(ncmat)  = max(xs*VXDELZ21,xs*VXDELZ22)
c$$$            xrlc(ncmat)   = VXDELE/xrsili
c$$$            xelosc(ncmat) = VXDELE*dedxsili
c$$$*
c$$$            ncmat         = ncmat+1
c$$$            rcmat(ncmat)  = VXDELR31
c$$$            zcmin(ncmat)  = min(xs*VXDELZ31,xs*VXDELZ32)
c$$$            zcmax(ncmat)  = max(xs*VXDELZ31,xs*VXDELZ32)
c$$$            xrlc(ncmat)   = VXDELE/xrsili
c$$$            xelosc(ncmat) = VXDELE*dedxsili
c$$$*
c$$$            ncmat         = ncmat+1
c$$$            rcmat(ncmat)  = VXDELR41
c$$$            zcmin(ncmat)  = min(xs*VXDELZ41,xs*VXDELZ42)
c$$$            zcmax(ncmat)  = max(xs*VXDELZ41,xs*VXDELZ42)
c$$$            xrlc(ncmat)   = VXDELE/xrsili
c$$$            xelosc(ncmat) = VXDELE*dedxsili
c$$$*
c$$$            ncmat         = ncmat+1
c$$$            rcmat(ncmat)  = VXDELR51
c$$$            zcmin(ncmat)  = min(xs*VXDELZ51,xs*VXDELZ52)
c$$$            zcmax(ncmat)  = max(xs*VXDELZ51,xs*VXDELZ52)
c$$$            xrlc(ncmat)   = VXDELE/xrsili
c$$$            xelosc(ncmat) = VXDELE*dedxsili
c$$$         enddo
*     support tube
c$$$         ncmat         = ncmat+1
c$$$         rcmat(ncmat)  = vxdsupr
c$$$         zcmin(ncmat)  = -vxdsupz
c$$$         zcmax(ncmat)  = vxdsupz
c$$$         xrlc(ncmat)   = vxdsupt/xrbery
c$$$         xelosc(ncmat) = vxdsupt*dedxbery
c$$$*     endplate
c$$$*
c$$$         do xs=-1.,1.,2.
c$$$            npmat         = npmat+1
c$$$            zpmat(npmat)  = xs*(VXDEPZ2+VXDEPZ1)/2.
c$$$            rpmin(npmat)  = VXDEPR1
c$$$            rpmax(npmat)  = VXDEPR2
c$$$            xrlp(npmat)   = VXDEPT/xralu
c$$$            xelosp(npmat) = VXDEPT*dedxalu
c$$$         enddo
*
*     cryostat barrel
c$$$         ncmat         = ncmat+1
c$$$         rcmat(ncmat)  = vxdcr2
c$$$         zcmin(ncmat)  = -vxdcrz
c$$$         zcmax(ncmat)  = vxdcrz
c$$$         xrlc(ncmat)   = vxdcrt2/xralu
c$$$         xelosc(ncmat) = vxdcrt2*dedxalu
*     cryostate endplate
*
c$$$         do xs=-1.,1.,2.
c$$$            npmat         = npmat+1
c$$$            zpmat(npmat)  = xs*(VXDCRZ + VXDCRT1)
c$$$            rpmin(npmat)  = VXDCR1
c$$$            rpmax(npmat)  = VXDCR2+VXDCRT1+VXDCRT2
c$$$            xrlp(npmat)   = VXDCRT2/xralu
c$$$            xelosp(npmat) = VXDCRT2*dedxalu
c$$$         enddo

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
c$$$
c$$$csja      if (lsit) then
c$$$         do i=1,nsit
c$$$            ncmat         = ncmat+1
c$$$            rcmat(ncmat)  = sitrad(i)
c$$$            zcmin(ncmat)  = -sitlen(i)
c$$$            zcmax(ncmat)  = sitlen(i)
c$$$            xrlc(ncmat)   = sitthk/xrsili *8.72/2.33
c$$$            xelosc(ncmat) = sitthk*dedxsili *8.72/2.33
c$$$*
c$$$            nexs          = nexs+1
c$$$            typx(nexs)    = 'CYLI'
c$$$            rzsurf(nexs)  = rcmat(ncmat)
c$$$            zrmin(nexs)   = zcmin(ncmat)
c$$$            zrmax(nexs)   = zcmax(ncmat)
c$$$         enddo
csja      endif


* TPC
c      if (ltpc) then
*     inner tube
c$$$         ncmat         = ncmat+1
c$$$c(KH)    rcmat(ncmat)  = RTPCINN
c$$$         rcmat(ncmat)  = RTPCINN+TPCTHBI/2.
c$$$         zcmin(ncmat)  = -TPCHLFZ
c$$$         zcmax(ncmat)  = TPCHLFZ
c$$$         xrlc(ncmat)   = TPCTHBI/xralu
c$$$         xelosc(ncmat) = TPCTHBI*dedxalu
c$$$*
c$$$         nexs          = nexs+1
c$$$         typx(nexs)    = 'CYLI'
c$$$         rzsurf(nexs)  = rcmat(ncmat)+1
c$$$         zrmin(nexs)   = zcmin(ncmat)
c$$$         zrmax(nexs)   = zcmax(ncmat)
c$$$
c$$$         WRITE(*,*) 'setmat.F ---->'
c$$$
c$$$         write(*,*) 'RTPCINN = ',RTPCINN
c$$$         WRITE(*,*) 'TPCHLFZ = ',TPCHLFZ
c$$$         WRITE(*,*) 'TPCTHBI = ',TPCTHBI


*     gas
c(KH)    xstep = (rtpcout-rtpcinn)/real(nsttpc)
c$$$         xstep = (rtpcout-tpcthbo-tpcacri-tpcthbi)/real(nsttpc)
c$$$         WRITE(*,*) 'nsttpc = ',NSTTPC
c$$$         WRITE(*,*) 'rtpcout = ',rtpcout
c$$$         WRITE(*,*) 'tpcacri = ',tpcacri
c$$$         WRITE(*,*) 'tpcthbi = ',tpcthbi
c$$$         do i=1,nsttpc
c$$$           ncmat         = ncmat+1
c$$$c(KH)      rcmat(ncmat)  = rtpcinn+(real(i)-0.5)*xstep
c$$$           rcmat(ncmat)  = rtpcinn+tpcthbi+(real(i)-0.5)*xstep
c$$$           zcmin(ncmat)  = -TPCHLFZ
c$$$           zcmax(ncmat)  = TPCHLFZ
c$$$           xrlc(ncmat)   = xstep/xrargon
c$$$           xelosc(ncmat) = xstep*dedxargon
c$$$        enddo
c$$$*
c$$$
c$$$*     outer tube
c$$$         ncmat         = ncmat+1
c$$$         rcmat(ncmat)  = RTPCOUT
c$$$         zcmin(ncmat)  = -TPCHLFZ
c$$$         zcmax(ncmat)  = TPCHLFZ
c$$$         xrlc(ncmat)   = TPCTHBO/xralu
c$$$         xelosc(ncmat) = TPCTHBO*dedxalu
c$$$*     endplate
c$$$         do xs=-1.,1.,2.
c$$$            npmat         = npmat+1
c$$$            zpmat(npmat)  = xs*(TPCHLFZ-TPCTHKE/2.)
c$$$            rpmin(npmat)  = RTPCINN
c$$$            rpmax(npmat)  = RTPCOUT
c$$$            xrlp(npmat)   = 0.35
c$$$            xelosp(npmat) = 0.35*xralu*dedxalu
c$$$         enddo

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

         DO I=1,NEXS
            IF (ITEXTS(I).EQ.0) THEN
               TYPX(I) = 'CYLI'
            ELSE
               TYPX(I) = 'PLAN'
            ENDIF
         ENDDO


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


