*  
* $Id: vxdgeom.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:50  hvogt
* Brahms V308
*
*  
      SUBROUTINE VXDGEOM
**********************************************************************
*
*     VXDGEOM: routine to define the geometry of the
*              CCD vertex detector
*
*
*
*     author: Ties Behnke
*     created: 4/12/98
*     updated dimensions for CCD vertex detector
*     as communicated by Chris Damerell Aug 9, 1999
*     updated again 12/99 to reflect 14mm beampipe
*
**********************************************************************
      IMPLICIT NONE
c#include "turnon.inc"
#include "include/vxdccd.inc"
#include "include/ccdpipe.inc"

      real beampipe
      external beampipe

c      IF ( .NOT. LVXDC ) GOTO 999

      IF ( BPCINR1 .EQ. 0. ) THEN
        WRITE(*,*) ' Trying to build CCD VTX detector without'
        WRITE(*,*) ' first building the beampipe, abort!'
        GOTO 999
      END IF

C Set up the CCD VXD

C Radiation length of silicon is taken to be 9.36 cm
C Radiation length of beryllium is taken to be 35.28 cm
c Radiation length of steel is taken to be 1.76 cm
c Radiation length of aluminium is taken to be 8.9 cm

c Thickness of active silicon layers (appropriate for CCD option)
      VXDCLT=0.4E-3 * 9.36
C thickness of support structure (Be)
      VXDSUP=0.8E-3 * 35.28
C thickness of electronics structure
      VXDELE=2.1E-3 * 9.36

C thickness of strip lines (use Kapton at the moment)
      VXSTRD = 0.33E-3 * 28.6

C     updated dimensions for CCD vertex detector
C     as communicated by Chris Damerell Dec 1 1998


C End plate + cables dimensions (treat as aluminium for now):
C  radial
      VXDEPR1=2.2
      VXDEPR2=6.5
c  z
      VXDEPZ1=13.5
      VXDEPZ2=13.5
*      VXDEPT=16.7E-3 * 8.9
      VXDEPT = 0.0

c cryostat dimensions:
c  two radial,z, thickness
      VXDCR1=2.32
      VXDCR2=9.0
      VXDCRZ=16.0
*      VXDCRT=4.4E-3 * 8.9
      VXDCRT1=1.0
      VXDCRT2=0.05
*      VXDCRT = 1.

c support structure radial, thickness
      VXDSUPR=6.5
      VXDSUPZ=13.5
      VXDSUPT=1.4E-3 * 35.28

C support structure, endplates
      VXSUPER1 = 2.32
      VXSUPER2 = VXDSUPR+VXDSUPT
      VXSUPEZ1 = 13.5
      VXSUPEZ2 = VXSUPEZ1 + VXDSUPT

C Local electronics:
      VXDELZ1=180.
      VXDELZ2=230.
      VXDELT=1.0
C strip lines from the detector to the local electronics
C     at the moment only implemented to the start of the
C     tungsten part of the mask

C     first check that the beampipe has been properly defined
      IF ( BPCINR3 .NE. 0. ) THEN

        VXDSTZ1(1) = BPCINZ1
        VXDSTZ2(1) = BPCINZ2
        VXDSTR1(1) = BPCINR1+0.1
        VXDSTR2(1) = BPCINR2+0.1

        VXDSTT(1)  = 0.33E-3 * 28.6

        VXDSTZ1(2) = BPCINZ2
        VXDSTZ2(2) = BPCINZ3
        VXDSTR1(2) = BPCINR2+0.1
        VXDSTR2(2) = BPCINR3+0.1

        VXDSTT(2)  = 0.33E-3 * 28.6

        VXDSTZ1(3) = BPCINZ3
        VXDSTZ2(3) = BPCINZ4
        VXDSTR1(3) = BPCINR3+0.1
        VXDSTR2(3) = BPCINR4+0.1

        VXDSTT(2)  = (0.33+0.33+0.50+0.66+0.82)*1E-3 * 28.6

      END IF

c Radii and z-dimensions of successive layers:
c radial

*      VXDRL1=1.2
*      VXDRL2=2.4
*      VXDRL3=3.6
*      VXDRL4=4.8
*      VXDRL5=6.0
c
C    updated for the new 14mm beampipe decided December 99
C    as default version

      VXDRL1=1.55
      VXDRL2=2.7
      VXDRL3=3.8
      VXDRL4=4.9
      VXDRL5=6.0
c
C     half length of detector in z
C
CCC      VXDZL1=2.5
      VXDZL1=5.0
      VXDZL2=12.5
      VXDZL3=12.5
      VXDZL4=12.5
      VXDZL5=12.5

C     gap between ladders in the middle of the detector
      VXDGA1=0.00
      VXDGA2=0.004
      VXDGA3=0.004
      VXDGA4=0.004
      VXDGA5=0.004

C     support shells: thickness in cm Be
      VXDSUD1 = 272E-4
      VXDSUD2 = 272E-4
      VXDSUD3 = 272E-4
      VXDSUD4 = 272E-4
      VXDSUD5 = 272E-4

C     z positions of electronics at the end of the detector
C     radii (inner, outer)
C     thickness

CCC      VXDELZ11 = 2.5
CCC      VXDELZ12 = 3.5
      VXDELZ11 = 5.0
      VXDELZ12 = 6.0
*      VXDELR11 = 1.2
*      VXDELR12 = 1.2
      VXDELR11 = 1.55
      VXDELR12 = 1.55

      VXDELDZ1 = VXDELE

C     strip lines

C     *** layer 1 ***

      VXDSTZ11 = VXDELZ12+VXDEPT
      VXDSTZ12 = VXDSTZ11+1.1+VXDEPT
      VXDSTR11 = 1.55
*      VXDSTR12 = 2.4

      vxdstr12 = 1.7
c      vxdstr12 = beampipe(vxdstz12)+0.05
csja FIXME: taken straight out of BRAHMS for CCD (turn on vtx=12)
      vxdstr12 = 1.54999995
      VXDSTD1  = VXSTRD

C     electronics

      VXDELZ21 = 12.5
      VXDELZ22 = 13.5
      VXDELR21 = 2.7
      VXDELR22 = 2.7
      VXDELDZ2 = VXDELE

C     strip lines

      VXDSTZ21 = 13.5+VXDSUPT+VXDEPT
      VXDSTZ22 = 13.6+VXDEPT
      VXDSTR21 = 2.7
      VXDSTR22 = 2.3
      VXDSTD2  = VXSTRD

C     electronics

      VXDELZ31 = 12.5
      VXDELZ32 = 13.5
      VXDELR31 = 3.8
      VXDELR32 = 3.8
      VXDELDZ3 = VXDELE

C     strip lines

      VXDSTZ31 = 13.5+VXDSUPT+VXDEPT
      VXDSTZ32 = 14.0+VXDEPT
      VXDSTR31 = 3.8
      VXDSTR32 = 2.3
      VXDSTD3  = VXSTRD

C     electronics

      VXDELZ41 = 12.5
      VXDELZ42 = 13.5
      VXDELR41 = 4.9
      VXDELR42 = 4.9
      VXDELDZ4 = VXDELE

C     strip lines

      VXDSTZ41 = 13.5+VXDSUPT+VXDEPT
      VXDSTZ42 = 14.5+VXDEPT
      VXDSTR41 = 4.9
      VXDSTR42 = 2.3
      VXDSTD4  = VXSTRD

C     electronics

      VXDELZ51 = 12.5
      VXDELZ52 = 13.5
      VXDELR51 = 6.0
      VXDELR52 = 6.0
      VXDELDZ5 = VXDELE

C     strip lines

      VXDSTZ51 = 13.5+VXDSUPT+VXDEPT
      VXDSTZ52 = 15.0+VXDEPT
      VXDSTR51 = 6.0
      VXDSTR52 = 2.30
      VXDSTD5  = VXSTRD


C Point resolution(cm):
      VXDCPNT=3.5E-4

  999 return
      end


