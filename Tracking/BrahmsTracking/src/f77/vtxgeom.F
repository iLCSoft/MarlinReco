
      subroutine vtxgeom
*******************************************************************************
*
*      This routine will be responsibe for the geometric information needed
*      for the fitting of track with vtx hits. The information will be taken from
*      Brahms v3-1-2.
*
*      WRITTEN BY                                  Steve Aplin
*******************************************************************************
      IMPLICIT NONE

#include "include/ipsmear.inc"

c     for 500 GeV machine:
C     Interaction point smearing in x,y,z (cm):
C     as of TESLA machine specifications from Fall 98
      SMRIPX=553.0 E-7
      SMRIPY=5.0 E-7
      SMRIPZ=0.04
    
      END
