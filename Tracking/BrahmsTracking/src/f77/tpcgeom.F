
*  
* $Id: tpcgeom.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2005/06/20 18:49:46  aplin
* change tpc geometery so that it is read from  same header file as cpp i.e. marlin_tpcgeom.h
*
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:50  hvogt
* Brahms V308
*
*  
      SUBROUTINE TPCGEOM
**********************************************************************
*
*     SUBROUTINE TPCGEOM Routine to define the
*                geometry of the time projection chamber
*
*
*     author Ties Behnke
*     created: 08/12/98
*     updated 15/1/99     K.M after consultation of Ron Settles
*
**********************************************************************
      IMPLICIT NONE

C#include "turnon.inc"
C#include "maskdims.inc"
#include "include/tpcdims.inc"
#include "include/tpcres.inc"
#include "include/coildims.inc"
#include "marlin_tpcgeom.h"

C----------------------------------------------------------------------
c TPC Dimensions:
C Total inner thickness of TPC material is 3%X_0 and 1X_0 of Al is 8.9 cm
c THickness of material in the Inner wall (Barrel) is 1%X_0
      TPCTHBI=0.01*8.9
c THickness of material in the Outer wall (Barrel) is 3%X_0
      TPCTHBO=0.03*8.9
c Endplate thickness - 30cm.
      TPCTHKE=30.
*     thinner endplates to accomodate the FCH and to not
*     interfere with the ECAL
      TPCTHKE=23.


C Inner Active radius
      TPCACRI=INNERRAD
C Outer Active radius
      TPCACRO=OUTERRAD
C Number of radial readout rings:      
      NRTPC = NPADROWS

C Inner Radius
      RTPCINN=32.
C Outer Radius      
      RTPCOUT=TPCACRO + 6.4
      
c     z-half-length of TPC:
      TPCHLFZ=273.

C Maximum drift length:
      ZDRIFT=MAXDRIFT
c Radial pad size
      TPCPADR=(TPCACRO-TPCACRI)/REAL(NRTPC)

      write(*,*) 'NRTPC ', NRTPC
      write(*,*) 'TPCACRO ', TPCACRO
      write(*,*) 'TPCPADR ', TPCPADR 
      write(*,*) 'RTPCOUT ', RTPCOUT

C----------------------------------------------------------------------
C TPC Response:
c Pad Resolution function (cm)
      TPCPRF=TPCPIXZ
c Radial Resolution (Pad radial size/sqrt(12))(cm)
      TPCRRS=TPCPADR/SQRT(12.)
c Transverse drift parameter (cm^(1/2))
      TPCDT0=600.E-4
c Longitudinal drift parameter (cm^(1/2)
      TPCDL0=380.E-4
c Omega tau/B  (Cyclotron frequency parameter) (T^(-1))
      TPCWTB=5.9
c Excitation energy for a single electron (GeV) (10*Z eV; Z=18 for Argon)
      TPCEION=IONPOTEN
C  Set up the transverse diffusion coefficient (Factor of 10 to convert
c  B-field  from kG to T)
      TPCDT=TPCDT0/SQRT(1.+(TPCWTB*BFIELD/10.)**2)
C For simple studies, use a general r-phi resolution and z-resolution (cm^2):
      TPRP=TPCRPRES
      TPZ=TPCZRES
      TPRPSQ=TPRP**2
      TPZSQ=TPZ**2
      TPRRSQ=TPCRRS**2

  999 RETURN
      END
