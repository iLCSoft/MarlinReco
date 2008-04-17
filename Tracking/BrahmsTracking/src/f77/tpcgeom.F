
*  
* $Id: tpcgeom.F,v 1.5 2008-04-17 14:36:29 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.4  2006/11/02 11:49:01  rasp
*
* Modified files adopted to FullLDCTracking processor
*
* Revision 1.3  2006/06/28 15:29:04  aplin
* The B-Field is now variable for LEPTracking via the gear xml file. The B-Field is specified in the TPCParameters as follows: <parameter name="tpcBField" type="double"> 4.0  </parameter>
*
* The value is passed internaly to the F77 code via the same function which passes the TPC geometry i.e. gettpcgeom(float* innerrad, float* outerrad, int* npadrows, float* maxdrift, float* tpcpixz, float* ionpoten, float* tpcrpres, float* tpczres, float* tpcbfield). It is set in setmat.F. tpcgeom.F had to be modified as it also uses gettpcgeom, although it does not make use of the B-Field.
*
* Revision 1.2  2005/08/08 07:09:13  aplin
* Made f77 tracking code use GEAR to define the geomtery of the TPC. LTPDRO now defines the maximum number of rows is used to define the size of arrays, this is limited to 224 due the use of 7 '32 bit' bit registers in trkfnd.F increased, though at present it is not likely that anybody would want more. The number of TPC padrows is defined at run time by NRTPC which should of course not exceed LTPDRO, although this is checked and the programe exits with a verbose error message. A wrapper function gettpcgeom is used to pass the GEAR TPC parameters from C++ to f77. MarlinUtil/include/marlin_tpcgeom.h have MarlinUtil/src/marlin_tpcgeom.cc consequently been removed as they are no longer needed.
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
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
c#include "marlin_tpcgeom.h"

      REAL INNERRAD, OUTERRAD, MAXDRIFT, TPCPIXZ, IONPOTEN 
      REAL TPCRPRES, TPCZRES, TPCBFIELD
      
      INTEGER NPADROWS

      INTEGER TRY
      INTEGER gettpcgeom
      external gettpcgeom



C CALL wrapper function to get tpc parameters from GEAR
      TRY = gettpcgeom(INNERRAD, OUTERRAD, NPADROWS, MAXDRIFT, TPCPIXZ, 
     +     IONPOTEN, TPCRPRES, TPCZRES, TPCBFIELD)


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
      TPCHLFZ = MAXDRIFT + 23.

C Maximum drift length:
      ZDRIFT=MAXDRIFT
c Radial pad size
      TPCPADR=(TPCACRO-TPCACRI)/REAL(NRTPC)
      BFIELD =  TPCBFIELD*10.0

C      write(*,*) 'NRTPC ', NRTPC
C      write(*,*) 'TPCACRO ', TPCACRO
C      write(*,*) 'TPCPADR ', TPCPADR 
C      write(*,*) 'RTPCOUT ', RTPCOUT

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
