*  
* $Id: sitgeom.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:51  hvogt
* Brahms V308
*
*  
      SUBROUTINE SITGEOM
**********************************************************************
*
*     SUBROUTINE SITGEOM Routine to define the
*                geometry of the sit
*
*
*     author Klaus Moenig
*     created: 15/1/99
*     last modified: 12/10/99 K. Moenig
*
**********************************************************************
      IMPLICIT NONE

c#include "turnon.inc"
#include "include/maskdims.inc"
#include "include/sitdims.inc"

C----------------------------------------------------------------------

c SIT Dimensions:
c THickness of material
      SITTHK = 3.0E-2
C Inner Radius
      SITRAD(1) = 16.
      SITRAD(2) = 30.
c Z-length of ITC:
*      SITLEN(1)= 36.
      SITLEN(1)= 38.
      SITLEN(2)= 66.
*          number of measurements
C For simple studies, use a general r-phi resolution and z-resolution (cm^2):
      SITRPRES=0.001
      SITZRES =0.005

  999 RETURN
      END


