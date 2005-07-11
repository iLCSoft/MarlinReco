*  
* $Id: delid.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:18  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      INTEGER FUNCTION DELID(BRADET)
C****************************************************************************
C* CONVERT BRAHMS DETECTOR ID TO DELPHI CONVENTION                          *
C*                                                                          *
C* Kristian Harder, October 1999                                            *
C****************************************************************************
      IMPLICIT NONE
#include "include/detnums.inc"
#include "include/fsfparams.inc"

      INTEGER BRADET

      IF (BRADET.EQ.ID_VTX) THEN
        DELID=I_VTX
      ELSE IF (BRADET.EQ.ID_VTX+1) THEN
        DELID=I_VTX1
      ELSE IF (BRADET.EQ.ID_VTX+2) THEN
        DELID=I_VTX2
      ELSE IF (BRADET.EQ.ID_VTX+3) THEN
        DELID=I_VTX3
      ELSE IF (BRADET.EQ.ID_VTX+4) THEN
        DELID=I_VTX4
      ELSE IF (BRADET.EQ.ID_VTX+5) THEN
        DELID=I_VTX5
      ELSE IF (BRADET.EQ.ID_VTXP) THEN
        DELID=I_VTX
      ELSE IF (BRADET.EQ.ID_VTXP1+1) THEN
        DELID=I_VTX1
      ELSE IF (BRADET.EQ.ID_VTXP1+2) THEN
        DELID=I_VTX2
      ELSE IF (BRADET.EQ.ID_VTXP1+3) THEN
        DELID=I_VTX3
      ELSE IF (BRADET.EQ.ID_VTXP2+2) THEN
        DELID=I_VTX4
      ELSE IF (BRADET.EQ.ID_FTD) THEN
        DELID=I_FTD
      ELSE IF (BRADET.EQ.ID_FTD+1) THEN
        DELID=I_FTD1
      ELSE IF (BRADET.EQ.ID_FTD+2) THEN
        DELID=I_FTD2
      ELSE IF (BRADET.EQ.ID_FTD+3) THEN
        DELID=I_FTD3
      ELSE IF (BRADET.EQ.ID_FTD+4) THEN
        DELID=I_FTD4
      ELSE IF (BRADET.EQ.ID_FTD+5) THEN
        DELID=I_FTD5
      ELSE IF (BRADET.EQ.ID_FTD+6) THEN
        DELID=I_FTD6
      ELSE IF (BRADET.EQ.ID_FTD+7) THEN
        DELID=I_FTD7
      ELSE IF (BRADET.EQ.ID_ITC) THEN
        DELID=I_ITC
      ELSE IF (BRADET.EQ.ID_SIT+1) THEN
        DELID=I_SIT1
      ELSE IF (BRADET.EQ.ID_SIT+2) THEN
        DELID=I_SIT2
      ELSE IF (BRADET.EQ.ID_TPC) THEN
        DELID=I_TPC
      ELSE IF (BRADET.EQ.ID_FCH) THEN
        DELID=I_FCH
      ELSE
        DELID=0
      END IF

      RETURN
      END

*
* patch brtrak provides an interface to the DELPHI track fitting routines
*

