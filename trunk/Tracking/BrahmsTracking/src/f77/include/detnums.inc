*  
* $Id: detnums.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:47  hvogt
* Brahms V308
*
*  
****************************************************************
*
*     numbering scheme for subdetectors:
*     n*1000 + m*100 * i*10 + j
*        where:   n  n=0  tracking detector
*                    n=1  calorimeter detectors
*                    n=2  muon chambers
*                 m       number of major detector system:
*                         VTX, ITC, TPC, etc
*                 i       division within one detector:
*                    i=1  barrel
*                    i=2  endcap
*                 j       subdetector components, e.g. layer
*                         number in VTX etc.
*
***************************************************************
*
* tracking detectors:

      INTEGER ID_VTX, ID_ITC, ID_SIT, ID_TPC
      PARAMETER ( ID_VTX=100, ID_ITC=300, ID_SIT=400, ID_TPC=500)

      INTEGER ID_VTXP, ID_VTXP1, ID_VTXP2, ID_VTXP3
      PARAMETER ( ID_VTXP=100)
      PARAMETER ( ID_VTXP1=ID_VTXP+10 , ID_VTXP2=ID_VTXP+20,
     $            ID_VTXP3=ID_VTXP+30)

* forward tracking detectors

      INTEGER ID_FTD, ID_FTD1, ID_FTD2, ID_FTD3, ID_FTD4, ID_FTD5
      INTEGER ID_FTD6, ID_FTD7
C(KH) PARAMETER ( ID_FTD=600)
      PARAMETER ( ID_FTD=200)
      PARAMETER ( ID_FTD1=ID_FTD+1, ID_FTD2=ID_FTD+2,
     &            ID_FTD3=ID_FTD+3, ID_FTD4=ID_FTD+4, ID_FTD5=ID_FTD+5,
     &            ID_FTD6=ID_FTD+6, ID_FTD7=ID_FTD+7)

* calorimeter detectors
      INTEGER ID_ECAL, ID_EBCAL, ID_EECAL
      INTEGER ID_HCAL, ID_HBCAL, ID_HECAL
      integer ID_YOK, ID_BYOK, ID_EYOK
      INTEGER ID_COIL
      PARAMETER ( ID_ECAL=1100)
      PARAMETER ( ID_EBCAL=ID_ECAL+10, ID_EECAL=ID_ECAL+20)
      PARAMETER ( ID_HCAL=1200)
      PARAMETER ( ID_HBCAL=ID_HCAL, ID_HECAL=ID_HCAL+100)
      PARAMETER ( ID_YOK=1400)
      PARAMETER ( ID_BYOK=ID_YOK+10, ID_EYOK=ID_YOK+20)
      PARAMETER ( ID_COIL=1900)

* forward calo type detectors
      INTEGER ID_LAT, ID_LCAL
      PARAMETER ( ID_LAT=1500)
      PARAMETER ( ID_LCAL=ID_LAT+100)

* forward chambers
      INTEGER ID_FCH
C(KH) PARAMETER ( ID_FCH=1500)
      PARAMETER ( ID_FCH=600)

* muon detectors
      INTEGER ID_BMU, ID_EMU
      PARAMETER ( ID_BMU=2110, ID_EMU=2130)
