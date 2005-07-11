*  
* $Id: fsfparams.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:49  hvogt
* Brahms V308
*
*  
#include "fkparm.inc"
#include "bkgmode.inc"
#include "fsfpar1.inc"
      integer maxpd
      parameter (maxpd=8)      ! number of primary detectors (first maxpd ones
                               ! in list below !!!!!
      integer detmax,tssize,maxtev
      parameter (detmax=23)    ! number of detectors we know of
      parameter (tssize=MXTS)  !(max)size of TS stucture.
      parameter (maxtev=detmax+20) !max number of TE we can handle per TS.
c --- detector number inside DELFIT/DELSEARCH !
c     put all primary detectors to the top of the list
      integer i_VTX,i_FTD,i_ITC,i_SIT1,i_SIT2,i_TPC,i_FCH
      integer i_VTX1,i_VTX2,i_VTX3,i_VTX4,i_VTX5,i_SCAN1,i_SCAN2,i_SCAN3
      integer i_SCAN4,i_FTD1,i_FTD2,i_FTD3,i_FTD4,i_FTD5,i_FTD6,i_FTD7
      parameter (i_TPC   =  1 )
      parameter (i_FTD   =  2 )
      parameter (i_ITC   =  3 )
      parameter (i_VTX   =  4 )
      parameter (i_SCAN1 =  5 )
      parameter (i_SCAN2 =  6 )
      parameter (i_SCAN3 =  7 )
      parameter (i_SCAN4 =  8 )
      parameter (i_FCH   =  9 )
      parameter (i_SIT1  = 10 )
      parameter (i_SIT2  = 11 )
      parameter (i_VTX1  = 12 )
      parameter (i_VTX2  = 13 )
      parameter (i_VTX3  = 14 )
      parameter (i_VTX4  = 15 )
      parameter (i_VTX5  = 16 )
      parameter (i_FTD1  = 17 )
      parameter (i_FTD2  = 18 )
      parameter (i_FTD3  = 19 )
      parameter (i_FTD4  = 20 )
      parameter (i_FTD5  = 21 )
      parameter (i_FTD6  = 22 )
      parameter (i_FTD7  = 23 )
c --- The order of the detectors in printout
      integer detorder(detmax)
      data detorder( 1)  /i_SCAN1/
      data detorder( 2)  /i_SCAN2/
      data detorder( 3)  /i_SCAN3/
      data detorder( 4)  /i_SCAN4/
      data detorder( 5)  /i_VTX1/
      data detorder( 6)  /i_VTX2/
      data detorder( 7)  /i_VTX3/
      data detorder( 8)  /i_VTX4/
      data detorder( 9)  /i_VTX5/
      data detorder(10)  /i_VTX/
      data detorder(11)  /i_SIT1/
      data detorder(12)  /i_SIT2/
      data detorder(13)  /i_ITC/
      data detorder(14)  /i_FTD1/
      data detorder(15)  /i_FTD2/
      data detorder(16)  /i_FTD3/
      data detorder(17)  /i_FTD4/
      data detorder(18)  /i_FTD5/
      data detorder(19)  /i_FTD6/
      data detorder(20)  /i_FTD7/
      data detorder(21)  /i_FTD/
      data detorder(22)  /i_TPC/
      data detorder(23)  /i_FCH/
c --- names of detectors
      character *5,detname(detmax)
      data detname(i_TPC)   /'TPC'/
      data detname(i_FTD)   /'FTD'/
      data detname(i_ITC)   /'ITC'/
      data detname(i_VTX)   /'VTX'/
      data detname(i_SCAN1) /'SCAN1'/
      data detname(i_SCAN2) /'SCAN2'/
      data detname(i_SCAN3) /'SCAN3'/
      data detname(i_SCAN4) /'SCAN4'/
      data detname(i_FCH)   /'FCH'/
      data detname(i_SIT1)  /'SIT1'/
      data detname(i_SIT2)  /'SIT2'/
      data detname(i_VTX1)  /'VTX1'/
      data detname(i_VTX2)  /'VTX2'/
      data detname(i_VTX3)  /'VTX3'/
      data detname(i_VTX4)  /'VTX4'/
      data detname(i_VTX5)  /'VTX5'/
      data detname(i_FTD1)  /'FTD1'/
      data detname(i_FTD2)  /'FTD2'/
      data detname(i_FTD3)  /'FTD3'/
      data detname(i_FTD4)  /'FTD4'/
      data detname(i_FTD5)  /'FTD5'/
      data detname(i_FTD6)  /'FTD6'/
      data detname(i_FTD7)  /'FTD7'/
c --- tanagra bits for detectors for DELSEARCH detector index
c --- attention: some bits are reserved (1=primary vertex,
c --- detid(i_TPC) must be equal to detid_shower) - see KEEP,FSFBITS.  -KH-
      integer detid(detmax)
      data detid(i_TPC)  / 5/
      data detid(i_FTD)  / 3/
      data detid(i_ITC)  / 4/
      data detid(i_VTX)  / 6/
      data detid(i_SCAN1)/ 0/
      data detid(i_SCAN2)/ 0/
      data detid(i_SCAN3)/ 0/
      data detid(i_SCAN4)/ 0/
      data detid(i_FCH)  / 7/
      data detid(i_SIT1) / 8/
      data detid(i_SIT2) / 9/
      data detid(i_VTX1) /10/
      data detid(i_VTX2) /11/
      data detid(i_VTX3) /12/
      data detid(i_VTX4) /13/
      data detid(i_VTX5) /14/
      data detid(i_FTD1) /15/
      data detid(i_FTD2) /16/
      data detid(i_FTD3) /17/
      data detid(i_FTD4) /18/
      data detid(i_FTD5) /19/
      data detid(i_FTD6) /20/
      data detid(i_FTD7) /21/
c --- magic number for pseudo-TEs
      integer    maxpseudoTE
      PARAMETER (maxpseudoTE = temax)

