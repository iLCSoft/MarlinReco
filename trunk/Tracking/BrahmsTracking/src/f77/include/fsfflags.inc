*  
* $Id: fsfflags.inc,v 1.1 2005-07-21 16:37:00 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:49  hvogt
* Brahms V308
*
*  
      integer FLAG_DATA,FLAG_LABL,FBIT_LABL
      integer FLAG_BEAMSPOT,FBIT_BEAMSPOT
      integer FLAG_USEUSED,FBIT_USEUSED
      integer FLAG_GIVESINGLTE,FBIT_GIVESINGLTE
      integer FLAG_OUTLAYER,FBIT_OUTLAYER
      integer FLAG_AVERAGETHETA,FBIT_AVERAGETHETA
      integer FLAG_PCONST,FBIT_PCONST
      integer FLAG_NOPLINV,FBIT_NOPLINV
      integer FLAG_SCANDET,FBIT_SCANDET
      integer FLAG_NOSCALING,FBIT_NOSCALING
      integer FLAG_LIFTTK,FBIT_LIFTTK
      integer FLAG_CUT100,FBIT_CUT100
      integer FLAG_EXTRAPOL,FBIT_EXTRAPOL
      parameter (FLAG_DATA        =0)
      PARAMETER (FLAG_LABL        =1   ,FBIT_LABL        =1)
      parameter (FLAG_BEAMSPOT    =4   ,FBIT_BEAMSPOT    =3)
      parameter (FLAG_USEUSED     =8   ,FBIT_USEUSED     =4)
      parameter (FLAG_GIVESINGLTE =16  ,FBIT_GIVESINGLTE =5)
      parameter (FLAG_OUTLAYER    =32  ,FBIT_OUTLAYER    =6)
      parameter (FLAG_AVERAGETHETA=64  ,FBIT_AVERAGETHETA=7)
      parameter (FLAG_PCONST      =128 ,FBIT_PCONST      =8)
      parameter (FLAG_NOPLINV     =256 ,FBIT_NOPLINV     =9)
      parameter (FLAG_SCANDET     =512 ,FBIT_SCANDET     =10)
      parameter (FLAG_NOSCALING   =1024,FBIT_NOSCALING   =11)
      parameter (FLAG_LIFTTK      =2048,FBIT_LIFTTK      =12)
      parameter (FLAG_CUT100      =4096,FBIT_CUT100      =13)
      parameter (FLAG_EXTRAPOL    =8192,FBIT_EXTRAPOL    =14)
*
      integer FLAG_DBL,FLAG_SUB
      parameter (FLAG_DBL=1,FLAG_SUB=2)

