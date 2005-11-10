*  
* $Id: tkperfc.inc,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
**************************************************************************
*
* performance measurement variables
*
**************************************************************************

      REAL    PFPMC(NTKMX)      ! momentum of all tree tracks
      REAL    PFCTMC(NTKMX)     ! |cos(theta)| of all tree tracks

      INTEGER PFLOPA(9,DETMAX)  ! local pattern recognition performance
      INTEGER PFDELS(9)         ! track search (DELSEARCH) performance
      INTEGER PFGLOB(9)         ! global (DELSEARCH+DELAMBI) performance
      INTEGER PFTETS(5,DETMAX)  ! TE contribution to TS
      INTEGER PFTETK(5,DETMAX)  ! TE contribution to TK

      REAL PFFAIL(5)         ! classification of inefficiencies
      REAL PFFAKE(5)         ! classification of fake tracks
      REAL PFSPLIT(13)       ! classification of track splittings

      COMMON /TKPF/ PFPMC,PFCTMC,PFLOPA,PFDELS,PFGLOB,PFTETS,PFTETK,
     >              PFFAIL,PFFAKE,PFSPLIT
      SAVE   /TKPF/


**************************************************************************
*
* TKSTEER patch
*
* steering of pattern recognition and track finding in BRAHMS
*
* Kristian Harder, September 1999
*
**************************************************************************
*
* routines:
* ---------
*
* GUTREV      -  GEANT user routine (replaces BRAHMS original version)
* TKINIT      -  initialisation routine (definition of FFREAD cards etc.)
* TKSETR      -  second initialisation routine to be called after FFGO
* TKTREV      -  main steering routine
* TKFIN       -  final routine (histogram writing etc.)
*
* TKNUMB      -  get number of hit/TE/TS/TK entries in TKBANK
* TKREAD      -  read hit/TE/TS/TK entry from TKBANK
* TKMKTE      -  create TE from hit list only (calls TKSTOR and TKEXTR)
* TKSTOR      -  store TE/TS/TK in TKBANK
* TKUPDA      -  replace hit in TKBANK by new data (dirty, do not use it!)
* TKEXTR      -  add extrapolation data to existing TE
* TKADDH      -  add left-over hits from a certain subdetector as TEs
* TKGLAB      -  get MC label for TE/TS/TK entry in TKBANK
* TKLIST      -  get component (hit,TE) list for TKBANK entry (TE,TS/TK)
* TKMKEX      -  create list of mutual exclusions for TE/TS list in TKBANK
* TKEXCL      -  check if two TE/TS are mutually exclusive
* TKLSEX      -  get list of all TE/TS excluded against a given TE/TS
*
* TKGHIT      -  general routine to create hit and GEANT track bank
* TKHVTX      -  TKGHIT calling routine for VTX hits
* TKHFTD      -  TKGHIT calling routine for FTD hits
* TKHITC      -  TKGHIT calling routine for ITC hits
* TKHSIT      -  TKGHIT calling routine for SIT hits
* TKHTPC      -  TKGHIT calling routine for TPC hits
* TKHPRE      -  TKGHIT calling routine for PRE hits
*
***************************************************************************


