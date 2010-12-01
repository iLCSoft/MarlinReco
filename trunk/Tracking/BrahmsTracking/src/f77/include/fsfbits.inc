*  
* $Id: fsfbits.inc,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:49  hvogt
* Brahms V308
*
*  
      INTEGER detid_pvtx,detid_shower,detid_zfudge,
     1        subid_vertex,subid_track
      PARAMETER (detid_pvtx   = 32)
      PARAMETER (detid_zfudge = 31)
      PARAMETER (detid_shower = 30)

c
c use Monte Carlo information for track merging?
c
