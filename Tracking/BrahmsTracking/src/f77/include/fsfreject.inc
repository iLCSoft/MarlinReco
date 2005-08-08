*  
* $Id: fsfreject.inc,v 1.1 2005-07-21 16:37:00 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:49  hvogt
* Brahms V308
*
*  
c --- The following strings are to be killed, because theta and z are very bad
c     Use bit combinations of detector ids here (array detid() ).
      INTEGER badcomb(30)
      DATA badcomb(1)  /  0  /            ! The Terminator :-)
c --- the following strings need at least 2 (hit?) TEs in VTX/FTD
c     (number of TEs greater than number of detectors, multi-hit VTX/FTD TEs
c      corrected)
      INTEGER badsingle(30)
      DATA badsingle(1)  / 0 /      ! The Terminator :-)

