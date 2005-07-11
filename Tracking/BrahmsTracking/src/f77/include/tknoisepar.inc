*  
* $Id: tknoisepar.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
***********************************************************
*
* background hits per bunch crossing for tracking detectors
*
***********************************************************

      INTEGER TNVXL1,TNVXL2,TNVXL3,TNVXL4,TNVXL5
      INTEGER TNVTX1,TNVTX2,TNVTX3,TNVF2A
      INTEGER TNFTD1,TNFTD2,TNFTD3,TNFTD4,TNFTD5,TNFTD6,TNFTD7
      INTEGER TNSIT1,TNSIT2
      INTEGER TNFCH
      INTEGER TNNBX,TNENERGY

      COMMON /TKNOISEPAR/
     >        TNVXL1,TNVXL2,TNVXL3,TNVXL4,TNVXL5,
     >        TNVTX1,TNVTX2,TNVTX3,TNVF2A,
     >        TNFTD1,TNFTD2,TNFTD3,TNFTD4,TNFTD5,TNFTD6,TNFTD7,
     >        TNSIT1,TNSIT2,
     >        TNFCH,
     >        TNNBX,TNENERGY
      SAVE /TKNOISEPAR/



