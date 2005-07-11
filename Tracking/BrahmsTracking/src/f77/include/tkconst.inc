*  
* $Id: tkconst.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
#include "tkcnst1.inc"
**************************************************************************
*
* array boundary definitions
*
**************************************************************************

* --- maximum number of track candidates per event
      INTEGER NTSMX
      PARAMETER (NTSMX=2000)

* --- maximum number of tracks per event
      INTEGER NTKMX
      PARAMETER (NTKMX=1000)

* --- maximum size of exclusion list
      INTEGER MXEXCL
      PARAMETER (MXEXCL=22000)

* --- number of words describing a hit
      INTEGER MXHT
      PARAMETER (MXHT=11)


