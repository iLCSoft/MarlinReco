*  
* $Id: my_link.inc,v 1.2 2006-02-08 15:03:28 owendt Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
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
      INTEGER NMAXLNK
*      PARAMETER(NMAXLNK=25000)
      PARAMETER(NMAXLNK=50000)
      INTEGER NLNK
      REAL RLNK
      COMMON/MY_LINK/NLNK,RLNK(NMAXLNK,16)

