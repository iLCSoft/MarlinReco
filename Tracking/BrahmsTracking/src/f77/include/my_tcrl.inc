*  
* $Id: my_tcrl.inc,v 1.2 2005-08-09 10:27:50 aplin Exp $
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
      INTEGER NMAXROW
      PARAMETER(NMAXROW=1000)
      INTEGER IHITROW(LTPDRO,NMAXROW),MY_JTCRNC(LTPDRO)
      COMMON/MY_TCRL/IHITROW,MY_JTCRNC

