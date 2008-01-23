*  
* $Id: my_tarc.inc,v 1.2 2008-01-23 17:22:22 aplin Exp $
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

#include "padrow.inc"

C      INTEGER N32BITREG2
C      PARAMETER(N32BITREG2=8)

      REAL MY_JTARIR,MY_JTARTL,MY_JTARP0,MY_JTARD0,MY_JTARZ0,MY_JTARC1
      REAL MY_JTARC2,MY_JTARER
      INTEGER MY_JTARPI,MY_JTAROF,MY_JTARNC,MY_JTARAN,MY_JTARTN
      COMMON/MY_TARC/
     & MY_JTARIR(MAXCND),MY_JTARTL(MAXCND),MY_JTARP0(MAXCND),
     & MY_JTARD0(MAXCND),MY_JTARZ0(MAXCND),MY_JTARC1(MAXCND),
     & MY_JTARC2(MAXCND),MY_JTARER(MAXCND),MY_JTARPI(MAXCND,N32BITREG),
     & MY_JTAROF(MAXCND),MY_JTARNC(MAXCND),MY_JTARAN(MAXCND),
     & MY_JTARTN(MAXCND)

