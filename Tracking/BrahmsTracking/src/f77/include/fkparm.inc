*  
* $Id: fkparm.inc,v 1.3 2008-06-26 09:35:56 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2007/11/23 08:51:15  rasp
* New version of parameters
*
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

       integer ntemx,mxts,mxte,mxtk,nexsmx,ncmamx,npmamx,nsrfmx,
     $         nsegmx,npgrmx, nplmamx
c(kh)  PARAMETER (NTEMX=150,MXTS=37,MXTE=50,MXTK=37,NEXSMX=50)
       PARAMETER (NTEMX=LTPDRO,MXTS=37,MXTE=50,MXTK=37,NEXSMX=50,
     & NCMAMX=100,NPMAMX=50,NPLMAMX=200)
       PARAMETER (NSRFMX=NTEMX+NEXSMX+NCMAMX+NPMAMX)
       PARAMETER (NSEGMX=20,NPGRMX=200)

C SJA
	integer p_flag
