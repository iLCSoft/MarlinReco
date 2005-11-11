*  
* $Id: maskdims.inc,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2003/11/28 10:25:09  buesserk
* Added bruser routines which are needed to read in background files.
* Karsten Buesser,
* Nov 28, 2003
*
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
      REAL MKZ1,MKZ2,MKZ2a,MKZ3,MKTHMIN,MKTHMAX,MKZMAX
      REAL MKTHCRT
      REAL MKZI1, MKZI2, MKRI1, MKTHI1, MKDI1, MKTHI2
      REAL MKZG1, MKZG2, MKRG1, MKTHG1
      REAL MKZO1, MKZO2, MKRO1, MKRO2, MKTHO1, MKTHO2
      REAL MKZG3, MKZG4, MKTHG3, MKTHG4,MKIR,MKOR
      COMMON/MASKDIMS/MKZ1,MKZ2,MKZ3,MKTHMIN,MKTHMAX,MKZMAX,
     &  MKTHCRT, MKZI1, MKZI2, MKRI1, MKTHI1, MKDI1,
     &  MKZG1, MKZG2, MKRG1, MKTHG1,MKZ2a,
     &  MKZO1, MKZO2, MKRO1, MKRO2, MKTHO1, MKTHO2,
     &  MKZG3, MKZG4, MKTHG3, MKTHG4, MKTHI2,MKIR,MKOR
      SAVE /MASKDIMS/

