*  
* $Id: allpipe.inc,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
      INTEGER NBP
      REAL BPP3(5), BPP4(5)
      REAL BPP5(3), BPP6(4), BPP7(3), BPP8(3), BPP9(5)
      REAL BPP10(3), BPP11(5), BPP12(3)
      REAL BPT0, BPT1, BPT2, BPT3, BPT4, BPT5
      REAL BPR(10,50), BPZ(10,50)
      INTEGER BPM(5,50)
      CHARACTER*4 BPT(2,50)
      COMMON / ALLPIPE / BPP3, BPP4,
     &  BPT0, BPT1, BPT2, BPT3, BPT4, BPT5,
     &  BPP5, BPP6, BPP7, BPP8, BPP9, BPP10, BPP11, BPP12,
     &  BPR, BPZ, NBP, BPM, BPT

