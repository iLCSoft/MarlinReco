*  
* $Id: rflags.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      INTEGER JDBDRF,JHISRF,JPRSRF,JPRERF,MDET, IDEBBC, IDEBTC
      LOGICAL FDETRF,FREPRO,FPASS0
      INTEGER NCHDET
      PARAMETER (MDET=17)
      PARAMETER (NCHDET=4)
      CHARACTER*(NCHDET) TNAMRD
      INTEGER JCMORF,JEBIRF,JFBIRF
      COMMON /RFLAGS/FDETRF(MDET),FREPRO,FPASS0,JDBDRF(MDET),
     &        JHISRF(MDET),JPRSRF(MDET),JPRERF(MDET),JCMORF,
     &        JEBIRF(2),JFBIRF(2),TNAMRD(MDET), IDEBBC, IDEBTC

