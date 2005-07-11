*  
* $Id: trfdct.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:36  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
      REAL TDZTOL,TZ0TOL,TD0TOL,RHOMIN
      REAL TDWMAX,TDZMAX,TDWSTR,TDZSTR
      INTEGER MAXGAP
      REAL DZNEXT,DPNEXT,CURTOL,DIPTOL
      INTEGER MAXTRY
      REAL SIGMAW,SFACTR,SIGMAZ,CHLIM1,CHLIM2
      REAL TDCHM1,TDCHM2
      COMMON /TRFDCT/ TDZTOL,TZ0TOL,TD0TOL,RHOMIN
     1               ,TDWMAX,TDZMAX,TDWSTR,TDZSTR
     2               ,MAXGAP,DZNEXT,DPNEXT,CURTOL,DIPTOL
     3               ,MAXTRY,SIGMAW,SFACTR,SIGMAZ,CHLIM1,CHLIM2
     4               ,TDCHM1,TDCHM2

