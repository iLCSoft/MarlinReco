*  
* $Id: fxprob.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      FUNCTION FXPROB(CH,N)
      implicit none
      REAL CH,FXPROB,PROB
      INTEGER N

*     Calls prob wihout crashing if CH too big or N too small
c MVL
      IF(CH.GT.0.0.and.CH.lT.300.0.and.N.GT.0)THEN
        FXPROB=max(PROB(CH,N),1.e-25)
      ELSE
        ! print*,'FXPROB: P set to 1e-25.'
        FXPROB=1.e-25
      ENDIF
      END


