*  
* $Id: fkdpv0.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:18  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:56  hvogt
* Brahms V308
*
*  
      SUBROUTINE FKDPV0(A,N)
      DOUBLE PRECISION A(*)
C      DO 1 I=1,N
C    1 A(I)=0.
      CALL VZERO(A,N*2)
      RETURN
      END
