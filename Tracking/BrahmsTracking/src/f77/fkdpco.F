*  
* $Id: fkdpco.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE FKDPCO(A,B,N)
      DOUBLE PRECISION A(*),B(*)
      DO 1 I=1,N
    1 B(I)=A(I)
      RETURN
      END
