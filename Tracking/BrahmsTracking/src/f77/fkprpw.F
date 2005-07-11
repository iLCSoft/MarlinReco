*  
* $Id: fkprpw.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:18  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
************************************************************************
*                                                                      *
      SUBROUTINE FKPRPW(WGI,DER,WGF)
*                                                                      *
*   propagation of weight matrix from a plane to another one           *
*                                                                      *
*   Input  : WGI  : initial weight matrix                              *
*            DER  : inverse derivatives                                *
*   Output : WGF  : final weight matrix (may overwrite DPARI)          *
*                                                                      *
************************************************************************
*
*   propagation of weight matrix from a plane to a plane
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DER(8),WGI(15),WGF(15)
*
*   computation of   W  =  DER transp * WG
*
      W55    =WGI(11)*DER(3)+WGI(12)*DER(6)+WGI(14)*DER(8)+WGI(15)
      W54    =WGI( 7)*DER(3)+WGI( 8)*DER(6)+WGI(10)*DER(8)+WGI(14)
      W53    =WGI( 4)*DER(3)+WGI( 5)*DER(6)+WGI( 9)*DER(8)+WGI(13)
      W52    =WGI( 2)*DER(3)+WGI( 3)*DER(6)+WGI( 8)*DER(8)+WGI(12)
      W51    =WGI( 1)*DER(3)+WGI( 2)*DER(6)+WGI( 7)*DER(8)+WGI(11)
*
      W44   =WGI(7)*DER(2)+WGI(8)*DER(5)+WGI(10)
      W43   =WGI(4)*DER(2)+WGI(5)*DER(5)+WGI( 9)
      W42   =WGI(2)*DER(2)+WGI(3)*DER(5)+WGI( 8)
      W41   =WGI(1)*DER(2)+WGI(2)*DER(5)+WGI( 7)
*
      W34   =WGI(7)*DER(1)+WGI(8)*DER(4)+WGI(9)+WGI(10)*DER(7)
      W33   =WGI(4)*DER(1)+WGI(5)*DER(4)+WGI(6)+WGI( 9)*DER(7)
      W32   =WGI(2)*DER(1)+WGI(3)*DER(4)+WGI(5)+WGI( 8)*DER(7)
      W31   =WGI(1)*DER(1)+WGI(2)*DER(4)+WGI(4)+WGI( 7)*DER(7)
*
*   computation of   W * DER
*
      WGF(1)=WGI(1)
      WGF(2)=WGI(2)
      WGF(3)=WGI(3)
*
      WGF( 4)=W31
      WGF( 5)=W32
      WGF( 6)=W31*DER(1)+W32*DER(4)+W33+W34*DER(7)
*
      WGF( 7)=W41
      WGF( 8)=W42
      WGF( 9)=W41*DER(1)+W42*DER(4)+W43+W44*DER(7)
      WGF(10)=W41*DER(2)+W42*DER(5)+W44
*
      WGF(11)=W51
      WGF(12)=W52
      WGF(13)=W51*DER(1)+W52*DER(4)+W53+W54*DER(7)
      WGF(14)=W51*DER(2)+W52*DER(5)+W54
      WGF(15)=W51*DER(3)+W52*DER(6)+W54*DER(8)+W55
*
      END
