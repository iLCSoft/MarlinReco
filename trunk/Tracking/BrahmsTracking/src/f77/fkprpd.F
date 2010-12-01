*  
* $Id: fkprpd.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE FKPRPD(DPARI,DER,DPARF)
*                                                                      *
*   propagation of deviations w.r.t. a reference trajectory,           *
*   from a plane to another one (linear approximation)                 *
*                                                                      *
*   Input  : DPARI  : initial deviations                               *
*            DER    : derivatives (see definition in FKXDCY)           *
*   Output : DPARF  : final deviations (may overwrite DPARI)           *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DPARI(5),DPARF(5),DER(8)
*
      DPARF(1)=DPARI(1)+DER(1)*DPARI(3)+DER(2)*DPARI(4)+DER(3)*DPARI(5)
      DPARF(2)=DPARI(2)+DER(4)*DPARI(3)+DER(5)*DPARI(4)+DER(6)*DPARI(5)
      DPARF(3)=DPARI(3)
      DPARF(4)=DER(7)*DPARI(3)+DPARI(4)+DER(8)*DPARI(5)
      DPARF(5)=DPARI(5)
*
      END
