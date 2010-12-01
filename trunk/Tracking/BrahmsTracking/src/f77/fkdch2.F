*  
* $Id: fkdch2.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
***********************************************************************
*                                                                     *
      DOUBLE PRECISION FUNCTION FKDCH2(WG,DP)
*                                                                     *
*   Contribution to chi2 incrementation by a variation DP of param.   *
*                                                                     *
***********************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WG(15),DP(5)
*
      FKDCH2=
     &      WG( 1)*DP(1)**2+WG( 3)*DP(2)**2+
     &      WG( 6)*DP(3)**2+WG(10)*DP(4)**2+WG(15)*DP(5)**2+
     &  2.*(WG( 2)*DP(1)*DP(2)+WG( 4)*DP(1)*DP(3)+
     &      WG( 5)*DP(2)*DP(3)+WG( 7)*DP(1)*DP(4)+
     &      WG( 8)*DP(2)*DP(4)+WG( 9)*DP(3)*DP(4)+
     &      WG(11)*DP(1)*DP(5)+WG(12)*DP(2)*DP(5)+
     &      WG(13)*DP(3)*DP(5)+WG(14)*DP(4)*DP(5))
*
      END
