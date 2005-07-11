*  
* $Id: fkmsct.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE FKMSCT (SINTH,COSTH,PINV,XRL,WG)
*                                                                     *
*   Transformation of weight matrix to account for multiple scatt.    *
*   at a given point                                                  *
*                                                                     *
*   Input  :  SINTH  :  sin(theta)                                    *
*             COSTH  :  cos(theta)                                    *
*             PINV   :  1/p                                           *
*             XRL    :  number of rad. lengths crossed at this point  *
*             WG     :  initial weight matrix                         *
*                                                                     *
*   Output :  WG     :  final weight matrix                           *
*                                                                     *
***********************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "include/fksave.inc"
      REAL SINTH,COSTH,PINV,XRL
      DIMENSION WG(15)
#include "include/fkfild.inc"
*
*      mean squared scattering angle in a direction perpendicular to
*      the trajectory (assuming large ratio momentum/mass)
*      Correct multiple scattering J.Wickens 12-05-92
       FACTOR=1.
       IF(XRL.GT.0.00001)FACTOR=1.0+0.038*ALOG(XRL)
       TH2=(0.0136*PINV*FACTOR)**2*XRL
*
*   matrix (2x2)  X = (cov.m.s.(theta,beta))**(-1) + TH2 * WG(theta,beta)
*   A = cov(theta,1/R)
      A=-CONSB*PINV*COSTH/SINTH**2
*
      WG4P=WG(4)+A*WG(11)
      WG5P=WG(5)+A*WG(12)
      WG6P=WG(6)+A*WG(13)
      WG9P=WG(9)+A*WG(14)
      WG13P=WG(13)+A*WG(15)
*
      XTT=1.+TH2*(WG6P+A*WG13P)
      XTB=TH2*WG9P
      XBB=SINTH**2+TH2*WG(10)
*
*   matrix  U = TH2 * X**(-1)
      Q=TH2/(XBB*XTT-XTB**2)
      UBB=Q*XTT
      UTB=-Q*XTB
      UTT=Q*XBB
*
*   computation of  WG - WG * (U ext. to 5x5) * WG
*
      S13=WG4P*UTT + WG(7)*UTB
      S14=WG4P*UTB + WG(7)*UBB
      S23=WG5P*UTT + WG(8)*UTB
      S24=WG5P*UTB + WG(8)*UBB
      S33=WG6P*UTT + WG(9)*UTB
      S34=WG6P*UTB + WG(9)*UBB
      S43=WG9P*UTT + WG(10)*UTB
      S44=WG9P*UTB + WG(10)*UBB
      S53=WG13P*UTT+ WG(14)*UTB
      S54=WG13P*UTB+ WG(14)*UBB
*
      WG(1)=WG(1)-S13*WG4P-S14*WG(7)
      WG(2)=WG(2)-S23*WG4P-S24*WG(7)
      WG(4)=WG(4)-S33*WG4P-S34*WG(7)
      WG(11)=WG(11)-S53*WG4P-S54*WG(7)
      WG(7)=WG(7)-S43*WG4P-S44*WG(7)
*
      WG(3)=WG(3)-S23*WG5P-S24*WG(8)
      WG(5)=WG(5)-S33*WG5P-S34*WG(8)
      WG(12)=WG(12)-S53*WG5P-S54*WG(8)
      WG(8)=WG(8)-S43*WG5P-S44*WG(8)
*
      WG(6)=WG(6)-S33*WG6P-S34*WG(9)
      WG(13)=WG(13)-S53*WG6P-S54*WG(9)
      WG(9)=WG(9)-S43*WG6P-S44*WG(9)
*
      WG(15)=WG(15)-S53*WG13P-S54*WG(14)
      WG(14)=WG(14)-S53*WG9P-S54*WG(10)
      WG(10)=WG(10)-S43*WG9P-S44*WG(10)
*
      END
