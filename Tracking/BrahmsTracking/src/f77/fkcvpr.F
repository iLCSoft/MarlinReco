*  
* $Id: fkcvpr.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
************************************************************************
*                                                                      *
      SUBROUTINE FKCVPR(LDEB,COV)
*                                                                      *
*     Level 2 debug output : compute and write on LDEB the errors and  *
*     non-diagonal elements of the correlation matrix, from the        *
*     covariance matrix COV (their abs. value should be less than 1.)  *
*     Computes also for each triplet i,j,k the quantity :              *
*        CORij**2 + CORjk**2 + CORki**2 - 2 * CORij * CORjk * CORki    *
*     (should be less than 1.), and write this quantity                *
*                                                                      *
************************************************************************
*
      DIMENSION COV(15)
      CHARACTER*5 TXT(5)
*
      DATA TXT /'r*Phi','  z  ','theta',' phi ',' 1/p '/
*
      ERR1 = SQRT(COV( 1))
      ERR2 = SQRT(COV( 3))
      ERR3 = SQRT(COV( 6))
      ERR4 = SQRT(COV(10))
      ERR5 = SQRT(COV(15))
      COR12 = COV( 2)/(ERR1*ERR2)
      COR13 = COV( 4)/(ERR1*ERR3)
      COR23 = COV( 5)/(ERR2*ERR3)
      COR14 = COV( 7)/(ERR1*ERR4)
      COR24 = COV( 8)/(ERR2*ERR4)
      COR34 = COV( 9)/(ERR3*ERR4)
      COR15 = COV(11)/(ERR1*ERR5)
      COR25 = COV(12)/(ERR2*ERR5)
      COR35 = COV(13)/(ERR3*ERR5)
      COR45 = COV(14)/(ERR4*ERR5)
*
      WRITE(LDEB,'(24X,5(A5,9X),/,'' errors'',9X,5F14.6,/,'' correl'',
     & 4X,A5,F14.4,/,11X,A5,2F14.4,/,11X,A5,3F14.4,/,11X,A5,4F14.4)')
     & (TXT(K),K=1,5),ERR1,ERR2,ERR3,ERR4,ERR5,
     & TXT(2),COR12,TXT(3),COR13,COR23,TXT(4),COR14,COR24,COR34,
     & TXT(5),COR15,COR25,COR35,COR45
*
      WRITE(LDEB,'('' check    123       124       125       134       1
     &35       145       234       235       245       345'',/,
     & 5X,10F10.6)')
     & COR12**2+COR13**2+COR23**2-2.*COR12*COR13*COR23,
     & COR12**2+COR14**2+COR24**2-2.*COR12*COR14*COR24,
     & COR12**2+COR15**2+COR25**2-2.*COR12*COR15*COR25,
     & COR13**2+COR14**2+COR34**2-2.*COR13*COR14*COR34,
     & COR13**2+COR15**2+COR35**2-2.*COR13*COR15*COR35,
     & COR14**2+COR15**2+COR45**2-2.*COR14*COR15*COR45,
     & COR23**2+COR24**2+COR34**2-2.*COR23*COR24*COR34,
     & COR23**2+COR25**2+COR35**2-2.*COR23*COR25*COR35,
     & COR24**2+COR25**2+COR45**2-2.*COR24*COR25*COR45,
     & COR34**2+COR35**2+COR45**2-2.*COR34*COR35*COR45
*
      END
