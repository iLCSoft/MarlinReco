*  
* $Id: fkmi4d.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE FKMI4D(WGT,COV,IERR)
*                                                                     *
*   Inversion of a (4x4) symmetric positive matrix (Cholesky method)  *
*   Internal computation in double precision                          *
*   Check on positivity (IERR=2 if WGT not positive)                  *
*   COV (output) may overwrite WGT (input)                            *
*                                                                     *
***********************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION WGT(10),COV(10)
*
      IERR=0
*
      IF(WGT(1).LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T11=1./DSQRT(WGT(1))
      S12=WGT( 2)*T11
      S13=WGT( 4)*T11
      S14=WGT( 7)*T11
*
      T22=WGT(3)-S12**2
      IF(T22.LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T22=1./DSQRT(T22)
      S23=(WGT( 5)-S12*S13)*T22
      S24=(WGT( 8)-S12*S14)*T22
*
      T33=WGT(6)-S13**2-S23**2
      IF(T33.LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T33=1./DSQRT(T33)
      S34=(WGT( 9)-S13*S14-S23*S24)*T33
*
      T44=WGT(10)-S14**2-S24**2-S34**2
      IF(T44.LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T44=1./DSQRT(T44)
*
      T34=-T33*(S34*T44)
      T23=-T22*(S23*T33)
      T24=-T22*(S23*T34+S24*T44)
      T12=-T11*(S12*T22)
      T13=-T11*(S12*T23+S13*T33)
      T14=-T11*(S12*T24+S13*T34+S14*T44)
*
      COV( 1)=T11*T11+T12*T12+T13*T13+T14*T14
      COV( 2)=T12*T22+T13*T23+T14*T24
      COV( 3)=T22*T22+T23*T23+T24*T24
      COV( 4)=T13*T33+T14*T34
      COV( 5)=T23*T33+T24*T34
      COV( 6)=T33*T33+T34*T34
      COV( 7)=T14*T44
      COV( 8)=T24*T44
      COV( 9)=T34*T44
      COV(10)=T44*T44
*
      END
