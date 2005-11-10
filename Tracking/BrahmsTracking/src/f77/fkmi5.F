*  
* $Id: fkmi5.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
***********************************************************************
*                                                                     *
      SUBROUTINE FKMI5(WGT,COV,IERR)
*                                                                     *
*   Inversion of a (5x5) symmetric positive matrix (Cholesky method)  *
*   Internal computation in double precision                          *
*   Check on positivity (IERR=2 if WGT not positive)                  *
*   COV (output) may overwrite WGT (input)                            *
*                                                                     *
***********************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION WGT(15),COV(15)
      REAL WGT,COV
*
      IERR=0
*
      IF(WGT(1).LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T11=1./SQRT(WGT(1))
      S12=WGT( 2)*T11
      S13=WGT( 4)*T11
      S14=WGT( 7)*T11
      S15=WGT(11)*T11
*
      T22=WGT(3)-S12**2
      IF(T22.LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T22=1./DSQRT(T22)
      S23=(WGT( 5)-S12*S13)*T22
      S24=(WGT( 8)-S12*S14)*T22
      S25=(WGT(12)-S12*S15)*T22
*
      T33=WGT(6)-S13**2-S23**2
      IF(T33.LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T33=1./DSQRT(T33)
      S34=(WGT( 9)-S13*S14-S23*S24)*T33
      S35=(WGT(13)-S13*S15-S23*S25)*T33
*
      T44=WGT(10)-S14**2-S24**2-S34**2
      IF(T44.LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T44=1./DSQRT(T44)
      S45=(WGT(14)-S14*S15-S24*S25-S34*S35)*T44
*
      T55=WGT(15)-S15**2-S25**2-S35**2-S45**2
      IF(T55.LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T55=1./DSQRT(T55)
*
      T45=-T44*(S45*T55)
      T34=-T33*(S34*T44)
      T35=-T33*(S34*T45+S35*T55)
      T23=-T22*(S23*T33)
      T24=-T22*(S23*T34+S24*T44)
      T25=-T22*(S23*T35+S24*T45+S25*T55)
      T12=-T11*(S12*T22)
      T13=-T11*(S12*T23+S13*T33)
      T14=-T11*(S12*T24+S13*T34+S14*T44)
      T15=-T11*(S12*T25+S13*T35+S14*T45+S15*T55)
*
      COV( 1)=T11*T11+T12*T12+T13*T13+T14*T14+T15*T15
      COV( 2)=T12*T22+T13*T23+T14*T24+T15*T25
      COV( 3)=T22*T22+T23*T23+T24*T24+T25*T25
      COV( 4)=T13*T33+T14*T34+T15*T35
      COV( 5)=T23*T33+T24*T34+T25*T35
      COV( 6)=T33*T33+T34*T34+T35*T35
      COV( 7)=T14*T44+T15*T45
      COV( 8)=T24*T44+T25*T45
      COV( 9)=T34*T44+T35*T45
      COV(10)=T44*T44+T45*T45
      COV(11)=T15*T55
      COV(12)=T25*T55
      COV(13)=T35*T55
      COV(14)=T45*T55
      COV(15)=T55*T55
*
      END
