*  
* $Id: fkmi3d.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE FKMI3D(WGT,COV,IERR)
*                                                                     *
*   Inversion of a (3x3) symmetric positive matrix (Cholesky method)  *
*   Internal computation in double precision                          *
*   Check on positivity (IERR=2 if WGT not positive)                  *
*   COV (output) may overwrite WGT (input)                            *
*                                                                     *
***********************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION WGT(6),COV(6)
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
*
      T22=WGT(3)-S12**2
      IF(T22.LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T22=1./DSQRT(T22)
      S23=(WGT( 5)-S12*S13)*T22
*
      T33=WGT(6)-S13**2-S23**2
      IF(T33.LE.0.) THEN
        IERR=2
        RETURN
      ENDIF
      T33=1./DSQRT(T33)
*
      T23=-T22*(S23*T33)
      T12=-T11*(S12*T22)
      T13=-T11*(S12*T23+S13*T33)
*
      COV( 1)=T11*T11+T12*T12+T13*T13
      COV( 2)=T12*T22+T13*T23
      COV( 3)=T22*T22+T23*T23
      COV( 4)=T13*T33
      COV( 5)=T23*T33
      COV( 6)=T33*T33
*
      END
