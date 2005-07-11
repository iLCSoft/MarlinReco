*  
* $Id: fsfic2.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE FSFIC2(A,AINV,IER)
*******************************************************************************
*
*
*
*     invert 2x2 packed triangular matrix
*
*
*
*******************************************************************************
*
      IMPLICIT NONE
      DOUBLE PRECISION A(3), AINV(3), DINV, D
      INTEGER IER
      IER     = 0
      DINV    = (A(2)*A(2)-A(1)*A(3))
      IF (DINV .EQ. 0) THEN
        IER = 1
        RETURN
      ENDIF
      D       =  1./DINV
      AINV(1) = -A(3)*D
      AINV(2) =  A(2)*D
      AINV(3) = -A(1)*D
      RETURN
      END
