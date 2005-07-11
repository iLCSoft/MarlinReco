*  
* $Id: smav.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      SUBROUTINE SMAV(R,A,V,N)
C -----------------------------------------------------------
C!    R(N) = A(N,N)*V(N)
C!    AUTHOR: MARTIN POPPE
C -----------------------------------------------------------
      SAVE
C
      DIMENSION R(*),A(*),V(*)
C---> ADDRESSE IN TRIAGULAR MATRIX, ROW II COLUMN KK
      IND(II,KK)=(II*II-II)/2+KK
C---> ADDRESSE IN TRIAGULAR MATRIX, ROW II COLUMN KK
      DO 30 I=1,N
        R(I) = 0.0
        DO 20 K=1,N
          IF(I.GE.K)THEN
            R(I)=R(I)+A(IND(I,K))*V(K)
          ELSE
            R(I)=R(I)+A(IND(K,I))*V(K)
          END IF
  20    CONTINUE
  30  CONTINUE
      RETURN
      END

C==================================================================
