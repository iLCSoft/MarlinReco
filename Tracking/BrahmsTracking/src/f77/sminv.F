*  
* $Id: sminv.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      SUBROUTINE SMINV(V,B,N,M,NRANK)
C     AUTHOR: V.BLOBEL
C
C     SUBROUTINE SMINV
C     ----------------
C!    OBTAIN SOLUTION OF A SYSTEM OF LINEAR EQUATIONS V *  X  =  B  WITH
C!    SYMMETRIC MATRIX V AND INVERSE (FOR M =  1)  OR  MATRIX  INVERSION
C!    ONLY (FOR M = 0)
C
C                   - - - -
C        CALL SMINV(V,B,N,M,NRANK)
C                   - -     -----
C
C           V = SYMMETRIC N-BY-N MATRIX IN SYMMETRIC STORAGE MODE
C               V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, . . .
C               REPLACED BY INVERSE MATRIX
C           B = N-VECTOR   (FOR M = 0 USE A DUMMY ARGUMENT)
C               REPLACED BY SOLUTION VECTOR
C           M = SEE ABOVE
C
C
C     METHOD OF SOLUTION IS BY ELIMINATION SELECTING THE  PIVOT  ON  THE
C     DIAGONAL EACH STAGE. THE RANK OF THE MATRIX IS RETURNED IN  NRANK.
C     FOR NRANK NE N, ALL REMAINING  ROWS  AND  COLS  OF  THE  RESULTING
C     MATRIX V AND THE CORRESPONDING ELEMENTS OF  B  ARE  SET  TO  ZERO.
C     SMINV USES A WORK ARRAY OF 2*N WORDS IN COMMON/MATCOM/. FOR N> 200
C     THE USER HAS TO DEFINE COMMON/MATCOM/ WITH 2*N WORDS.
C ----------------------------------------------------------------------
      SAVE
C
      REAL V(*),B(*),VKK,D,E
      REAL EPS
      COMMON/MATCOM/DR(2,200)
      DATA EPS/1.E-6/
C
C     CONSTRUCT TABLE
C
      DO 10 I=1,N
   10 DR(1,I)=1.0
      NI=N
      GOTO 14
      ENTRY XMINV(V,B,N,M,NRANK)
      NI=0
      DO 12 I=1,N
      IF(DR(1,I).EQ.0.0) GOTO 12
      NI=NI+1
   12 CONTINUE
   14 II=0
      DO 16 I=1,N
      II=II+I
   16 DR(2,I)=ABS(V(II))
C
C     LOOP BEGIN
C
      NRANK=N-NI
      DO 60 I=1,NI
C
C     SEARCH FOR PIVOT AND TEST FOR LINEARITY AND ZERO MATRIX
C
      K=0
      JJ=0
      VKK=0.0
      DO 20 J=1,N
      JJ=JJ+J
      IF(DR(1,J).EQ.0.0) GOTO 20
      IF(ABS(V(JJ)).LE.VKK) GOTO 20
      IF(ABS(V(JJ)).LT.EPS*DR(2,J)) GOTO 20
      VKK=ABS(V(JJ))
      K=J
      KK=JJ
   20 CONTINUE
      IF(K.EQ.0) GOTO 80
C
C     PREPARATION FOR ELIMINATION
C
      NRANK=NRANK+1
      DR(1,K)=0.0
      D=1.0/V(KK)
      V(KK)=-D
      IF(M.EQ.1) B(K)=B(K)*D
      JK=KK-K
      JL=0
C
C          ELIMINATION
C
      DO 50 J=1,N
      IF(J-K) 24,22,26
   22 JK=KK
      JL=JL+J
      GOTO 50
   24 JK=JK+1
      GOTO 28
   26 JK=JK+J-1
   28 E=V(JK)
      V(JK)=D*E
      IF(M.EQ.1) B(J)=B(J)-B(K)*E
      LK=KK-K
      DO 40 L=1,J
      JL=JL+1
      IF(L-K) 34,32,36
   32 LK=KK
      GOTO 40
   34 LK=LK+1
      GOTO 38
   36 LK=LK+L-1
   38 V(JL)=V(JL)-V(LK)*E
   40 CONTINUE
   50 CONTINUE
   60 CONTINUE
C
C          CHANGE SIGN
C
      IJ=0
      DO 70 I=1,N
      DO 70 J=1,I
      IJ=IJ+1
   70 V(IJ)=-V(IJ)
      GOTO 100
C
C          CLEAR REST OF MATRIX
C
   80 IJ=0
      DO 90 I=1,N
      IF(M.EQ.1.AND.DR(1,I).NE.0.0) B(I)=0.0
      DO 90 J=1,I
      IJ=IJ+1
      IF(DR(1,I)+DR(1,J).NE.0.0) V(IJ)=0.0
   90 V(IJ)=-V(IJ)
  100 RETURN
      END
C--------------------------------------------------------------
