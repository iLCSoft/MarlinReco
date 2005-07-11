CDECK  ID>, SYISRT.
********************************
      SUBROUTINE SYISRT(N,Q,III)
*  Just sort III with energy
*  S.Y
********************************
      implicit none
      INTEGER N,III(*),I,J
      INTEGER IMAXS
      REAL Q(5,*),PP(500)
      REAL RMAXS

      DO I=1,N
         PP(I)=Q(4,I)
      ENDDO
      DO I=1,N
         IMAXS=1
         RMAXS=PP(1)
         DO J=2,N
            IF(PP(J).GT.RMAXS) THEN
               RMAXS = PP(J)
               IMAXS = J
            ENDIF
         ENDDO
         III(I)=IMAXS
         PP(IMAXS)=-1.E+09
      ENDDO

      END
