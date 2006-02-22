*  
* $Id: tkaddh.F,v 1.2 2005-07-21 16:37:00 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1  2005/07/20 17:46:38  aplin
* The itedat bank has been added with read and write interface. Though
* only the write function is used (tkstor). The read will have to be
* done through a statement function.
*
* Revision 1.2  2004/07/22 15:16:03  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKADDH
**********************************************************************
*
*     subroutine TKADDH
*
*     add TEs for all unused SI detector hits to TE bank
*
*     Kristian Harder, November 1999
*
**********************************************************************
      IMPLICIT NONE
#include "include/bkgmode.inc"
#include "include/tkconst.inc"
#include "include/tkffread.inc"
#include "include/fkparm.inc"
#include "include/tkbank.inc"
C#include "gcbank.inc"
#include "include/detnums.inc"
      REAL PI,TWOPI
#include "include/fkpipar.inc"

      INTEGER I,J,IUSED(NHMAX),INDEX,IERR
      INTEGER ITELOC(MXTE),IMPNT
      REAL    RTELOC(MXTE),POS(3),COV(3),PHI
      EQUIVALENCE (ITELOC,RTELOC)

#include "include/tkfunc1.inc"
#include "include/tkfunc.inc"

C.....init marker array
      CALL VZERO(IUSED,NHMAX)

C.....loop over all TEs of this subdetector
      DO I=1,NTE
        IF (ITE(1,I).EQ.ID_VTX) THEN

C.........mark all constituent hits as used in this TE

C     this needs to be uncommented after ITEDAT included in the C++/F77 interface

           DO J=ITEDAT(1,I),ITEDAT(1,I)+ITEDAT(2,I)-1
              IUSED(IHITTE(J))=I
           END DO

        END IF
      END DO

C.....now loop over all hits and create TEs for previously unused ones
      DO I=1,NTKHIT
        IF ((ITKHIT(5,I)/100.EQ.ID_VTX/100).OR.
     >      (ITKHIT(5,I)/100.EQ.ID_FTD/100).OR.
     >      (ITKHIT(5,I)/100.EQ.ID_SIT/100)) THEN

C.........create TE if either hit is unused or we want to save all hits anyway
          IF ((IUSED(I).EQ.0).OR.(IPATRC(1).EQ.1)) THEN

C.........calculate position and covariance matrix for desired detector type
            IF (ITKHIT(9,I).EQ.3) THEN
              ! cylindrical detector. convert hit pos to cyl. coordinates
              POS(1)=SQRT(RTKHIT(1,I)**2+RTKHIT(2,I)**2)
              POS(2)=ATAN2(RTKHIT(2,I),RTKHIT(1,I))
              IF (POS(2).LT.0.0) POS(2)=POS(2)+2.*PI
              POS(2)=POS(2)*POS(1)
              POS(3)=RTKHIT(3,I)
              COV(1)=RTKHIT(10,I)**2
              COV(2)=0.0
              COV(3)=RTKHIT(11,I)**2
              IMPNT=1     ! cylinder coordinates
     >            +4+8+16 ! 3 coordinates of reference point measured
     >            +32+64  ! covariance matrix for r and rphi
            ELSE
              ! forward detector (planar)
              POS(1)=RTKHIT(1,I)
              POS(2)=RTKHIT(2,I)
              POS(3)=RTKHIT(3,I)
              PHI=ATAN2(POS(2),POS(1))
              COV(1)= COS(PHI)**2*RTKHIT(11,I)**2
     >               +SIN(PHI)**2*RTKHIT(10,I)**2
              COV(2)= SIN(PHI)*COS(PHI)
     >               *(RTKHIT(11,I)**2-RTKHIT(10,I)**2)
              COV(3)= SIN(PHI)**2*RTKHIT(11,I)**2
     >               +COS(PHI)**2*RTKHIT(10,I)**2
              IMPNT= 0     ! cartesian coordinates
     >             +4+8+16 ! 3 coordinates of reference point measured
     >             +32+64  ! covariance matrix for x and y
            END IF

C...........save in local TE array
            ITELOC( 1)=ITKHIT(5,I) ! module identifier
            IF (ITKHIT(5,I)-ITKHIT(5,I)/100*100.GT.0) THEN
              ITELOC(2)=2**(ITKHIT(5,I)-ITKHIT(5,I)/100*100-1) ! submodule ID (layered det.)
            ELSE
              ITELOC(2)=0                         ! submodule ID (monolithic)
            END IF
            ITELOC( 3)=0        ! reserved; but must be zero
            ITELOC( 4)=IMPNT    ! measurement code
            ITELOC( 5)=20       ! pointer to mass info
            ITELOC( 6)=3        ! unknown charge
            ITELOC( 7)=0        ! no degrees of freedom
            RTELOC( 8)=0.0      ! chi^2 of the fit
            RTELOC( 9)=0.1      ! length of track element
            RTELOC(10)=POS(1)   ! r or x
            RTELOC(11)=POS(2)   ! r*phi or y
            RTELOC(12)=POS(3)   ! z
            RTELOC(13)=0.0      ! theta angle
            RTELOC(14)=0.0      ! phi angle
            RTELOC(15)=0.0      ! 1/p at reference point
            RTELOC(16)=0.0      ! dE/dx or something else
            RTELOC(17)=COV(1)   ! cov. matrix
            RTELOC(18)=COV(2)   ! cov. matrix
            RTELOC(19)=COV(3)   ! cov. matrix
            RTELOC(20)=0.0      ! mass information

C...........store as TE
            CALL TKSTOR('TE',INDEX,RTELOC,1,I,0,IERR)
            IF (IERR.NE.0) THEN
              WRITE(6,*) 'TKADDH: HIT NOT STORED AS TE.'
              CALL TKVETO(1)
            END IF
          END IF
        END IF
      END DO

  999 RETURN
      END


