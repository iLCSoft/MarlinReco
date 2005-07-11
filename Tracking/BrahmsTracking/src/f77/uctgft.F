*  
* $Id: uctgft.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
        SUBROUTINE UCTGFT
        IMPLICIT NONE
C****************************************************************************
c* Version number: 0.0      Date of last modification  26.08.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
C* Main modifiers from ALEPH to BRAHMS code:                                *
c* Name          Institute                        Contact email             *
c* G.A. Blair    Royal Holloway, Univ. of London  blair@ppu1.ph.rhbnc.ac.uk *
c*                                                                          *
c****************************************************************************
C----------------------------------------------------------------
C! Change to new helix convention in bank T1FT
C
C   Subroutine that changes the T1FT bank from the TASSO
C   convention of parameters to the mathematical convention,
C   and puts in the covariance matrix instead of the inverse of it.
C   This is done by refitting the track with UFITMS, which gives
C   its output in the desired convention.
C
C   AUTHOR:  LLUIS GARRIDO   14/6/1988
C   Modified: R. Johnson     19/5/1991   (use T1FT instead of TGFT)
C             D. Casper      09/7/1997   Check for garbage fits
C             D. Casper      28/1/1998   Test error return of UFITMS
C
C   Called by TFITTK
C
C   CONVERSION:
C          R_TASSO=-R_MAT
C         D0_TASSO= SIGN(R_MAT)*D0_MAT
C----------------------------------------------------------------
C
C
C----------------------------------------------------------------
C My Stuff:
#include "include/fndparms.inc"
#include "include/tpco.inc"
#include "include/my_tcal.inc"
#include "include/my_t1ft.inc"
#include "include/my_t1tl.inc"
#include "include/my_t1cl.inc"
#include "include/tpcdims.inc"
#include "include/tpcmatpar.inc"

       INTEGER J,NTR,NOFF,NTPC,K,LIST2(2),LIST3(2),IERR,IERR2
       INTEGER LIST1(200),NMULT

       REAL CURSAVE,D0SAVE,CHI2
C----------------------------------------------------------------
      REAL VV0,EE0,COV
      DIMENSION VV0(6),EE0(21),COV(21)
      LOGICAL FIRST

      REAL PROB
      EXTERNAL PROB
       SAVE
      DATA FIRST/.TRUE./
C
C#include "bmacro.h"
C
C       IF (FIRST) THEN
C         FIRST=.FALSE.
C         NT1FT=NAMIND('T1FT')
C         NT1TL=NAMIND('T1TL')
C         NT1CL=NAMIND('T1CL')
C       ENDIF
C
C       KT1FT=IW(NT1FT)
c GAB:
C       IF (KT1FT.EQ.0) RETURN
C       NTR=LROWS(KT1FT)
       IF (NT1FT.EQ.0) RETURN
       NTR=NT1FT

       IF(NTR.EQ.0) RETURN
C GAB:
C       KT1TL=IW(NT1TL)
C       KT1CL=IW(NT1CL)
C       IF (KT1TL.EQ.0 .OR. KT1CL.EQ.0) THEN
C         CALL RERROR('UCTGFT',1,'Missing T1TL and T1CL banks!')
C         CALL BDROP(IW,'T1FTT1TLT1CL')
C         RETURN
C       ENDIF

       DO 1 J=1,NTR
C GAB:
C         IT1FT=KROW(KT1FT,J)
C         CURSAVE = RW(IT1FT+JT1FIR)
C         D0SAVE  = RW(IT1FT+JT1FD0)
C         NOFF=ITABL(KT1TL,J,JT1TIO)
C         NTPC=ITABL(KT1TL,J,JT1TN1)
         CURSAVE = MY_JT1FIR(J)
         D0SAVE  = MY_JT1FD0(J)
         NOFF=MY_JT1TIO(J)
         NTPC=MY_JT1TN1(J)
c GAB:
         DO K=1,NTPC
          LIST1(K)=MY_T1CLT(NOFF+K)
         ENDDO

C         CALL UFITMS(NTPC,0,0,IW(KT1CL+LMHLEN+NOFF+1),LIST2,LIST3,2,1
C     &        ,29.6,.027,VV0,EE0,COV,CHI2,IERR)
c
         NMULT=1

         CALL UFITMS(NTPC,0,0,LIST1,LIST2,LIST3,2,NMULT
     &        ,RMS,RADLEN,VV0,EE0,COV,CHI2,IERR)

C  Check for garbage
         CALL TFCHEK(VV0,.FALSE.,IERR2)

         IF (IERR.NE.0 .OR. IERR2.NE.0) THEN
C  Do a naive conversion and hope for the best...
            VV0(1) = -CURSAVE
            VV0(4) = SIGN(1.,VV0(1))*D0SAVE
         ENDIF
C GAB:
C         DO 41 K=JT1FIR,JT1FZ0
C           RW(IT1FT+K)=VV0(K)
C   41    CONTINUE
C         DO 42 K=1,15
C           RW(IT1FT+JT1FEM-1+K)=COV(K)
C   42    CONTINUE
C         RW(IT1FT+JT1FCD)=CHI2
C         IW(IT1FT+JT1FDF)=2*NTPC-5

         MY_JT1FIR(J)=VV0(1)
         MY_JT1FTL(J)=VV0(2)
         MY_JT1FP0(J)=VV0(3)
         MY_JT1FD0(J)=VV0(4)
         MY_JT1FZ0(J)=VV0(5)
   41    CONTINUE
         DO 42 K=1,15
           MY_JT1FEM(J,K)=COV(K)
   42    CONTINUE
         MY_JT1FCD(J)=CHI2
         MY_JT1FDF(J)=2*NTPC-5

    1  CONTINUE


       RETURN
       END

C------------------------------------------------------------------------
