CDECK  ID>, SYCONE.
**************************************************************************
      SUBROUTINE SYCONE(NJREQ,NTRAK,ITKMAX,PTRAK,R,EPS,NJMAX,NJET,PJET,
     +                  IPASS,IJMUL,IERR)
*.
*...SYCONE Like PXCONE, but forcing NJET=NJREQ.
*          Fixes Epsilon and looks for largest value of R that gives NJREQ.
*.
*. INPUT     :
*. OUTPUT    :
*.
*. COMMON    :
*. SEQUENCE  :
*. CALLS     : PXCONE VZERO
*. CALLED    :
*.
*. REPORT CONDITIONS
*.
*. AUTHOR    : D.R.Ward
*. VERSION   : 1.02
*. CREATED   :  5-Aug-94
*. LAST MOD  : 09-Feb-1996
*.
*. Modification Log.
*. 09-Feb-96  S.Yamashita    get from WWCONE in WWpackage WW106
*.**********************************************************************
      implicit none
      INTEGER NTRAK,ITKMAX,NJMAX,NJET,IPASS(*),IJMUL(*),IERR,NIT,I,NJREQ
      REAL    PTRAK(ITKMAX,*),R,EPS,PJET(ITKMAX,*),R0,R1,R2,ESUM
      ESUM=0.
      DO 10 I=1,NTRAK
         ESUM=ESUM+PTRAK(4,I)
  10  CONTINUE
      R0=0.
      R1=0.8
      R2=1.6
      NIT=0
      R=R1
      EPS=ESUM/17.4
  20  CALL PXCONE(NTRAK,ITKMAX,PTRAK,R,EPS,NJMAX,NJET,PJET,
     +            IPASS,IJMUL,IERR)
*      Print *,'NIT=',nit,' R=',R,' NJET=',njet
      IF(NJET.GT.NJREQ-1) THEN
         R0=R1
         R1=(R1+R2)/2.
      ELSE
         R2=R1
         R1=(R1+R0)/2.
      ENDIF
      NIT=NIT+1
      R=R1
      IF(NIT.LE.5 .OR. (NIT.LE.10 .AND. NJET.NE.NJREQ)) GO TO 20
      R=R0
      CALL VZERO(IPASS,NTRAK)
      CALL PXCONE(NTRAK,ITKMAX,PTRAK,R,EPS,NJMAX,NJET,PJET,
     +            IPASS,IJMUL,IERR)
      IERR=0
      IF(NJET.NE.NJREQ) IERR=1
*      Print *,'Finally: R=',R,'  Njet=',Njet, ' Ierr=',ierr
      END

*********************** E N D   O F    C O D E ********************************
*
*+PATCH,SYJDOC.
*INTRODUCTION:
*
*  SYJJET
*   is a set of routine which is desined to get better association of
*   jet and particles especially for high mass Higgs and WW/ZZ SM processes.
*
*   It needs to link with PX-lib, CKERN-lib (Stan's code for Cambridge)
*
*    1. SOURCE CODE    :  /u/ws/satoru/SYJJET/syjjet.car   (SnOPAL)
*                         /u/ws/satoru/SYJJET/syjjet.car   (shift)
*    2. Cradle File    :  /u/ws/satoru/SYJJET/SYJJET.CRA   (SnOPAL)
*                         /u/ws/satoru/SYJJET/SYJJET.CRA   (shift)
*    3. Library        :  /u/ws/satoru/syjjet.a     (SnOPAL)
*                         /u/ws/satoru/syjjet.a     (shift)
*  ** This library already includes CKERN. So you just need to link with
*   PX library as usual.
*
*
*   This program can handle so-called traditional jet-finders and also
*   newly developed "jet-associator".
*
*   The idea of the improvement of jet-particle association is as follows.
*   1. Make core jet
*   2. using core-jet direction etc... re-association is done for all
*      particles choosing the closest jet in terms of angle, Jade-parameter
*      etc...
*   Hence there are variety of combination can be made for
*   possible core formation, jet-association, jet-merge etc...
*
*   I've chosen a sort of good association and can be done simplly with
*   call one routine
*  SYJJET4 for 4 jet formation.
*
*   This package includes a kernel program
*  SYJKRN
*
*   This routine can handle various jet finders and re-association scheme.
*
*    Currently available traditional jet-finders are
*   1. JADE E0,E,P,P0
*   2. Geneva
*   3. Durham
*   4. Cambridge
*   5. Lund
*   6. Cone
*
*    All can be used for
*   A) fixed Ycut(Xcut/R/Eps)
*   B) fixed number of jets
*
*    Reassociation can be done with various "distance-parameters" such as
*    1. angle
*    2. JADE-type
*    3. Invariant mass
*    4. Durham-type
*        etc...
*
*    Also Jet Merge can be done with similar paramter for jet-jet "distance".
*
*    5jet, variable n-jet, etc... almost all things can be done with this
*    program.
*
* ============================================================================
*USAGE:
*
*   **** 4 J E T *****
* Just for 4-Jet, I describe here how to call it.
*
* O. link:
*      your program and this library (syjjet.a) and px-library.
*
* A. CALL SYJJET4(NPAR,PPAR,JAS,PJET,Y34,Y45,IERR)
*
**  Inputs:
**    NPAR  :       I : Number of "MT-particles"
**    PPAR(5,*)     R : 5-momentum for each "MT-particles"
**  Outputs:
**    JAS(*)        I : Jet association like YASSO
**    PJET(5,*)     R : Jet 5-momentum
**    Y34 Y45       R : Durham Y34/Y45
**    IERR          I : 0 for O.K.
*
*   If all calculations are O.K. (normal) IERR is 0.
*   Y34 and Y45 with normal Durham are automatically back as Y34 and Y45.
*   JAS(i) (i=1 to NPAR) should be between 1 and 4.
*   This JAS is the results of re-association.
*   PJET is the jet energy with energy ordering.
*
*  The internal procedure for the jet-finding are as follows.
*  1. The program starts from 4-core-jets formation just with "particle" with
*   more than 1.2 GeV threshold.
*  2. The event is forced to be 4-jets.
*  3. Then re-association starts.
*   The "distance" used in this version is
*   JADE-E0 type parameter between "particle" and "core-jet".
*   All "particles" including that not used in core-formation (having lower
*   energy than 1.2 GeV) are re-assigned to jet having smallest
