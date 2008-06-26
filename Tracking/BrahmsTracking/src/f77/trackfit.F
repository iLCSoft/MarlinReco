*  
* $Id: trackfit.F,v 1.4 2008-06-26 09:35:56 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.3  2007/09/05 09:47:29  rasp
* Updated version
*
* Revision 1.1  2006/11/02 12:04:59  rasp
*
* New fitting routine used in SiliconTracking and FullLDCTracking processors
*
* Revision 1.3  2005/07/29 15:06:17  aplin
* Initialised number of tracks found to 0 in tpcrun.
*
* Revision 1.2  2005/07/20 17:46:38  aplin
* The itedat bank has been added with read and write interface. Though
* only the write function is used (tkstor). The read will have to be
* done through a statement function.
*
* Revision 1.1.1.1  2005/07/11 15:15:48  aplin
* Initial version
*
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.2  2004/07/22 15:21:02  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:51  hvogt
* Brahms V308
*
*  
      SUBROUTINE TRACKFIT(NHITS,IDET,ITYPE,XHIT,YHIT,ZHIT,
     *                    RCOV,ZCOV,REF,
     *                    IERR,RFIT,RFITE,CHI2FIT,NDFFIT,NOUTL,
     *                    IDOUTL,FITCODE)
**********************************************************************
*
*     subroutine TRACKFIT (Interface to MarlinReco)
*
*     fits track given a list of measurement points
*
*     input:  NHITS - number of hits
*             IDET  - array of hit detector identifiers
*             ITYPE - array of hit type : 
*                     3 : hit in cyllindrical-type detector (VTX,SIT,TPC); 
*                     2 : hit in planar detector (FTD)
*             XHIT,YHIT,ZHIT - arrays of hit cartesian coordinates
*             RCOV - Rphi_reso**2 for each hit
*                    where Rphi is the spatial resolution in R-phi/R
*                    (for planar detectors Rphi_reso is supposed to be equal R_reso)
*             ZCOV - Z_reso**2 for each hit
*                    where Z_reso is the spatial resolution in z
*                    (ignored for planar detector)
*             REF    - reference track parameters in TANAGRA format.
*                      only used if REF(1).NE.0.0
*     output: IERR    - 0 if successful
*             RFIT    - fitted track parameters at reference point 
*                       as defined by TANAGRA convention
*             RFITE   - covariance matrix for fitted parameters as 
*                       defined by TANAGRA convention
*             CHI2FIT - fit chi2
*             NDFFIT  - number of degrees of freedom
*             NOUTL   - number of outliers
*             IDOUTL  - array of outliers indecies
*             FITCODE - measurement code, 
*                       0 : cartesian coordinates of the reference point
*                           are given (X,Y,Z)
*                       1 : coordinates of the reference point are given in 
*                           the form (R,R*phi,Z)     
*
*     Alexei Raspereza, MPI Munich, September 2006
*
**********************************************************************
      IMPLICIT NONE
#include "include/padrow.inc"
#include "include/bkgmode.inc"
#include "include/tkconst.inc"
#include "include/fkparm.inc"
#include "include/fkcons.inc"
#include "include/fktkx.inc"
#include "include/tkbank.inc"
#include "include/tkffread.inc"
#include "include/fkdebug.inc"
c#include "gcbank.inc"
#include "include/tpcdims.inc"
#include "include/detnums.inc"
      REAL PI,TWOPI
#include "include/fkpipar.inc"
#include "include/fkexts.inc"
C#include "include/fkcom.inc"
C.....max number of hits per track
      INTEGER NHTMAX
      PARAMETER (NHTMAX=LTPDRO)
      INTEGER IH,NHITS
      REAL XHIT(*),YHIT(*),ZHIT(*)
      INTEGER ITYPE(*),IDET(*)
      REAL RCOV(*),ZCOV(*)
      INTEGER NOUTL
      INTEGER IDOUTL(*)
      REAL RFIT(6)
      REAL RFITE(15)
      REAL CHI2FIT
      INTEGER NDFFIT
      INTEGER K,ID,IERR1,IERR2,ITSREF(MXTS)
      INTEGER IDDET,N,N2,IERR,SUBMOD,INDEX
      INTEGER I,J,IMPNT,NOUT,IDOU(NHTMAX),IFLTE(NHTMAX),IERX(NEXSMX)
      INTEGER DELDET,NTRY,TKNUMB,ILOC,IDTE(NTEMX)
      REAL    POS(3),COV(3),PHI,REF(*),RTSREF(MXTS),PARAM(6)
      INTEGER ITELOC(MXTE,NHTMAX),ITKLOC(MXTK),ITEFIT(MXTE)
      REAL    RTELOC(MXTE,NHTMAX),RTKLOC(MXTK),RTEFIT(MXTE)
      EQUIVALENCE (ITELOC,RTELOC)
      EQUIVALENCE (ITKLOC,RTKLOC)
      EQUIVALENCE (ITEFIT,RTEFIT)
      EQUIVALENCE (ITSREF,RTSREF)
      REAL    HIT(MXHT)
      INTEGER IHIT(MXHT)
      EQUIVALENCE (HIT,IHIT)
      INTEGER ILIST(NHTMAX),LIST3(NHTMAX)
      REAL    RLIST(NHTMAX),ZLIST(NHTMAX)
      LOGICAL FIRST
      DATA    FIRST /.TRUE./
      SAVE    FIRST
      REAL RHMIN
      real    rmax,rmin,zmax,zmin,pmax,pmin
      integer hid
      data    hid /100000/
      save    hid

      INTEGER IXPNT
      REAL DISTMIN
      REAL DELTADIST
      REAL REFRAD
      INTEGER IX
      INTEGER FITCODE
      INTEGER IMTRK
      PARAMETER (IMTRK=0      ! cylinder coordinates bit will be added when used
     >                +4+8+16 ! 3 coordinates of reference point measured
     >                 +32+64 ! covariance matrix for r and rphi
     >               +128+256 ! track angles theta, phi measured
     >                   +512)! 1/p measured

#include "include/tkfunc1.inc"
#include "include/tkfunc.inc"
#include "include/fkddes.inc"

c      write(*,*) "hello I'm about to run trackfit @$@#$#@$@#$@#$@#$"

      N = NHITS

C.....check range
      IERR=0
      IF (N.GT.NHTMAX) THEN
         WRITE(6,*) 'TKMKTE: too many hits. increase NHTMAX to ',N
         N=NHTMAX
      END IF
      IF (N.LT.3) THEN
         WRITE(6,*) 'TKMKTE: too few hits'
 
         IERR=-2
         GOTO 999
      END IF


C.....sort hits according to increasing r
      J=0

C      WRITE(*,*) 'NUMBER OF HITS = ',N
      DO I=1,N
         J=J+1
C         WRITE(*,*) XHIT(I),' ',YHIT(I),' ',ZHIT(I),' ',ITYPE(I)
         RLIST(J)=SQRT(XHIT(I)**2+YHIT(I)**2)
         ZLIST(J)=ZHIT(I)
      END DO
      N=J
      CALL SORTZV(RLIST,ILIST,N,1,0,0)


c$$$      WRITE(*,*) ' NCMAT = ', NCMAT
c$$$      DO I=1,NCMAT
c$$$        WRITE(*,*) I,'XRLC = ',XRLC(I)
c$$$      ENDDO

C.....fill all associated hits into local TE structure for TK2FIT
      ILOC=0
      DO I=1,N

C.......read hit
         IH = ILIST(I)
         HIT(1) = XHIT(IH)
         HIT(2) = YHIT(IH)
         HIT(3) = ZHIT(IH)
         LIST3(I) = IH

C.......calculate position and covariance matrix for desired detector type

         RHMIN = 1.0E+10

        IF (ITYPE(IH).EQ.3) THEN
           
          ! cylindrical or conical detector. convert hit pos to cyl. coordinates
          POS(1)=SQRT(HIT(1)**2+HIT(2)**2)
          POS(2)=ATAN2(HIT(2),HIT(1))
          IF (POS(1).LE.RHMIN) RHMIN = POS(1)          
          IF (POS(2).LT.0.0) POS(2)=POS(2)+2.*PI
          POS(2)=POS(2)*POS(1)
          POS(3)=HIT(3)
          
          COV(1)=RCOV(IH)
          COV(2)=0.0
          COV(3)=ZCOV(IH)
          IMPNT=1     ! cylinder coordinates
     >        +4+8+16 ! 3 coordinates of reference point measured
     >        +32+64  ! covariance matrix for r and rphi
c          WRITE(*,23) ILOC,IDET(IH),HIT(1),HIT(2),HIT(3),
c     >                              COV(1),COV(3)
        ELSE
          ! forward detector (planar)
          POS(1)=SQRT(HIT(1)**2+HIT(2)**2)
          IF (POS(1).LE.RHMIN) RHMIN = POS(1)          
          POS(1)=HIT(1)
          POS(2)=HIT(2)
          POS(3)=HIT(3)
          
          PHI=ATAN2(POS(2),POS(1))
          COV(1)= RCOV(IH)
          COV(2)= 0.0
          COV(3)= RCOV(IH)
          IMPNT= 0     ! cartesian coordinates
     >         +4+8+16 ! 3 coordinates of reference point measured
     >         +32+64  ! covariance matrix for x and y
c          WRITE(*,24) ILOC,IDET(IH),HIT(1),HIT(2),HIT(3),
c     >                              COV(1),COV(3)
        END IF
 23     FORMAT(I4,' VTX HIT : IDET = ',I2,
     >                      ' R = ',F6.2,' P = ',F6.2,' Z = ',F6.2,
     >                      ' RP = ',F10.8,'  RZ = ',F10.8)
 24     FORMAT(I4,' FTD HIT : IDET = ',I2,
     >                      ' X = ',F6.2,' Y = ',F6.2,' Z = ',F6.2,
     >                      ' RX = ',F10.8,'  RY = ',F10.8)


C.......save in local TE array
        ILOC=ILOC+1
        ITELOC( 1,ILOC)=IDET(IH) ! module identifier
        ITELOC( 2,ILOC)=0     ! submodule id
        ITELOC( 3,ILOC)=0     ! reserved; but must be zero
        ITELOC( 4,ILOC)=IMPNT ! measurement code
        ITELOC( 5,ILOC)=20    ! pointer to mass info
        ITELOC( 6,ILOC)=3     ! unknown charge
        ITELOC( 7,ILOC)=0     ! no degrees of freedom
        RTELOC( 8,ILOC)=0.0   ! chi^2 of the fit
        RTELOC( 9,ILOC)=0.1   ! length of track element
        RTELOC(10,ILOC)=POS(1) ! r or x
        RTELOC(11,ILOC)=POS(2) ! r*phi or y
        RTELOC(12,ILOC)=POS(3) ! z
        RTELOC(13,ILOC)=0.0   ! theta angle
        RTELOC(14,ILOC)=0.0   ! phi angle
        RTELOC(15,ILOC)=0.0   ! 1/p at reference point
        RTELOC(16,ILOC)=0.0   ! dE/dx or something else
        RTELOC(17,ILOC)=COV(1) ! cov. matrix
        RTELOC(18,ILOC)=COV(2) ! cov. matrix
        RTELOC(19,ILOC)=COV(3) ! cov. matrix
        RTELOC(20,ILOC)=0.0   ! mass information

      END DO

C.....reference trajectory
      
 100  IF (REF(1).EQ.0.0) THEN
C       create one using FKREF
CSJAWRT
c        write(*,*) "hello I'm about to call fkref @$@#$#@$@#$@#$@#$"
        CALL FKREF(ILOC,RTELOC,RTSREF,ITSREF)
      ELSE
C       use parameters given by calling function
        CALL VZERO(ITSREF,MXTS)
        ITSREF(3)=1+4+8
        CALL UCOPY(REF,RTSREF(17),6)
      END IF

CSAJWRT      
c      do i=1, MXTS
c      write(*,*) "RTSREF(",i,") = ", RTSREF(i) 
c      end do 

CSAJWRT
c      do i=1, MXTE
c      write(*,*) "RTELOC(",i,1,") = ", RTELOC(i,1) 
c      end do 


C.....track fit
      IF (ILOC.GT.NTEMX) THEN
         IF (TKSTDBG.GT.0) THEN
            WRITE(6,*) 'TKMKTE: INTERNAL ERROR. NTEMX=',NTEMX,
     >                 ' MUST BE SET TO AT LEAST ',ILOC
         END IF
         ILOC=NTEMX
      END IF
      DO I=1,ILOC
        IDTE(I)=I
      END DO
      CALL FK3TRK(3,ILOC,IDTE,RTELOC,RTSREF,RTKLOC,
     >            NOUT,IDOU,IERR1,IERR2,IFLTE,IERX)
      
C      write(*,*) 'HERE WE ARE ---',NOUT

      IERR=IERR1+10*IERR2

C.....debugging: compare reference trajectory, fit result and true parameters
cccIF (REF(1).NE.0.0.AND.TKSTDBG.GT.1) THEN
ccc  CALL VZERO(PARAM,6)
ccc  IF (TKMCTR(1,IHIT(6))**2+TKMCTR(2,IHIT(6))**2.GT.0.0) THEN
ccc    PARAM(4)=ATAN2(SQRT(TKMCTR(1,IHIT(6))**2
ccc         +TKMCTR(2,IHIT(6))**2),TKMCTR(3,IHIT(6)))
ccc    IF (PARAM(4).LT.0.0) PARAM(4)=PARAM(4)+2.*PI
ccc    PARAM(5)=ATAN2(TKMCTR(2,IHIT(6)),TKMCTR(1,IHIT(6)))
ccc    IF (PARAM(5).LT.0.0) PARAM(5)=PARAM(5)+2.*PI
ccc    PARAM(6)=1./SQRT(TKMCTR(1,IHIT(6))**2+TKMCTR(2,IHIT(6))**2
ccc         +TKMCTR(3,IHIT(6))**2)
ccc  END IF
ccc  WRITE(6,'(A,6F10.4)') 'TK: true params : ',(PARAM(I),I=1,6)
ccc  WRITE(6,'(A,6F10.4)') 'TK: ref traj    : ',(RTSREF(I),I=17,22)
ccc  IF (IERR.EQ.0) THEN
ccc    WRITE(6,'(A,6F10.4)')'TK: fit result  : ',(RTKLOC(I),I=17,22)
ccc      WRITE(6,'(A,6F10.4)') 'TK: proj result : ',(PARAMX(I,2),I=1,6)
ccc      WRITE(6,'(A,15F10.4)') 'TK: proj result : ',(COVX(I,2),I=7,21)
ccc  ELSE
ccc    WRITE(6,'(A,6F10.4)') 'TK: error codes : ',IERR1,IERR2
ccc  END IF
ccc  WRITE(6,'(A)') 'TK: ------------------------------------------'
cccEND IF

C.....checks
      IF (IERR.NE.0) THEN
        IF (TKSTDBG.GT.0) THEN
          WRITE(6,*) 'TKMKTE: track fit failed with error ',IERR
        ENDIF
C----------------------
c        zmin=999.0
c        zmax=-999.0
c        rmin=999.0
c        rmax=-999.0
c        pmin=999.0
c        pmax=-999.0
c        do i=1,n
c          if (rteloc(10,i).lt.rmin) rmin=rteloc(10,i)
c          if (rteloc(10,i).gt.rmax) rmax=rteloc(10,i)
c          if (rteloc(12,i).lt.zmin) zmin=rteloc(12,i)
c          if (rteloc(12,i).gt.zmax) zmax=rteloc(12,i)
c          if (rteloc(11,i)/rteloc(10,i).lt.pmin)
c     >         pmin=rteloc(11,i)/rteloc(10,i)
c          if (rteloc(11,i)/rteloc(10,i).gt.pmax)
c     >         pmax=rteloc(11,i)/rteloc(10,i)
c        end do
c        call hbook2(hid,'r vs. phi',100,pmin,pmax,100,rmin,rmax,0.0)
c        call hbook2(hid+1,'r vs. z',100,zmin,zmax,100,rmin,rmax,0.0)
c        do i=1,n
c          call hf2(hid,rteloc(11,i)/rteloc(10,i),rteloc(10,i),1.0)
c          call hf2(hid+1,rteloc(12,i),rteloc(10,i),1.0)
c        end do
c        hid=hid+10
C----------------------
        ILOC=NINT(0.7*ILOC)
        IF ((ILOC.GT.N/3).AND.(ILOC.GE.3)) THEN
ccc           IF (TKSTDBG.GT.0) WRITE(6,*) 'TKMKTE: fit failed. retry.'
           GOTO 100
        END IF
ccc           WRITE(6,*) 'TKMKTE: fit failed. no retry.'
        IF (TKSTDBG.GT.0) THEN
ccc           WRITE(6,*) 'TKMKTE: fit failed. no retry.'
ccc           WRITE(6,*) '  track properties of last hit: '
ccc           PMAX=SQRT(TKMCTR(1,IHIT(6))**2+TKMCTR(2,IHIT(6))**2
ccc     >          +TKMCTR(3,IHIT(6))**2)
ccc           WRITE(6,*) '    momentum: ',PMAX
ccc           IF (PMAX.GT.PFPCUT) WRITE(6,*) '    TKMKTE: good p failed'
ccc           WRITE(6,*) '    theta   : ',ATAN2(SQRT(TKMCTR(1,IHIT(6))**2
ccc     >          +TKMCTR(2,IHIT(6))**2),TKMCTR(3,IHIT(6)))
        END IF
        GOTO 999
      END IF

      IF (RTKLOC(17).EQ.0.0) THEN
        STOP 'TKMKTE: reference point at r=0.0 from FK3TRK. Fail.'
      END IF

      if (RTKLOC(20).lt.0. .or. RTKLOC(20).gt.pi .or.
     1    RTKLOC(21).lt.0. .or. RTKLOC(21).gt.2*pi .or.
     1    RTKLOC(10).lt.0.) then
        if (TKSTDBG.gt.0) then
          print *,'WARNING in TKMKTE: fitted theta/phi/length is off'
          print *,'TS dump: ',(ITKLOC(k),k=1,8),(RTKLOC(k),k=9,22)
          print *,'Drop it !'
        endif
        GOTO 999
      endif

      IF (NOUT.GT.NHMAX) THEN
        WRITE(6,*) 'TKMKTE: IDOU ARRAY BOUNDARY VIOLATION. CRASH.'
        STOP
      END IF

      IF (NOUT.GT.0.AND.TKSTDBG.GT.0) THEN
        WRITE(6,*) 'TKMKTE: WARNING - track fit rejected ',NOUT,
     >             ' out of ',N,' hits.'
      END IF

C.....get submodule code
c$$$      SUBMOD=0
c$$$      DO I=1,N
c$$$        IF (IFLTE(I).EQ.0) THEN
c$$$          IF (LIST2(I)-IDDET.GE.1) THEN
c$$$            SUBMOD=IOR(SUBMOD,2**(LIST2(I)-IDDET-1))
c$$$          ELSE
c$$$            SUBMOD=1
c$$$          END IF
c$$$        END IF
c$$$      END DO

C.....create TE structure
      ITEFIT( 1)=IDDET          ! subdetector ID
      ITEFIT( 2)=SUBMOD         ! submodule code
      ITEFIT( 3)=0              ! reserved; must be zero
      ITEFIT( 4)=IMTRK+IAND(ITKLOC(3),1) ! data descriptor/measurement code
      ITEFIT( 5)=32             ! pointer to the mass information
      ITEFIT( 6)=ITKLOC(6)      ! charge
      ITEFIT( 7)=ITKLOC(8)      ! number of degrees of freedom
      RTEFIT( 8)=RTKLOC(9)      ! chi^2
      RTEFIT( 9)=RTKLOC(10)     ! length of track element
      RTEFIT(10)=RTKLOC(17)     ! x or R     )
      RTEFIT(11)=RTKLOC(18)     ! y or R*phi ) of the TE reference point
      RTEFIT(12)=RTKLOC(19)     ! z          )
      RTEFIT(13)=RTKLOC(20)     ! theta angle
      RTEFIT(14)=RTKLOC(21)     ! phi angle
      RTEFIT(15)=RTKLOC(22)     ! 1/p at reference point
      RTEFIT(16)=    0.0        ! dE/dx measurement (not yet...)
      ! covariance matrix
      CALL UCOPY(RTKLOC(23),RTEFIT(17),15)
      CALL UCOPY(RTKLOC(17),RFIT(1),6)
      CALL UCOPY(RTKLOC(23),RFITE(1),15)
      CHI2FIT = RTKLOC(9)
      NDFFIT  = ITKLOC(8)
      FITCODE = ITKLOC(3)

      NOUTL = NOUT
C      WRITE(*,*) 'NOUT = ',NOUT
      
      DO I=1,NOUTL
        IDOUTL(I)=LIST3(IDOU(I))
C        WRITE(*,*) IDOUTL(I)
      END DO

      

      IXPNT = 2
      DISTMIN = 1E+10
      REFRAD = SQRT(REF(1)**2+REF(2)**2)
      DELTADIST = 0

      DO 63 IX=1,NEXS
         IF (IFLX(IX).EQ.0) THEN
            IF (TYPX(IX).EQ.'CYLI'.AND.ABS(DATX(1,IX)).GT.REFRAD) THEN               
               DELTADIST = ABS(DATX(1,IX)-REFRAD)
c               DELTADIST = ABS(DATX(1,IX)-RHMIN)  
c            ELSE
c               DELTADIST = ABS(DATX(3,IX)-REF(3))
c               REFRAD = SQRT(DATX(1,IX)**2+DATX(2,IX))
c               DELTADIST = ABS(REFRAD-RHMIN)
            ENDIF
            IF (DELTADIST.LT.DISTMIN) THEN
               DISTMIN = DELTADIST
               IXPNT = IX
            ENDIF
         ENDIF
 63   CONTINUE

      IF (TYPX(IXPNT).EQ.'CYLI') THEN
         FITCODE = 1
      ELSE
         FITCODE = 0
      ENDIF

      DO IX=1,6
         RFIT(IX) = DATX(IX,IXPNT)
      ENDDO

      DO IX=1,15
         RFITE(IX) = DATX(IX+6,IXPNT)
      ENDDO


c      WRITE(*,*) RTKLOC(17),RTKLOC(18),RTKLOC(19)
c      WRITE(*,*) REF(1),REF(2),REF(3)

C.....update hit lists
      N2=0
c$$$      DO I=1,N
c$$$        IF (IFLTE(I).EQ.0) THEN
c$$$          N2=N2+1
c$$$          LIST1(N2)=LIST1(I)
c$$$          LIST2(N2)=LIST2(I)
c$$$        END IF
c$$$      END DO

C.....store TE
C      CALL TKSTOR('TE',INDEX,ITEFIT,N2,LIST1,LIST2,IERR)



      IF (IERR.NE.0) THEN
        WRITE(6,*) 'TKMKTE: TKSTOR FAILED.'
        CALL TKVETO(1)
        GOTO 999
      END IF

C.....debug information
      IF (TKSTDBG.GT.1) THEN
        WRITE(6,*) 'TKMKTE: last TE stored at index ',INDEX
      END IF

C.....before returning, fill histograms with tracks
 999  IF (TKSTHST.GT.0) THEN
c        CALL HCDIR('//PAWC/TRACKING/TKSTEER',' ')
        IF (FIRST) THEN
          FIRST=.FALSE.
c          CALL HBOOK2(1,'good tracks',100,-160.,160.,100,-160.,160.,0.)
c          CALL HBOOK2(2,'bad tracks', 100,-160.,160.,100,-160.,160.,0.)
        END IF
        IF (IERR.EQ.0) THEN
          ID=1
        ELSE
          ID=2
        END IF
        DO I=1,ILOC
c          CALL HF2(ID,RTELOC(10,I)*COS(RTELOC(11,I)/RTELOC(10,I)),
c     >         RTELOC(10,I)*SIN(RTELOC(11,I)/RTELOC(10,I)),1.0)
        END DO
      END IF

      call flush
      RETURN
      END


