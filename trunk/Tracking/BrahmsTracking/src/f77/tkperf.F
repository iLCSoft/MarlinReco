*  
* $Id: tkperf.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.2  2004/07/22 15:21:02  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:51  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKPERF
**********************************************************************
*
*     subroutine TKPERF
*
*     performance analysis at the end of each event
*
*     Kristian Harder, November 1999
*
**********************************************************************
      IMPLICIT NONE
#include "include/fsfparams.inc"
c#include "include/gcflag.inc"
c#include "include/gcbank.inc"
#include "include/fsfpi.inc"
#include "include/tkconst.inc"
#include "include/tkbank.inc"
#include "include/detnums.inc"
#include "include/tkperfc.inc"
#include "include/brcntl.inc"
#include "include/tkffread.inc"
#include "include/tkfunc1.inc"

      CHARACTER*15 USE
      INTEGER I,J,K,L,M,ID,TKNHIT,NHTREE,NHREC,IERR,IDET,SIHITS
      INTEGER ICOUNT(DETMAX,0:NTKMX),TENUM,DELID,TELAB,IDGEN
      INTEGER IHCOUNT(DETMAX,0:NTKMX)
      REAL    FAKCT,FAKP,TRUEP
      LOGICAL TKCOUNT,TKFCOUNT,LSI

      INTEGER BESTTS(0:NTKMX),BESTTK(0:NTKMX)
      INTEGER NUMTS(0:NTKMX),NUMTK(0:NTKMX)
      INTEGER LENTS(0:NTKMX),LENTK(0:NTKMX)

      INTEGER MXSORT
      PARAMETER (MXSORT=1000)
      INTEGER ISORT(MXSORT),INDEX(MXSORT)

      INTEGER MLIST
      PARAMETER (MLIST=25)
      INTEGER NLIST,ILIST(MLIST),ILAB(MLIST)
      REAL    RLIST(MLIST)

#include "include/tkfunc.inc"

C.....skip this if event was vetoed for technical reasons
      I=-1
      CALL TKVETO(I)
      IF (I.GT.0) GOTO 999

      write(*,*), "tkperf called: tkveto passed"

CSJA      IF (TKSTHST.GT.0) CALL HCDIR('//PAWC/TRACKING/TKSTEER',' ')
CSJA
CSJAC.....------------------------------------------------------------
CSJAC.....=== CALCULATE CUT VALUES FOR MC TRACKS (USED IN TKCOUNT) ===
CSJAC.....------------------------------------------------------------
CSJA
CSJAC.....calculate p, and |cos(theta)| for each MC track
CSJA      DO I=1,TKNTRK
CSJA        PFPMC(I)=SQRT(TKMCTR(1,I)**2+TKMCTR(2,I)**2+TKMCTR(3,I)**2)
CSJA
CSJAC.......the usual case
CSJA        IF (PFPMC(I).GT.0.0) THEN
CSJA          PFCTMC(I)=TKMCTR(3,I)/PFPMC(I)
CSJA
CSJAC.......secondary particles not treated correctly within BRAHMS are vetoed
CSJA        ELSE
CSJA          PFCTMC(I)=10.0
CSJA        END IF
CSJA      END DO
CSJA
CSJA
C.....--------------------------------
C.....=== LOCAL PATREC PERFORMANCE ===
C.....--------------------------------

C.....count number of hits each track produces in each subdetector
      CALL VZERO(IHCOUNT,DETMAX*(TKNTRK+1))
      DO I=1,NTKHIT
        IDET=DELID(ITKHIT(5,I))
        ! use generic detector ID if layer IDs are unknown to DELSEARCH (FCH!)
        IF (IDET.EQ.0) THEN
          IDET=DELID(ITKHIT(5,I)/100*100)
        END IF
        IHCOUNT(IDET,IABS(ITKHIT(6,I)))=
     >       IHCOUNT(IDET,IABS(ITKHIT(6,I)))+1
        ! check generic detector IDs for discrete detectors
        IDGEN=ITKHIT(5,I)/100*100
        ! patrec combines VTX and FTD and SIT into 'VTX'
        IF (IDGEN.EQ.ID_FTD) IDGEN=ID_VTX
        IF (IDGEN.EQ.ID_SIT) IDGEN=ID_VTX
        IF ((IDGEN.EQ.ID_VTX).AND.
     >      (ITKHIT(5,I)/100*100.NE.ITKHIT(5,I))) THEN
          IHCOUNT(DELID(IDGEN),IABS(ITKHIT(6,I)))=
     >         IHCOUNT(DELID(IDGEN),IABS(ITKHIT(6,I)))+1
        END IF
      END DO

CSJAC.....for each subdetector, count number of tracks with certain number of hits
CSJA      DO I=1,DETMAX
CSJA        DO J=1,TKNTRK
CSJA          IF (TKCOUNT(J)) THEN
CSJA            IF (IHCOUNT(I,J).GE.1) THEN
CSJA              PFLOPA(1,I)=PFLOPA(1,I)+1 ! this is only approximately true
CSJA              PFLOPA(2,I)=PFLOPA(2,I)+1
CSJA              IF (IHCOUNT(I,J).GE.3) THEN
CSJA                PFLOPA(3,I)=PFLOPA(3,I)+1
CSJA              END IF
CSJA            END IF
CSJA          END IF
CSJA        END DO
CSJA      END DO
CSJA
CSJAC.....count number of TEs for each MC track and each subdetector
CSJAC     count only TEs that fulfil all cuts
CSJA      CALL VZERO(ICOUNT,DETMAX*(TKNTRK+1))
CSJA      DO I=1,NTE
CSJA        IDET=DELID(ITE(1,I))
CSJA        IF (IDET.NE.0.AND.ITEDAT(5,I).NE.0) THEN
CSJA          IF (TKCOUNT(IABS(ITEDAT(5,I)))) THEN
CSJA            ICOUNT(IDET,IABS(ITEDAT(5,I)))=
CSJA     >            ICOUNT(IDET,IABS(ITEDAT(5,I)))+1
CSJA          END IF
CSJA        ELSE
CSJA          IF (TKFCOUNT(RTE(13,I),RTE(15,I),ITE(4,I))) THEN
CSJA            ICOUNT(IDET,0)=ICOUNT(IDET,0)+1
CSJA          END IF
CSJA        END IF
CSJA      END DO
CSJA
CSJAC.....count number of tracks that have valid TEs in each subdetector
CSJAC     and the number of fake tracks in each subdetector
CSJA      DO I=1,DETMAX
CSJA        ! fake TEs in this detector
CSJA        PFLOPA(9,I)=PFLOPA(9,I)+ICOUNT(I,0)
CSJA        DO J=1,TKNTRK
CSJA          ! reconstructed TEs in this detector
CSJA          IF (ICOUNT(I,J).GT.0) THEN
CSJA            PFLOPA(3+MIN(ICOUNT(I,J),5),I)=
CSJA     >           PFLOPA(3+MIN(ICOUNT(I,J),5),I)+1
CSJA            IF (TKSTHST.GT.0) THEN
CSJA*
CSJA              if(i.eq.DELID(ID_FCH).and.LOG10(PFPMC(J)).gt.1.0) then
CSJA                 call hf1(87,LOG10(PFPMC(J)),1.0)
CSJA                 call hf1(18,ABS(PFCTMC(J)),1.0)
CSJA              endif                             !IBozovic
CSJA*
CSJA              CALL HF1(88000+I,LOG10(PFPMC(J)),1.0)
CSJA              CALL HF1(88100+I,ABS(PFCTMC(J)),1.0)
CSJA            END IF
CSJA          END IF
CSJA        END DO
CSJA      END DO
CSJA
CSJA
CSJAC.....--------------------------------------------------------
CSJAC.....=== FIND LONGEST CORRECT TS AND TK FOR EACH MC TRACK ===
CSJAC.....--------------------------------------------------------
CSJA
CSJAC.....count number of good TS for each MC track that fulfils cuts
CSJAC     and count all fake TS that fulfil cuts
CSJA      CALL VZERO(NUMTS,TKNTRK+1)
CSJA      DO I=1,NTS
CSJA        IF (ITSDAT(5,I).EQ.0) THEN
CSJA          IF (.NOT.TKFCOUNT(RTS(20,I),RTS(22,I),640)) GOTO 10
CSJA        ELSE
CSJA          IF (.NOT.TKCOUNT(IABS(ITSDAT(5,I)))) GOTO 10
CSJA        END IF
CSJA        IF (NUMTS(IABS(ITSDAT(5,I))).EQ.0) THEN
CSJA          NUMTS(IABS(ITSDAT(5,I)))=1
CSJA          BESTTS(IABS(ITSDAT(5,I)))=I
CSJA          LENTS(IABS(ITSDAT(5,I)))=TKNHIT('TS',I)
CSJA        ELSE IF (TKNHIT('TS',I).GT.LENTS(IABS(ITSDAT(5,I)))) THEN
CSJA          NUMTS(IABS(ITSDAT(5,I)))=NUMTS(IABS(ITSDAT(5,I)))+1
CSJA          BESTTS(IABS(ITSDAT(5,I)))=I
CSJA          LENTS(IABS(ITSDAT(5,I)))=TKNHIT('TS',I)
CSJA        ELSE
CSJA          NUMTS(IABS(ITSDAT(5,I)))=NUMTS(IABS(ITSDAT(5,I)))+1
CSJA        END IF
CSJA 10     CONTINUE
CSJA      END DO
CSJA
CSJAC.....count number of good TK for each MC track that fulfils cuts
CSJAC     and count all fake TK that fulfil cuts
CSJA      CALL VZERO(NUMTK,TKNTRK+1)
CSJA      DO I=1,NTK
CSJA        IF (ITKDAT(3,I).EQ.0) THEN
CSJA          IF (.NOT.TKFCOUNT(RTK(20,I),RTK(22,I),640)) GOTO 20
CSJA        ELSE
CSJA          IF (.NOT.TKCOUNT(IABS(ITKDAT(3,I)))) GOTO 20
CSJA        END IF
CSJA        IF (NUMTK(IABS(ITKDAT(3,I))).EQ.0) THEN
CSJA          NUMTK(IABS(ITKDAT(3,I)))=1
CSJA          BESTTK(IABS(ITKDAT(3,I)))=I
CSJA          LENTK(IABS(ITKDAT(3,I)))=TKNHIT('TK',I)
CSJA        ELSE IF (TKNHIT('TK',I).GT.LENTK(IABS(ITKDAT(3,I)))) THEN
CSJA          NUMTK(IABS(ITKDAT(3,I)))=NUMTK(IABS(ITKDAT(3,I)))+1
CSJA          BESTTK(IABS(ITKDAT(3,I)))=I
CSJA          LENTK(IABS(ITKDAT(3,I)))=TKNHIT('TK',I)
CSJA        ELSE
CSJA          NUMTK(IABS(ITKDAT(3,I)))=NUMTK(IABS(ITKDAT(3,I)))+1
CSJA        END IF
CSJA 20     CONTINUE
CSJA      END DO
CSJA
CSJAC.....--------------------------------
CSJAC.....=== TRACK SEARCH PERFORMANCE ===
CSJAC.....--------------------------------
CSJA
CSJAC.....count number of tracks with at least one and at least three hits
CSJAC     and how many TS have been found for each reconstructed track
CSJA      DO I=1,TKNTRK
CSJA        IF (TKCOUNT(I)) THEN
CSJA          IF (ITKMCT(10,I).GE.1) THEN
CSJA            PFDELS(1)=PFDELS(1)+1      ! this is only approximately true
CSJA            PFDELS(2)=PFDELS(2)+1
CSJA            IF (ITKMCT(10,I).GE.3) THEN
CSJA              PFDELS(3)=PFDELS(3)+1
CSJA            END IF
CSJA          END IF
CSJA          IF (NUMTS(I).GT.0) THEN
CSJA            PFDELS(3+MIN(NUMTS(I),5))=PFDELS(3+MIN(NUMTS(I),5))+1
CSJA          END IF
CSJA        END IF
CSJA      END DO
CSJA
CSJAC.....number of fake TS
CSJA      PFDELS(9)=PFDELS(9)+NUMTS(0)
CSJA
CSJAC.....--------------------------
CSJAC.....=== GLOBAL PERFORMANCE ===
CSJAC.....--------------------------
CSJA
CSJAC.....count number of tracks with at least one and at least three hits
CSJAC     and how many TK have been found for each reconstructed track
CSJA      DO I=1,TKNTRK
CSJA        IF (TKCOUNT(I)) THEN
CSJA          IF (ITKMCT(10,I).GE.1) THEN
CSJA            PFGLOB(1)=PFGLOB(1)+1          ! this is only approximately true
CSJA            PFGLOB(2)=PFGLOB(2)+1
CSJA            IF (ITKMCT(10,I).GE.3) THEN
CSJA              PFGLOB(3)=PFGLOB(3)+1
CSJA            END IF
CSJA          END IF
CSJA          IF (NUMTK(I).GT.0) THEN
CSJA            PFGLOB(3+MIN(NUMTK(I),5))=PFGLOB(3+MIN(NUMTK(I),5))+1
CSJA          END IF
CSJA        END IF
CSJA      END DO
CSJA
CSJAC.....number of fake TK
CSJA      PFGLOB(9)=PFGLOB(9)+NUMTK(0)
CSJA
CSJAC.....----------------------------------------------
CSJAC.....=== SUBDETECTOR CONTRIBUTION TO LONGEST TS ===
CSJAC.....----------------------------------------------
CSJA
CSJAC.....reminder: number of tracks passing each detector is in PFLOPA
CSJA
CSJAC.....loop over all TEs of track and check if they are matching, wrong or noise
CSJA      DO I=1,TKNTRK
CSJA        IF (TKCOUNT(I).AND.(NUMTS(I).GT.0)) THEN
CSJA          K=0
CSJA          DO J=1,ITSDAT(2,BESTTS(I))
CSJA            TENUM=ITSTEL(ITSDAT(1,BESTTS(I))+J-1)
CSJA            TELAB=IABS(ITEDAT(5,TENUM))
CSJA            IF (TELAB.EQ.I) THEN
CSJA              PFTETS(1,DELID(ITE(1,TENUM)))=
CSJA     >             PFTETS(1,DELID(ITE(1,TENUM)))+1
CSJA              ID=ITE(1,TENUM)/100*100
CSJA              IF (ID.EQ.ID_VTX.OR.ID.EQ.ID_SIT.OR.ID.EQ.ID_FTD) THEN
CSJA                IF (K.EQ.0) THEN
CSJA                  K=1
CSJA                  PFTETS(1,DELID(ID_VTX))=PFTETS(1,DELID(ID_VTX))+1
CSJA                END IF
CSJA              END IF
CSJA            ELSE IF (TELAB.EQ.0) THEN
CSJA              PFTETS(3,DELID(ITE(1,TENUM)))=
CSJA     >             PFTETS(3,DELID(ITE(1,TENUM)))+1
CSJA            ELSE
CSJA              PFTETS(2,DELID(ITE(1,TENUM)))=
CSJA     >             PFTETS(2,DELID(ITE(1,TENUM)))+1
CSJA            END IF
CSJA          END DO
CSJA        END IF
CSJA      END DO
CSJA
CSJAC.....loop over all fake TS and check for good/noise TEs in there
CSJA      DO I=1,NTS
CSJA        IF (ITSDAT(5,I).EQ.0) THEN
CSJA          IF (TKFCOUNT(RTS(20,I),RTS(22,I),640)) THEN
CSJA            DO J=1,ITSDAT(2,I)
CSJA              TENUM=ITSTEL(ITSDAT(1,I)+J-1)
CSJA              TELAB=IABS(ITEDAT(5,TENUM))
CSJA              IF (TELAB.GT.0) THEN
CSJA                PFTETS(4,DELID(ITE(1,TENUM)))=
CSJA     >               PFTETS(4,DELID(ITE(1,TENUM)))+1
CSJA              ELSE
CSJA                PFTETS(5,DELID(ITE(1,TENUM)))=
CSJA     >               PFTETS(5,DELID(ITE(1,TENUM)))+1
CSJA              END IF
CSJA            END DO
CSJA          END IF
CSJA        END IF
CSJA      END DO
CSJA
CSJA
CSJAC.....----------------------------------------------
CSJAC.....=== SUBDETECTOR CONTRIBUTION TO LONGEST TK ===
CSJAC.....----------------------------------------------
CSJA
CSJAC.....reminder: number of tracks passing each detector is in PFLOPA
CSJA
CSJAC.....loop over all TEs of track and check if they are matching, wrong or noise
CSJA      DO I=1,TKNTRK
CSJA        IF (TKCOUNT(I).AND.(NUMTK(I).GT.0)) THEN
CSJA          K=0
CSJA          DO J=1,ITKDAT(2,BESTTK(I))
CSJA            TENUM=ITKTEL(ITKDAT(1,BESTTK(I))+J-1)
CSJA            TELAB=IABS(ITEDAT(5,TENUM))
CSJA            IF (TELAB.EQ.I) THEN
CSJA              PFTETK(1,DELID(ITE(1,TENUM)))=
CSJA     >             PFTETK(1,DELID(ITE(1,TENUM)))+1
CSJA              ID=ITE(1,TENUM)/100*100
CSJA              IF (ID.EQ.ID_VTX.OR.ID.EQ.ID_SIT.OR.ID.EQ.ID_FTD) THEN
CSJA                IF (K.EQ.0) THEN
CSJA                  K=1
CSJA                  PFTETK(1,DELID(ID_VTX))=PFTETK(1,DELID(ID_VTX))+1
CSJA                END IF
CSJA              END IF
CSJA            ELSE IF (TELAB.EQ.0) THEN
CSJA              PFTETK(3,DELID(ITE(1,TENUM)))=
CSJA     >             PFTETK(3,DELID(ITE(1,TENUM)))+1
CSJA            ELSE
CSJA              PFTETK(2,DELID(ITE(1,TENUM)))=
CSJA     >             PFTETK(2,DELID(ITE(1,TENUM)))+1
CSJA            END IF
CSJA          END DO
CSJA        END IF
CSJA      END DO
CSJA
CSJAC.....loop over all fake TK and check for good/noise TEs in there
CSJA      DO I=1,NTK
CSJA        IF (ITKDAT(3,I).EQ.0) THEN
CSJA          IF (TKFCOUNT(RTK(20,I),RTK(22,I),640)) THEN
CSJA            DO J=1,ITKDAT(2,I)
CSJA              TENUM=ITKTEL(ITKDAT(1,I)+J-1)
CSJA              TELAB=IABS(ITEDAT(5,TENUM))
CSJA              IF (TELAB.GT.0) THEN
CSJA                PFTETK(4,DELID(ITE(1,TENUM)))=
CSJA     >               PFTETK(4,DELID(ITE(1,TENUM)))+1
CSJA              ELSE
CSJA                PFTETK(5,DELID(ITE(1,TENUM)))=
CSJA     >               PFTETK(5,DELID(ITE(1,TENUM)))+1
CSJA              END IF
CSJA            END DO
CSJA          END IF
CSJA        END IF
CSJA      END DO
CSJA
CSJA
CSJAC.....------------------
CSJAC.....=== HISTOGRAMS ===
CSJAC.....------------------
CSJA
CSJA      IF (TKSTHST.GT.0) THEN
CSJA
CSJA        CALL HCDIR('//PAWC/TRACKING/TKSTEER',' ')
CSJAC.......all generated tracks and all correctly reconstructed tracks
CSJA        DO I=1,TKNTRK
CSJA          IF (TKCOUNT(I)) THEN
CSJA            ! all generated tracks
CSJA            CALL HF1(99000,LOG10(PFPMC(I)),1.0)
CSJA            CALL HF1(99001,PFCTMC(I),1.0)
CSJA            ! generated tracks with at least three hits in detector
CSJA            IF (ITKMCT(10,I).GE.3) THEN
CSJA              CALL HF1(99100,LOG10(PFPMC(I)),1.0)
CSJA              CALL HF1(99101,PFCTMC(I),1.0)
CSJA            END IF
CSJA            ! reconstructed tracks
CSJA            IF (NUMTK(I).GT.0) THEN
CSJA              CALL HF1(99200,LOG10(PFPMC(I)),1.)
CSJA              CALL HF1(99201,PFCTMC(I),1.)
CSJA              ! resolution histograms
CSJA              IF (RTK(22,BESTTK(I)).NE.0.0) THEN
CSJA                FAKP=1./RTK(22,BESTTK(I))
CSJA              ELSE
CSJA                FAKP=1.E6
CSJA              END IF
CSJA              TRUEP=PFPMC(I)
CSJA              IF (TRUEP.NE.0) THEN
CSJA                CALL HF1(99350,(TRUEP-FAKP)/TRUEP,1.0)
CSJA                CALL HF1(99351,(1./TRUEP-1./FAKP)*TRUEP,1.0)
CSJA                CALL HF2(99352,(TRUEP-FAKP)/TRUEP,PFCTMC(I),1.0)
CSJA              END IF
CSJA            END IF
CSJA          END IF
CSJA        END DO
CSJA
CSJAC.......fake tracks
CSJA        DO I=1,NTK
CSJA          IF (ITKDAT(3,I).EQ.0) THEN
CSJA            FAKCT=COS(RTK(20,I))
CSJA            IF (RTK(22,I).NE.0.0) THEN
CSJA              FAKP=ABS(1./RTK(22,I))
CSJA            ELSE
CSJA              FAKP=1.E6
CSJA            END IF
CSJA            IF (TKFCOUNT(RTK(20,I),RTK(22,I),640)) THEN
CSJA              CALL HF1(99300,LOG10(FAKP),1.0)
CSJA              CALL HF1(99301,FAKCT,1.0)
CSJA            END IF
CSJA          END IF
CSJA        END DO
CSJA
CSJAC.......TE contribution to reconstructed tracks (longest piece)
CSJA        DO I=1,TKNTRK
CSJA          IF (TKCOUNT(I)) THEN
CSJA            ! check for each detector if it was hit by this track
CSJA            DO J=1,DETMAX
CSJA              IF (IHCOUNT(J,I).GT.0) THEN
CSJA                CALL HF1(99400+J,LOG10(PFPMC(I)),1.0)
CSJA                CALL HF1(99500+J,ABS(PFCTMC(I)),1.0)
CSJA                IF (IHCOUNT(J,I).GE.3) THEN
CSJA                  if(j.eq.delid(ID_FCH).and.LOG10(PFPMC(I)).gt.1.) then
CSJA                    call HF1(107,LOG10(PFPMC(I)),1.0)
CSJA                    call HF1(17,ABS(PFCTMC(I)),1.0)
CSJA                  endif                            !IBozovic
CSJA                  CALL HF1(99600+J,LOG10(PFPMC(I)),1.0)
CSJA                  CALL HF1(99700+J,ABS(PFCTMC(I)),1.0)
CSJA                END IF
CSJA              END IF
CSJA            END DO
CSJA            ! now check for each detector whether it contributed a TE
CSJA            IF (NUMTK(I).GT.0) THEN
CSJA              K=0
CSJA              DO J=1,ITKDAT(2,BESTTK(I))
CSJA                TENUM=ITKTEL(ITKDAT(1,BESTTK(I))+J-1)
CSJA                IF (IABS(ITEDAT(5,TENUM)).EQ.I) THEN
CSJA                  ID=DELID(ITE(1,TENUM))
CSJA                  if(id.eq.delid(ID_FCH).and.LOG10(PFPMC(I)).gt.1.) then
CSJA                     call HF1(97,LOG10(PFPMC(I)),1.0)
CSJA                     call HF1(19,ABS(PFCTMC(I)),1.0)
CSJA                  endif                               !IBozovic
CSJA                  CALL HF1(99800+ID,LOG10(PFPMC(I)),1.0)
CSJA                  CALL HF1(99900+ID,ABS(PFCTMC(I)),1.0)
CSJA                  ID=ITE(1,TENUM)/100*100
CSJA                  IF (ID.EQ.ID_VTX.OR.ID.EQ.ID_SIT.OR.ID.EQ.ID_FTD) THEN
CSJA                    ! count this TE also for generic VTX tracks
CSJA                    ! but take care you do this only once per track
CSJA                    IF (K.EQ.0) THEN
CSJA                      K=1
CSJA                      CALL HF1(99800+DELID(ID_VTX),LOG10(PFPMC(I)),1.0)
CSJA                      CALL HF1(99900+DELID(ID_VTX),ABS(PFCTMC(I)),1.0)
CSJA                    END IF
CSJA                  END IF
CSJA                END IF
CSJA              END DO
CSJA            END IF
CSJA          END IF
CSJA        END DO
CSJA
CSJA        CALL HCDIR('//PAWC',' ')
CSJA      END IF
CSJA
CSJA
CSJAC.....-------------------------------------------------------------------------
CSJAC.....=== LIST ALL MC tracks,TE,TS,TK and classify anything that went wrong ===
CSJAC.....-------------------------------------------------------------------------
CSJA
CSJAC.....skip list if not in debug mode
CSJA      IF (TKSTDBG.EQ.0) GOTO 999
CSJA
CSJAC.....sort MC tracks with respect to GEANT track ID
CSJA      DO I=1,MIN(TKNTRK,MXSORT)
CSJA        ISORT(I)=ITKMCT(9,I)
CSJA      END DO
CSJA      CALL SORTZV(ISORT,INDEX,MIN(TKNTRK,MXSORT),-1,0,0)
CSJA
CSJAC.....list MC tracks
CSJA      WRITE(6,*)
CSJA      WRITE(6,'(A50)')  '========================='//
CSJA     >                  '========================='
CSJA      WRITE(6,*) 'total number of tree tracks: ',TKNTRK
CSJA      WRITE(6,'(A50)')  '-------------------------'//
CSJA     >                  '-------------------------'
CSJA      WRITE(6,'(A3,2X,A8,2X,A8,2X,A4,2X,A12)')
CSJA     >     'no.','GEANT ID','momentum','hits','quality cuts'
CSJA      WRITE(6,'(A50)')  '-------------------------'//
CSJA     >                  '-------------------------'
CSJA      DO K=1,TKNTRK
CSJA        I=K
CSJA        IF (K.LE.MXSORT) I=INDEX(K)
CSJA        IF (TKCOUNT(I)) THEN
CSJA          USE='YES'
CSJA        ELSE
CSJA          USE='NO'
CSJA        END IF
CSJA        WRITE(6,'(I3,2X,I8,2X,F8.3,2X,I4,2X,A12)')
CSJA     >       I,ITKMCT(9,I),
CSJA     >       SQRT(TKMCTR(1,I)**2+TKMCTR(2,I)**2+TKMCTR(3,I)**2),
CSJA     >       ITKMCT(10,I),USE
CSJA      END DO
CSJA      WRITE(6,'(A50)')  '-------------------------'//
CSJA     >                  '-------------------------'
CSJA


CSJAC.....sort TEs with respect to GEANT track ID and detector ID
CSJA      DO I=1,MIN(NTE,MXSORT)
CSJA        IF (ITEDAT(5,I).EQ.0) THEN
CSJA          ISORT(I)=0
CSJA        ELSE
CSJA          ISORT(I)=ITKMCT(9,IABS(ITEDAT(5,I)))
CSJA        END IF
CSJA        ISORT(I)=ISORT(I)*100+DELID(ITE(1,I))
CSJA      END DO
CSJA      CALL SORTZV(ISORT,INDEX,MIN(NTE,MXSORT),-1,0,0)
CSJA
C.....print sorted list of TEs
      WRITE(6,*)
      WRITE(6,'(A40)')  '===================='//
     >                  '===================='
      WRITE(6,*) 'total number of TE: ',NTE
      WRITE(6,'(A40)')  '--------------------'//
     >                  '--------------------'
      WRITE(6,'(A3,2X,A5,2X,A8,2X,A4,2X,A12)')
     >     'TE','det. ','GEANT ID','hits','quality cuts'
      WRITE(6,'(A40)')  '--------------------'//
     >                  '--------------------'
      DO K=1,NTE
         I=K
CSJA        IF (K.LE.MXSORT) I=INDEX(K)

CSJA        IF (ITEDAT(5,I).NE.0) THEN
CSJA          IF (TKCOUNT(IABS(ITEDAT(5,I)))) THEN
CSJA            USE='YES'
CSJA          ELSE
CSJA            USE='NO'
CSJA          END IF
CSJA          ID=ITKMCT(9,IABS(ITEDAT(5,I)))
CSJA        ELSE
CSJA          IF (TKFCOUNT(RTE(13,I),RTE(15,I),ITE(4,I))) THEN
CSJA            USE='YES'
CSJA          ELSE
CSJA            USE='NO'
CSJA          END IF
CSJA          ID=0
CSJA        END IF

         IF (DELID(ITE(1,I)).NE.0) THEN
            WRITE(6,'(I3,2X,A5,2X,I8,2X,I4,2X,A12)')
     >           I,DETNAME(DELID(ITE(1,I))),ID,ITEDAT(2,I),USE
         ELSE
            WRITE(6,'(I3,2X,I5,2X,I8,2X,I4,2X,A12)')
     >           I,ITE(1,I),ID,ITEDAT(2,I),USE
         END IF

      END DO
      WRITE(6,'(A40)')  '===================='//
     >                  '===================='




CSJAC.....sort TS with respect to GEANT track ID
CSJA      DO I=1,MIN(NTS,MXSORT)
CSJA        IF (ITSDAT(5,I).EQ.0) THEN
CSJA          ISORT(I)=0
CSJA        ELSE
CSJA          ISORT(I)=ITKMCT(9,IABS(ITSDAT(5,I)))
CSJA        END IF
CSJA      END DO
CSJA      CALL SORTZV(ISORT,INDEX,MIN(NTS,MXSORT),-1,0,0)
CSJA

C.....print sorted list of TS
      WRITE(6,*)
      WRITE(6,'(A56)')  '============================'//
     >                  '============================'
      WRITE(6,*) 'total number of TS: ',NTS
      WRITE(6,'(A56)')  '----------------------------'//
     >                  '----------------------------'
      WRITE(6,'(A3,2X,A5,2X,A8,2(2X,A9),2X,A12)')
     >    'TS','det. ','GEANT ID','true hits','used hits','quality cuts'

      DO K=1,NTS
        I=K
CSJA        IF (I.LE.MXSORT) I=INDEX(K)
        WRITE(6,'(A56)')  '----------------------------'//
     >                    '----------------------------'
CSJA        IF (ITSDAT(5,I).NE.0) THEN
CSJA          IF (TKCOUNT(IABS(ITSDAT(5,I)))) THEN
CSJA            USE='YES'
CSJA          ELSE
CSJA            USE='NO'
CSJA          END IF
CSJA          ID=ITKMCT(9,IABS(ITSDAT(5,I)))
CSJA          NHTREE=ITKMCT(10,IABS(ITSDAT(5,I)))
CSJA        ELSE
CSJA          IF (TKFCOUNT(RTS(20,I),RTS(22,I),640)) THEN
CSJA            USE='YES'
CSJA          ELSE
CSJA            USE='NO'
CSJA          END IF
CSJA          ID=0
CSJA          NHTREE=0
CSJA        END IF

        WRITE(6,'(I3,2X,A5,2X,I8,2(2X,I9),2X,A12)')
     >       I,' ',
     >       ID,NHTREE,TKNHIT('TS',I),USE

        DO J=1,DETMAX
          IF (J.EQ.I_VTX.OR.J.EQ.I_FTD) GOTO 350
          IF (ITSDAT(5,I).EQ.0) THEN
             NHTREE=0
          ELSE
             NHTREE=IHCOUNT(J,IABS(ITSDAT(5,I)))
          END IF

CSJA          ID=IABS(ITSDAT(5,I))
CSJA          IF (ID.NE.0) ID=ITKMCT(9,ID)

          NHREC=0
          DO L=0,ITSDAT(2,I)-1
             IF (DELID(ITE(1,ITSTEL(ITSDAT(1,I)+L))).EQ.J) THEN
CSJA              ID=ITEDAT(5,ITSTEL(ITSDAT(1,I)+L))
CSJA              IF (ID.NE.0) ID=ITKMCT(9,IABS(ID))
             NHREC=ITEDAT(2,ITSTEL(ITSDAT(1,I)+L))
          END IF
          END DO
          IF (NHTREE.GT.0.OR.NHREC.GT.0) THEN
            WRITE(6,'(5X,A5,2X,I8,2(2X,I9))') DETNAME(J),ID,NHTREE,NHREC
          END IF
 350      CONTINUE
        END DO
      END DO
      WRITE(6,'(A56)')  '============================'//
     >                  '============================'




CSJAC.....sort TK with respect to GEANT track ID
CSJA      DO I=1,MIN(NTK,MXSORT)
CSJA        IF (ITKDAT(3,I).EQ.0) THEN
CSJA          ISORT(I)=0
CSJA        ELSE
CSJA          ISORT(I)=ITKMCT(9,IABS(ITKDAT(3,I)))
CSJA        END IF
CSJA      END DO
CSJA      CALL SORTZV(ISORT,INDEX,MIN(NTK,MXSORT),-1,0,0)
CSJA
CSJAC.....print sorted list of TK
      WRITE(6,*)
      WRITE(6,'(A56)')  '============================'//
     >                  '============================'
      WRITE(6,*) 'total number of TK: ',NTK
      WRITE(6,'(A56)')  '----------------------------'//
     >                  '----------------------------'
      WRITE(6,'(A3,2X,A5,2X,A8,2(2X,A9),2X,A12)')
     >    'TK','det. ','GEANT ID','true hits','used hits','quality cuts'

      DO K=1,NTK
        I=K
CSJA        IF (I.LE.MXSORT) I=INDEX(K)
        WRITE(6,'(A56)')  '----------------------------'//
     >                    '----------------------------'
CSJA        IF (ITKDAT(3,I).NE.0) THEN
CSJA          IF (TKCOUNT(IABS(ITKDAT(3,I)))) THEN
CSJA            USE='YES'
CSJA          ELSE
CSJA            USE='NO'
CSJA          END IF
CSJA          ID=ITKMCT(9,IABS(ITKDAT(3,I)))
CSJA          NHTREE=ITKMCT(10,IABS(ITKDAT(3,I)))
CSJA        ELSE
CSJA          IF (TKFCOUNT(RTK(20,I),RTK(22,I),640)) THEN
CSJA            USE='YES'
CSJA          ELSE
CSJA            USE='NO'
CSJA          END IF
CSJA          ID=0
CSJA          NHTREE=0
CSJA        END IF

        WRITE(6,'(I3,2X,A5,2X,I8,2(2X,I9),2X,A12)')
     >       I,' ',
     >       ID,NHTREE,TKNHIT('TK',I),USE

        DO J=1,DETMAX
          IF (J.EQ.I_VTX.OR.J.EQ.I_FTD) GOTO 360
          IF (ITKDAT(3,I).EQ.0) THEN
            ID=0
            NHTREE=0
          ELSE
CSJA            ID=ITKMCT(9,IABS(ITKDAT(3,I)))
            NHTREE=IHCOUNT(J,IABS(ITKDAT(3,I)))
          END IF
          NHREC=0
          DO L=0,ITKDAT(2,I)-1
            IF (DELID(ITE(1,ITKTEL(ITKDAT(1,I)+L))).EQ.J) THEN
c              ID=ITEDAT(5,ITKTEL(ITKDAT(1,I)+L))
c              IF (ID.NE.0) ID=ITKMCT(9,IABS(ID))
              NHREC=ITEDAT(2,ITKTEL(ITKDAT(1,I)+L))
            END IF
          END DO
          IF (NHTREE.GT.0.OR.NHREC.GT.0) THEN
            WRITE(6,'(5X,A5,2X,I8,2(2X,I9))') DETNAME(J),ID,NHTREE,NHREC
          END IF
 360      CONTINUE
        END DO
CSJA
CSJAC.......now classify fake TK
CSJA        IF (ITKDAT(3,I).EQ.0.AND.USE.EQ.'YES') THEN
CSJA          ! get TPC TE label and classify SI hits
CSJA          ID=-1
CSJA          NLIST=0
CSJA          SIHITS=0
CSJA          DO L=0,ITKDAT(2,I)-1
CSJA            IF (ITE(1,ITKTEL(ITKDAT(1,I)+L)).EQ.ID_TPC) THEN
CSJA              ID=IABS(ITEDAT(5,ITKTEL(ITKDAT(1,I)+L)))
CSJA            END IF
CSJA            IDGEN=ITE(1,ITKTEL(ITKDAT(1,I)+L))/100*100
CSJA            IF (IDGEN.EQ.ID_VTX.OR.IDGEN.EQ.ID_SIT
CSJA     >                                 .OR.IDGEN.EQ.ID_FTD) THEN
CSJA              SIHITS=SIHITS+1
CSJA              DO M=1,NLIST
CSJA                IF (ILAB(M).EQ.ITEDAT(5,ITKTEL(ITKDAT(1,I)+L))) THEN
CSJA                  ILIST(M)=ILIST(M)+1
CSJA                  GOTO 1234
CSJA                END IF
CSJA              END DO
CSJA              IF (NLIST.LT.MLIST) NLIST=NLIST+1
CSJA              ILAB(NLIST)=ITEDAT(5,ITKTEL(ITKDAT(1,I)+L))
CSJA              ILIST(NLIST)=1
CSJA 1234         CONTINUE
CSJA            END IF
CSJA          END DO
CSJA          LSI=.FALSE.
CSJA          DO M=1,NLIST
CSJA            IF (FLOAT(ILIST(M))/FLOAT(SIHITS).GT.0.7) THEN
CSJA               LSI=.TRUE.
CSJA            END IF
CSJA          END DO
CSJA
CSJA          IF (ITKDAT(2,I).EQ.1.AND.ID.EQ.0) THEN
CSJA            ! TPC noise TE and nothing else
CSJA            PFFAKE(1)=PFFAKE(1)+1
CSJA          ELSE IF (ID.EQ.0.AND.LSI) THEN
CSJA            ! good SI track ending in TPC noise
CSJA            PFFAKE(2)=PFFAKE(2)+1
CSJA          ELSE IF (ID.NE.0.AND..NOT.LSI) THEN
CSJA            ! good TPC track (or no TPC at all) with SI noise
CSJA            PFFAKE(3)=PFFAKE(3)+1
CSJA          ELSE IF (ID.GT.0.AND.LSI) THEN
CSJA            ! wrong merging
CSJA            PFFAKE(4)=PFFAKE(4)+1
CSJA          ELSE
CSJA            PFFAKE(5)=PFFAKE(5)+1
CSJA          END IF
CSJA        END IF
CSJA      END DO
CSJA      WRITE(6,'(A56)')  '============================'//
CSJA     >                  '============================'
CSJA
CSJAC.....find first track that has not been reconstructed (if any)
CSJA      I=1
CSJA      ID=0
CSJA 321  IF (TKCOUNT(I).AND.NUMTK(I).EQ.0) THEN
CSJA        ID=I
CSJA      ELSE
CSJA        I=I+1
CSJA        IF (I.LE.TKNTRK) GOTO 321
CSJA      END IF
CSJA
CSJAC.....print list of missed tracks (purpose: check out WHY they were missed)
CSJAC     this is only done if there is at least one missed track in this event
CSJA      IF (ID.GT.0) THEN
CSJA        WRITE(6,*)
CSJA        WRITE(6,'(A40)')  '===================='//
CSJA     >                    '===================='
CSJA        WRITE(6,*) 'Tracks that have not been reconstructed:'
CSJA        WRITE(6,'(A40)')  '--------------------'//
CSJA     >                    '--------------------'
CSJA        DO I=ID,TKNTRK
CSJA          IF (TKCOUNT(I).AND.NUMTK(I).EQ.0) THEN
CSJA            WRITE(6,*) 'MC tree track number ',I
CSJA            WRITE(6,*) '  momentum  : ',PFPMC(I),' GeV/c'
CSJA            WRITE(6,*) '  cos(theta): ',PFCTMC(I)
CSJA            WRITE(6,*) '  theta     : ',ACOS(PFCTMC(I))*180./PI,' deg'
CSJA            WRITE(6,*) '  no.of hits: ',ITKMCT(10,I)
CSJA            DO J=1,DETMAX
CSJA              IF (IHCOUNT(J,I).GT.0) THEN
CSJA                WRITE(6,*) '   ',DETNAME(J),': ',IHCOUNT(J,I)
CSJA              END IF
CSJA              ICOUNT(J,0)=0
CSJA            END DO
CSJA            TENUM=0
CSJA            SIHITS=0
CSJA            DO J=1,NTE
CSJA              IF (IABS(ITEDAT(5,J)).EQ.I) THEN
CSJA                IF (DELID(ITE(1,J)).GT.0) THEN
CSJA                  ICOUNT(DELID(ITE(1,J)),0)=ICOUNT(DELID(ITE(1,J)),0)+1
CSJA                  ID=ITE(1,J)/100*100
CSJA                  IF (ID.EQ.ID_VTX.OR.ID.EQ.ID_SIT.OR.ID.EQ.ID_FTD) THEN
CSJA                    ! do not count VTX track TE as hit here
CSJA                    IF (ID.NE.ITE(1,J)) THEN
CSJA                       SIHITS=SIHITS+1
CSJA                    END IF
CSJA                  END IF
CSJA                END IF
CSJA                TENUM=TENUM+1
CSJA              END IF
CSJA            END DO
CSJA            WRITE(6,*) '  no. of TEs: ',TENUM
CSJA            DO J=1,DETMAX
CSJA              IF (ICOUNT(J,0).GT.0) THEN
CSJA                WRITE(6,*) '   ',DETNAME(J),': ',ICOUNT(J,0)
CSJA              END IF
CSJA            END DO
CSJA            WRITE(6,*) '  vertex x  : ',TKMCTR(5,I),' cm'
CSJA            WRITE(6,*) '         y  : ',TKMCTR(6,I),' cm'
CSJA            WRITE(6,*) '  vertex r  : ',
CSJA     >                 SQRT(TKMCTR(5,I)**2+TKMCTR(6,I)**2),' cm'
CSJA            WRITE(6,*) '         z  : ',TKMCTR(7,I),' cm'
CSJA            WRITE(6,*) '  part. code: ',ITKMCT(8,I)
CSJA            WRITE(6,*) '  GEANT trk : ',ITKMCT(9,I)
CSJA            WRITE(6,*) '  HEPEVT trk: ',ITKMCT(11,I)
CSJA            WRITE(6,'(A40)')  '--------------------'//
CSJA     >                        '--------------------'
CSJA
CSJAc           now classify missing tracks
CSJA            IF (IHCOUNT(DELID(ID_TPC),I).LT.3.AND.SIHITS.LT.3) THEN
CSJA               ! too few hits in both primary detectors
CSJA               PFFAIL(1)=PFFAIL(1)+1
CSJA            ELSE IF (IHCOUNT(DELID(ID_TPC),I).LT.3.AND.
CSJA     >               ICOUNT(DELID(ID_VTX),0).EQ.0) THEN
CSJA               ! too few TPC hits, and SI patrec failed
CSJA               PFFAIL(2)=PFFAIL(2)+1
CSJA            ELSE IF (SIHITS.LT.3.AND.ICOUNT(DELID(ID_TPC),0).EQ.0) THEN
CSJA               ! too few SI hits, and TPC patrec failed
CSJA               PFFAIL(3)=PFFAIL(3)+1
CSJA            ELSE IF (ICOUNT(DELID(ID_TPC),0).EQ.0.AND.
CSJA     >               ICOUNT(DELID(ID_VTX),0).EQ.0) THEN
CSJA               ! both patrecs failed
CSJA               PFFAIL(4)=PFFAIL(4)+1
CSJA            ELSE
CSJA               ! others, most likely global track search problems (fit?)
CSJA               PFFAIL(5)=PFFAIL(5)+1
CSJA            END IF
CSJA
CSJA          END IF
CSJA        END DO
CSJA        WRITE(6,'(A40)')  '===================='//
CSJA     >                    '===================='
CSJA      END IF
CSJA
CSJAC.....find first track that has been split into several TK (if any)
CSJA      I=1
CSJA      ID=0
CSJA 322  IF (TKCOUNT(I).AND.NUMTK(I).GT.1) THEN
CSJA        ID=I
CSJA      ELSE
CSJA        I=I+1
CSJA        IF (I.LE.TKNTRK) GOTO 322
CSJA      END IF
CSJA
CSJAC.....print list of split tracks (i.e. all tracks with more than one TK)
CSJAC     this is only done if there is at least one split track in this event
CSJA      IF (ID.GT.0) THEN
CSJA        WRITE(6,*)
CSJA        WRITE(6,'(A40)')  '===================='//
CSJA     >                    '===================='
CSJA        WRITE(6,*) 'Tracks that have been split up into several:'
CSJA        WRITE(6,'(A40)')  '--------------------'//
CSJA     >                    '--------------------'
CSJA        DO I=ID,TKNTRK
CSJA          IF (TKCOUNT(I).AND.NUMTK(I).GT.1) THEN
CSJA            WRITE(6,*) 'MC tree track number ',I
CSJA            WRITE(6,*) '  momentum  : ',PFPMC(I),' GeV/c'
CSJA            WRITE(6,*) '  cos(theta): ',PFCTMC(I)
CSJA            WRITE(6,*) '  theta     : ',ACOS(PFCTMC(I))*180./PI,' deg'
CSJA            WRITE(6,*) '  no.of hits: ',ITKMCT(10,I)
CSJA            DO J=1,DETMAX
CSJA              IF (IHCOUNT(J,I).GT.0) THEN
CSJA                WRITE(6,*) '   ',DETNAME(J),': ',IHCOUNT(J,I)
CSJA              END IF
CSJA              ICOUNT(J,0)=0
CSJA            END DO
CSJA            TENUM=0
CSJA            DO J=1,NTE
CSJA              IF (IABS(ITEDAT(5,J)).EQ.I) THEN
CSJA                IF (DELID(ITE(1,J)).GT.0) THEN
CSJA                  ICOUNT(DELID(ITE(1,J)),0)=ICOUNT(DELID(ITE(1,J)),0)+1
CSJA                END IF
CSJA                TENUM=TENUM+1
CSJA              END IF
CSJA            END DO
CSJA            WRITE(6,*) '  no. of TEs: ',TENUM
CSJA            DO J=1,DETMAX
CSJA              IF (ICOUNT(J,0).GT.0) THEN
CSJA                WRITE(6,*) '   ',DETNAME(J),': ',ICOUNT(J,0)
CSJA              END IF
CSJA            END DO
CSJA            WRITE(6,*) '  vertex x  : ',TKMCTR(5,I),' cm'
CSJA            WRITE(6,*) '         y  : ',TKMCTR(6,I),' cm'
CSJA            WRITE(6,*) '  vertex r  : ',
CSJA     >                 SQRT(TKMCTR(5,I)**2+TKMCTR(6,I)**2),' cm'
CSJA            WRITE(6,*) '         z  : ',TKMCTR(7,I),' cm'
CSJA            WRITE(6,*) '  part. code: ',ITKMCT(8,I)
CSJA            WRITE(6,*) '  GEANT trk : ',ITKMCT(9,I)
CSJA            WRITE(6,*) '  HEPEVT trk: ',ITKMCT(11,I)
CSJA            WRITE(6,*) '  => number of TKs: ',NUMTK(I)
CSJA
CSJAC...........list TK composition and collect information on participating TEs
CSJA            NLIST=0
CSJA            DO J=1,NTK
CSJA              IF (IABS(ITKDAT(3,J)).EQ.I) THEN
CSJA                WRITE(6,*) '    TK number ',J,': TEs from'
CSJA                DO K=ITKDAT(1,J),ITKDAT(1,J)+ITKDAT(2,J)-1
CSJA                  ! does this TE have the right label?
CSJA                  IF (IABS(ITEDAT(5,ITKTEL(K))).EQ.I) THEN
CSJA                    IF (NLIST.LT.MLIST) NLIST=NLIST+1
CSJA                    ILIST(NLIST)=ITKTEL(K)
CSJA                    ILAB(NLIST)=J
CSJA                    IF (IAND(ITE(4,ITKTEL(K)),1).EQ.1) THEN
CSJA                       RLIST(NLIST)=RTE(10,ITKTEL(K))
CSJA                    ELSE
CSJA                       RLIST(NLIST)=SQRT(RTE(10,ITKTEL(K))**2+
CSJA     >                                   RTE(11,ITKTEL(K))**2)
CSJA                    END IF
CSJA                  END IF
CSJA                  ! print information
CSJA                  IF (ITE(1,ITKTEL(K)).EQ.ID_TPC) THEN
CSJA                    WRITE(6,*) '       ',
CSJA     >                    DETNAME(DELID(ITE(1,ITKTEL(K)))),' (label ',
CSJA     >                    ITEDAT(5,ITKTEL(K)),'), inner radius ',
CSJA     >                    RTE(10,ITKTEL(K))
CSJA                  ELSE
CSJA                    WRITE(6,*) '       ',
CSJA     >                    DETNAME(DELID(ITE(1,ITKTEL(K)))),
CSJA     >                    ' (label ',ITEDAT(5,ITKTEL(K)),')'
CSJA                  END IF
CSJA                END DO
CSJA              END IF
CSJA            END DO
CSJA            WRITE(6,'(A40)')  '--------------------'//
CSJA     >                        '--------------------'
CSJA
CSJAC...........now use TE info collected above to classify track splitting
CSJA            ! sort TE by radius
CSJA            CALL SORTZV(RLIST,ISORT,NLIST,+1,0,0)
CSJA            ! scan through all TEs and check whether they are in different TK
CSJA            DO J=2,NLIST
CSJA              IF (ILAB(ISORT(J-1)).NE.ILAB(ISORT(J))) THEN
CSJA                ! we have a splitting here. check generic detector IDs
CSJA                K=MAX(ITE(1,ILIST(ISORT(J-1)))/100*100,
CSJA     >                ITE(1,ILIST(ISORT(J  )))/100*100)
CSJA                L=MIN(ITE(1,ILIST(ISORT(J-1)))/100*100,
CSJA     >                ITE(1,ILIST(ISORT(J  )))/100*100)
CSJA                IF (K.EQ.ID_VTX.AND.L.EQ.K) THEN
CSJA                   ! within VTX
CSJA                   PFSPLIT(1)=PFSPLIT(1)+1
CSJA                ELSE IF (K.EQ.MAX(ID_VTX,ID_SIT).AND.
CSJA     >                   L.EQ.MIN(ID_VTX,ID_SIT)) THEN
CSJA                   ! between VTX and SIT
CSJA                   PFSPLIT(2)=PFSPLIT(2)+1
CSJA                ELSE IF (K.EQ.MAX(ID_VTX,ID_FTD).AND.
CSJA     >                   L.EQ.MIN(ID_VTX,ID_FTD)) THEN
CSJA                   ! between VTX and FTD
CSJA                   PFSPLIT(3)=PFSPLIT(3)+1
CSJA                ELSE IF (K.EQ.ID_SIT.AND.L.EQ.K) THEN
CSJA                   ! within SIT
CSJA                   PFSPLIT(4)=PFSPLIT(4)+1
CSJA                ELSE IF (K.EQ.ID_FTD.AND.L.EQ.K) THEN
CSJA                   ! within FTD
CSJA                   PFSPLIT(5)=PFSPLIT(5)+1
CSJA                ELSE IF (K.EQ.MAX(ID_FTD,ID_SIT).AND.
CSJA     >                   L.EQ.MIN(ID_FTD,ID_SIT)) THEN
CSJA                   ! between FTD and SIT
CSJA                   PFSPLIT(6)=PFSPLIT(6)+1
CSJA                ELSE IF (K.EQ.MAX(ID_TPC,ID_SIT).AND.
CSJA     >                   L.EQ.MIN(ID_TPC,ID_SIT)) THEN
CSJA                   ! between TPC and SIT
CSJA                   PFSPLIT(7)=PFSPLIT(7)+1
CSJA                ELSE IF (K.EQ.MAX(ID_TPC,ID_FTD).AND.
CSJA     >                   L.EQ.MIN(ID_TPC,ID_FTD)) THEN
CSJA                   ! between TPC and FTD
CSJA                   PFSPLIT(8)=PFSPLIT(8)+1
CSJA                ELSE IF (K.EQ.MAX(ID_FCH,ID_FTD).AND.
CSJA     >                   L.EQ.MIN(ID_FCH,ID_FTD)) THEN
CSJA                   ! between FCH and FTD
CSJA                   PFSPLIT(9)=PFSPLIT(9)+1
CSJA                ELSE IF (K.EQ.MAX(ID_FCH,ID_TPC).AND.
CSJA     >                   L.EQ.MIN(ID_FCH,ID_TPC)) THEN
CSJA                   ! between FCH and TPC
CSJA                   PFSPLIT(10)=PFSPLIT(10)+1
CSJA                ELSE IF (K.EQ.ID_TPC.AND.L.EQ.K) THEN
CSJA                   ! within TPC
CSJA                   PFSPLIT(11)=PFSPLIT(11)+1
CSJA                ELSE IF (K.EQ.ID_FCH.AND.L.EQ.K) THEN
CSJA                   ! within FCH
CSJA                   PFSPLIT(12)=PFSPLIT(12)+1
CSJA                ELSE
CSJA                   ! jump over subdetector without contrib. TE failed
CSJA                   PFSPLIT(13)=PFSPLIT(13)+1
CSJA                END IF
CSJA              END IF
CSJA            END DO
CSJA
CSJA          END IF
      END DO
      WRITE(6,'(A40)')  '===================='//
     >     '===================='
CSJA      END IF

 999  RETURN
      END


