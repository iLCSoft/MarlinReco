*  
* $Id: fxambi.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
       SUBROUTINE FXAMBI(NACT,INDOK,IFAIL,IFXDBG)
**********************************************************************
*      Steering of processor to remove ambiguities after fit.
*      J.Wickens 25/3/96
*
*      Processor expects TKR to be treated to be DEACTIVE. Signals
*      acceptance of TKR by reactivating the bank.
*
*      IFXDBG - debug flag
*      Output:
*      NACT   - number of good TKR
*      INDOK()- list of identifiers
*      IFAIL  - Bitted word. Bit set if error occurs in this event
* ???? needs updates ????
*               2 Error in TLIST/TGET call for list of TKR
*               3 Overflow of TKR storage
*               4 Failure to get corresponding TS header for TKR
*               5 More than 1 TS header for TKR
*               6 Failure in TLIST call to get TER in TKR
*               7 No TER found or too many TER in TKR
*               8 Failure to get TKR data in TGET call
*               9 Failure in TLOG call to check for excluded TER
*              10 More than MAXEXLI excluded/linked TER for this Event
*              11 Failure to get TER data in TGET call/ label problem
*              12 Failure to reactivate selected TKR bank
*              13 Tlist failure in saving list of accepted TKR
*              18 FXSOLV: Unknown detector code.
*              19 FXSOLV: Wrong Arguments.
*              20 FXSOLV: To many exclusions/links.
*              21 FXSOLV: 1VD removal created EmptyTK.
*
**********************************************************************
*
      IMPLICIT NONE
#include "include/bkgmode.inc"
#include "include/tkconst.inc"
#include "include/fxapar.inc"

      INTEGER IFAIL,IFXDBG,ifxflags,ifxcerr
      INTEGER IWARN,IERR,I,J,K,NDAT,id,newtk,iret,ip

      INTEGER TKNUMB,DELDET,JBIT

      INTEGER NACT,INDOK(MAXTKR)

*
      IFAIL=0
      NACT=0
      CALL VZERO (INDOK,MAXTKR)
* --- ZERO LOCAL VARIABLES
      ITEEXLIUSED=0
      IWARN=0
      NTKR=0
      NRET=0
      CALL VZERO(IDTKR, MAXTKR)
      CALL VZERO(NTERUS,MAXTKR)
      CALL VZERO(MODUSE,MAXTKR)
      CALL VZERO(IRETLI,MAXTKR)
      CALL VZERO(IDTERU,MXTEPT*MAXTKR)
      CALL VZERO(IMODTE,MXTEPT*MAXTKR)
      CALL VZERO(ITKDAT,LTSBAN*MAXTKR)
      CALL VZERO(ILABTE,MXTEPT*MAXTKR)
*
* --- Get number of TKR to fill internal banks. All errors fatal.
*
      NTKR = TKNUMB ('TK',0)
      IF(NTKR.LT.1)THEN
        IF(IFXDBG.GT.0)WRITE(6,1001)
        GO TO 998
      ENDIF
      IF(NTKR.GT.MAXTKR)THEN
        CALL SBIT1(IFAIL,3)
        GO TO 998
      ENDIF
*
* --- Loop over TKR
*
      DO I=1,NTKR
*
        IDTKR(I) = I
*
* --- get TEs which are in that TK
*
        CALL TKLIST ('TK',I,NTERUS(I),IDTERU(1,I),IERR)
        IF(IERR.NE.0)THEN
          CALL SBIT1(IFAIL,6)
          GO TO 998
        ENDIF
        IF(NTERUS(I).LT.1.OR.NTERUS(I).GT.MXTEPT)THEN
          CALL SBIT1(IFAIL,7)
          GO TO 998
        ENDIF
*
* --- get TK data
*
        CALL TKREAD('TK',0,I,TKDAT(1,I),ITKDAT(1,I),IERR)
        IF(IERR.NE.0)THEN
          CALL SBIT1(IFAIL,8)
          GO TO 998
        ENDIF

        DO J=1,NTERUS(I)
*
* --- get logical exclusions
*
          CALL TKLSEX('TE',IDTERU(J,I),MXDAT,NDAT,IDAT,IERR)
          IF(IERR.NE.0)THEN
            CALL SBIT1(IFAIL,9)
          ELSE
            IF(ITEEXLIUSED+NDAT.GT.MAXEXLI)THEN
              CALL SBIT1(IFAIL,10)
            ENDIF
            ITKTEEXLI(3,j,i)=min(NDAT,MAXEXLI-ITEEXLIUSED)
            ITKTEEXLI(4,j,i)=ITEEXLIUSED
            DO K=1,min(NDAT,MAXEXLI-ITEEXLIUSED)
              ITEEXLIIDS(ITEEXLIUSED+K)=IDAT(K)
            ENDDO
            ITEEXLIUSED=ITEEXLIUSED+min(NDAT,MAXEXLI-ITEEXLIUSED)
          ENDIF
*
* --- get TE data
*
          CALL TKREAD('TE',0,IDTERU(J,I),DAT,IDAT,IERR)
          IF(IERR.NE.0)THEN
            CALL SBIT1(IFAIL,11)
          ELSE
* --- transform into DELPHI like identifiers !!!
            IMODTE(J,I) = DELDET(IDAT(1))
            CALL SBIT1(MODUSE(I),IMODTE(J,I))
* --- get label
            CALL TKGLAB('TE',IDTERU(J,I),ILABTE(J,I),IERR)
            IF (IERR.NE.0)THEN
              CALL SBIT1(IFAIL,11)
            ENDIF
          ENDIF
*
        ENDDO
      ENDDO
* --- debug output
      IF(IFXDBG.GT.1)THEN
        WRITE(6,2000)
        CALL FXDUMP(IFXDBG)
      ENDIF
*
c(KH)--begin
      ifxflags=0
c(KH)--end
      call sbyt(ifxdbg,ifxflags,1,2) ! Copy debug level
*
* --- Resolve ambiguities
*
      CALL FXSOLV(NTKR,IDTKR,NTERUS,IDTERU,IMODTE,ITKTEEXLI,ITEEXLIIDS,
     /            ILABTE,ITKDAT,MXTEPT,MAXTKR,LTSBAN,IRETLI,
     /            NRET,IFXFLAGS,IFXCERR)

      call cbyt(IFXCERR, 1,IWARN,1,16) ! Bits  1 to 16 are warnings/infos
      call cbyt(IFXCERR,17,IFAIL,18,5) ! Bits 17 to 32 are errors

      IF(IFXDBG.GT.0)THEN
        print*,'FXAMBI: IFXCERR:',ifxcerr,' IWARN:',iwarn,
     /         ' IFAIL:',ifail
        IF(IFXDBG.GT.1)THEN
          WRITE(6,2001)
          CALL FXDUMP(IFXDBG)
        ENDIF
      ENDIF
*
* --- Reactivate accepted TKR via reactivation of associated TS
*     for cases where the track is not modified.
*
      DO I=1,NRET
        ID=IRETLI(I)
        IF(IDTKR(ID).GT.0)THEN
          NACT=NACT+1
          INDOK(NACT)=IDTKR(ID)
        ENDIF
      ENDDO
*
* --- here we create new TKs !!! ????
      CALL FXAREC(NEWTK,INDOK(NACT+1),IFXDBG,IERR)
      IF(IERR.NE.0)CALL SBYT(IERR,IFAIL,14,4)
      NACT=NACT+NEWTK

      IF(IFXDBG.GT.1)THEN
        WRITE(6,2002)NACT,(INDOK(K),K=1,NACT)
      ENDIF
*
*      Statistics
*
      IP = 1
      DO IRET=1,NRET
        I=IRETLI(IRET)
        ISTAT(0)=ISTAT(0)+1
        DO J = 1,NWSTAT
          ISTAT(J)=ISTAT(J)+JBIT(MODUSE(I),J)
        ENDDO
      ENDDO
*
 999  CONTINUE
      IF(IWARN.NE.0)THEN
        IF (IFXDBG.GT.0) WRITE(6,1002)IWARN
        DO I=1,NWCODE
          IF(JBIT(IWARN,I).EQ.1)KSUMWRN(I)=KSUMWRN(I)+1
        ENDDO
      ENDIF
      IF(IFAIL.NE.0)THEN
        IF (IFXDBG.GT.0) WRITE(6,1000)IFAIL
        DO 995 I=1,NECODE
          IF(JBIT(IFAIL,I).EQ.1)KSUMER(I)=KSUMER(I)+1
 995    CONTINUE
      ENDIF
*
 1000 FORMAT(1X,'FXAMBI: IFAIL= ',I8)
 1001 FORMAT(1X,'FXAMBI: NO TRACKS ON INPUT')
 1002 FORMAT(1X,'FXAMBI: IWARN= ',I8)
 2000 FORMAT(1X,' '/1X,'FXDUMP:Before call to FXSOLV')
 2001 FORMAT(1X,' '/1X,'FXDUMP:After call to FXSOLV')
 2002 FORMAT(1X,' '/1X,'Tracks retained : Total ',I5,
     /     ' With indices:'/(1X,30I4))
 998  RETURN
      END

