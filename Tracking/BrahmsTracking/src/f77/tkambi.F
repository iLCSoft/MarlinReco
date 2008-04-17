*  
* $Id: tkambi.F,v 1.2 2008-04-17 14:36:29 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1  2005/11/03 15:16:14  aplin
* Added the Trackstring creation and the biulding of full Track candiates (TK's) which have passed the Delphi Ambiguity resolver fxambi. The material description of the vtx detector, as for the TPC, is hard coded in setmat. Presently the VTX and SIT resolutions are hard coded in LEPTrackingProcessor. The debug output has been reduced and can be controlled via TKSTDBG etc. in tkinit.F. delsolve contains the delphi ambuguity resolver written in C and is contained in the directory named C. The Tk's are written back into the C++ side in tktrev. The corresponding Tk bank structure analogous to the TE bank structure has been added in tktkbank whilst the access wrapper functions are contained in LEPTracking.
*
* Revision 1.2  2004/07/22 15:16:36  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TKAMBI
C****************************************************************************
c* Calling routine for DELAMBI                                              *
c*                                                                          *
c* Kristian Harder, September 1999                                          *
c*                                                                          *
c*                                                Kristian.Harder@cern.ch   *
c****************************************************************************
      IMPLICIT NONE
#include "include/bkgmode.inc"
#include "include/tkffread.inc"
#include "include/tkconst.inc"
#include "include/fkparm.inc"
#include "include/tkbank.inc"
c#include "gcbank.inc"

      INTEGER I,IFAIL,IUSE(NTKMX)
      INTEGER NACT,INDOK(NTKMX)
#include "include/tkfunc1.inc"
#include "include/tkfunc.inc"

c      write(*,*) "---------------------"
c      write(*,*) "About to start TKAMBI" 
c      write(*,*) "---------------------"

      NACT = 0
      CALL VZERO (INDOK,NTKMX)

c --- copy the complete list of TEs in TSs into the same array for TKs.
c     this is possible for this simple ambiguity resolver, since it does not
c     change the TE composition of the TK with respect to the TS. It just
c     deletes some TSs, with the consequence that some entries in the newly
c     created (TEs in TK)-list are not referenced by any TK. this does not
c     bother us, though, since the array size is fixed anyway.
      CALL UCOPY(ITSTEL,ITKTEL,NTSTEL)
      NTKTEL=NTSTEL

c --- copy all used TSs into TK array
      DO I=1,NTS
        NTK=NTK+1
        IF (NTK.GT.NTKMX) THEN
          WRITE(6,*) 'TKAMBI: NTKMX TOO SMALL'
          NTK=NTKMX
          CALL TKVETO(1)
        END IF
c ---     copy TS structure itself and convert measurement code
        CALL UCOPY(RTS(1,I),RTK(1,NTK),MIN0(MXTS,MXTK))
        ITK(3,NTK)=IAND(ITS(3,I),2)/2
c ---     copy ITSDAT
        ITKDAT(1,NTK)=ITSDAT(1,I)
        ITKDAT(2,NTK)=ITSDAT(2,I)
        ITKDAT(3,NTK)=ITSDAT(5,I)
      END DO


c      write(*,*) "-------------------------"
c      write(*,*) "About to call FXAMBI" 
c      write(*,*) "-------------------------"


c --- now call ambiguity processor !!!
      CALL FXAMBI(NACT,INDOK,IFAIL,DELADBG)


      IF (IFAIL.NE.0.AND.TKSTDBG.GT.0) THEN
        PRINT *,'TKAMBI: error from FXAMBI is IERR: ',IFAIL
      ENDIF

c --- clean up TK array
      IF (NACT.LT.NTK) THEN
        IF (TKSTDBG.GT.0) WRITE(6,*) 'TKAMBI: purge TK bank.'

        CALL VZERO(IUSE,NTK)
        DO I=1,NACT
          IF (TKSTDBG.GT.1)
     >      WRITE(6,*) 'TKAMBI: mark TK ', INDOK(I),' as used.'
          IUSE(INDOK(I))=1
        END DO

        I=1
 10     IF (IUSE(I).GT.0) THEN
          IF (TKSTDBG.GT.1) WRITE(6,*) 'TKAMBI: TK ',I,' is used.'
          ! check next TK
          I=I+1
          IF (I.LE.NTK) GOTO 10
        ELSE
          IF (TKSTDBG.GT.1) WRITE(6,*) 'TKAMBI: TK ',I,' replaced.'
          ! kill this TK by copying last TK in array to current position
          CALL UCOPY(RTK(1,NTK),RTK(1,I),MXTK)
          CALL UCOPY(ITKDAT(1,NTK),ITKDAT(1,I),3)
          IUSE(I)=IUSE(NTK)
          NTK=NTK-1
          ! now we have to check this (formerly last, now current) TK
          IF (I.LE.NTK) GOTO 10
        END IF
      END IF

      RETURN
      END


