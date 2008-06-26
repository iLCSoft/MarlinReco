*  
* $Id: tktrev.F,v 1.14 2008-06-26 09:35:56 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.13  2008/04/17 14:36:29  aplin
* reduced printoutlevel in f77 and moved cout to debug streamlog_out for C++
*
* Revision 1.12  2008/01/24 17:21:50  aplin
* updated the error codes for the case of to many TPC padrows in the GEAR file
* error code returned by TKTREV = 99 this should be caught by LEPTracking
* and an exception thown.
*
* Revision 1.11  2008/01/23 17:22:22  aplin
*
*      The TPC Tracking code has been altered to cope with more than 224 pad rows
*      the max number of pad rows can be changed in padrow.inc by changing
*
*      The maximum number of TPC Padrows that the tracking code
*      can accept must be a multiple of 32
*
*      Note this is not the number of TPC Padrows in the dectector
*      but the maximum allowed because of the predifined array sizes
*      and a set of 32 bit registers
*
*      Currently it is set to 512
*
*      By increasing the value of N32BITREG in padrow.inc you can increase the max
*      number of pad rows
*
* Revision 1.10  2007/12/12 13:40:35  aplin
* Removed the Silicon track matching in LEPTracking. This is now done solely by FullLDCTracking in conjuction with SiliconTracking by Alexei Raspereza.
*
* Revision 1.9  2007/09/05 09:47:29  rasp
* Updated version
*
* Revision 1.7  2006/11/02 11:49:01  rasp
*
* Modified files adopted to FullLDCTracking processor
*
* Revision 1.6  2005/11/03 15:16:14  aplin
* Added the Trackstring creation and the biulding of full Track candiates (TK's) which have passed the Delphi Ambiguity resolver fxambi. The material description of the vtx detector, as for the TPC, is hard coded in setmat. Presently the VTX and SIT resolutions are hard coded in LEPTrackingProcessor. The debug output has been reduced and can be controlled via TKSTDBG etc. in tkinit.F. delsolve contains the delphi ambuguity resolver written in C and is contained in the directory named C. The Tk's are written back into the C++ side in tktrev. The corresponding Tk bank structure analogous to the TE bank structure has been added in tktkbank whilst the access wrapper functions are contained in LEPTracking.
*
* Revision 1.5  2005/08/03 19:14:06  aplin
* removed the combining of TEs to TSs for now
*
* Revision 1.4  2005/08/02 10:33:15  aplin
* included setup of VTX geometry
*
* Revision 1.3  2005/07/29 13:04:28  aplin
* protected expection for now hits in a subsystem
*
* Revision 1.2  2005/07/21 16:37:00  aplin
* The itedat bank has been added with read and write interface.
*
* TKADDH has now been added so that the silicon hits are included as
* TE's.
*
* Added FSFSTR from delsearch to steer the creation of
* TS's
*
* The following include files need to be added for FSFSTR
* include/fsfdebug.inc
* include/fsfreject.inc
* include/fsfflags.inc
*
* FSFSTR is now called by TKTREV although all action is commented out.
*
* Revision 1.1  2005/07/19 17:02:02  aplin
* Removed using namespace from .h files. Commented out remove hit from tkbank, as well as adding methods to give information about the number and positions of subdetector hit collections with in the bank structure. Common block CNTPC no longer needed so commented out in tpchitbank.h, and else where it was accessed. Added protocalls for cfortran.h for the new tkbank methods. Replaced TPCRUN as the main tracking steering routine by TKTREV which itself calls TPCRUN.
*
* Revision 1.2  2004/07/22 15:21:03  aplin
* added tkfunc1.inc
*
* Revision 1.1.1.1  2003/05/23 13:17:51  hvogt
* Brahms V308
*
*  
      INTEGER FUNCTION TKTREV()
**********************************************************************
*
*     subroutine TKTREV (changed to function to provide return value)
*     
*     called by GUTREV
*
*     Kristian Harder, September 1999
*
**********************************************************************
      IMPLICIT NONE
#include "include/bkgmode.inc"
C#include "include/gcflag.inc"
C#include "include/gcbank.inc"
#include "include/fkparm.inc"
#include "include/tkffread.inc"
#include "include/brcntl.inc"
#include "include/fkdebug.inc"
#include "include/tkconst.inc"
#include "include/tkbank.inc"
#include "include/fkexts.inc"
#include "include/detnums.inc"
#include "include/chgtrak.inc"
C#include "include/brtime.inc"
#include "include/tknoisepar.inc"
#include "include/fkddes.inc"
      INTEGER I,J,K,L,IERROR,IKINE, try, ierr
      REAL TSTART,TEND

      REAL cov(15)

      EXTERNAL numofsubdethits
      INTEGER  numofsubdethits
      EXTERNAL subdetfirsthitindex
      INTEGER  subdetfirsthitindex

      INTEGER IDATA(37)
      REAL RDATA(37)

      EQUIVALENCE(IDATA,RDATA)

      integer tkmktkcpp
      external tkmktkcpp

      integer addtetktkcpp
      external addtetktkcpp

      integer writetkitkdatcpp
      external writetkitkdatcpp

#include "include/tkfunc1.inc"
#include "include/tkfunc.inc"
 
C====================================================================




c      IF (TKSTDBG.GT.0) WRITE(6,*) 'TKTREV: event number ',IEVENT

CCC FIXME: the TPC will be intialised here with the correct geometric and physical constants
CCC        this should be done in a more apropriate plays for instance by filling the common 
CCC        block tpcdims on the MARLIN side

C.....Initalise TKTREV as error this will only be set to 0 if TKTREV exits normally      
      TKTREV = 1
      IERROR = 0

      CALL TKINIT
      CALL VTXGEOM
      CALL TPCGEOM

CCC FIXME: SJA: the database of materials is set up in SETMAT for the TPC only all other systems are
CCC             commented out.
C      CALL SETMAT

      CALL TKSETR

c$$$      WRITE(*,*) 'EXTRAPOLATIONS ---> '
c$$$      DO I=1,NEXS
c$$$         WRITE(*,*) I,' ',TYPX(I),' ',RZSURF(I),' ',ZRMIN(I),' ',
c$$$     +              ZRMAX(I)
c$$$      ENDDO
c$$$
c$$$      WRITE(*,*) 'CYLLINDRICAL SHAPES ---> '
c$$$      DO I=1,NCMAT
c$$$         WRITE(*,*) I,' ',RCMAT(I),' ',ZCMIN(I),' ',ZCMAX(I),' ',
c$$$     +             XRLC(I),' ',XELOSC(I)
c$$$      ENDDO
c$$$
c$$$      WRITE(*,*) 'PLANAR SHAPES ---> '
c$$$      DO I=1,NPMAT
c$$$         WRITE(*,*) I,' ',ZPMAT(I),' ',RPMIN(I),' ',RPMAX(I),' ',
c$$$     +             XRLP(I),' ',XELOSP(I)
c$$$      ENDDO



C.....init CHGTRAK common block (still needed in BRAHMS)
c MVL       added also into BRINI in GUTREV
       DO I=1,NMAXCHT
          NTKHITS(I)=0
       END DO

C.....init MC track, TE,TS, and TK arrays
      TKNTRK=0
      NTE=0
      NTS=0
      NTK=0

C.....init exclusion list
      NEXCL=0
      EXCLERR=0

C.....init hit bank index and sub-indices for subdetector IDs
      NTKHIT=0
      DO I=1,3000
        IHPOINT(I)=0
        IHNUMB(I)=0
      END DO

      IHNUMB(ID_TPC)  = numofsubdethits("TPC")

      IF(IHNUMB(ID_TPC).gt.0) IHPOINT(ID_TPC)=subdetfirsthitindex("TPC")
      
      if(TKSTDBG.gt.0) write(*,*) "IHNUMB(ID_TPC)  = " , IHNUMB(ID_TPC) 
      if(TKSTDBG.gt.0) write(*,*) "IHPOINT(ID_TPC) = " , IHPOINT(ID_TPC)


      IHNUMB(ID_VTX)  = numofsubdethits("VTX")

      IF(IHNUMB(ID_VTX).gt.0) IHPOINT(ID_VTX)=subdetfirsthitindex("VTX")
      
      if(TKSTDBG.gt.0) write(*,*) "IHNUMB(ID_VTX)  = " , IHNUMB(ID_VTX) 
      if(TKSTDBG.gt.0) write(*,*) "IHPOINT(ID_VTX) = " , IHPOINT(ID_VTX)

      NTKHIT = IHNUMB(ID_TPC) + IHNUMB(ID_VTX)

      IHNUMB(ID_SIT)  = numofsubdethits("SIT")

      IF(IHNUMB(ID_SIT).gt.0) IHPOINT(ID_SIT)=subdetfirsthitindex("SIT")
      
      if(TKSTDBG.gt.0) write(*,*) "IHNUMB(ID_SIT)  = " , IHNUMB(ID_SIT) 
      if(TKSTDBG.gt.0) write(*,*) "IHPOINT(ID_SIT) = " , IHPOINT(ID_SIT)

      NTKHIT = IHNUMB(ID_TPC) + IHNUMB(ID_VTX) + IHNUMB(ID_SIT)

C.....init TE hit list index
      NHITTE=0

C.....init TS and TK list-of-contained-TE indices
      NTSTEL=0
      NTKTEL=0

C.....clear event veto
      CALL TKVETO(0)

      CALL TIMEX(TSTART)

*     book ZEBRA banks for tracking hits

C This is no longer needed in the Marlin implementation
C      CALL TKBBANK

C.....add noise hits for VTX (only for realistic pattern recognition)
C      IF (IPATRC(1).EQ.1) CALL SIPNOS

C This is no longer needed in the Marlin implementation
C.....fill hits from all tracking subsystems into hit bank and create
C     database of true track parameters
c      CALL TKHVTX
c      CALL TKHFTD
c      CALL TKHITC
c      CALL TKHSIT
c      CALL TKHTPC
c      CALL TKHFCH

*     if anything went wrong - skip the event

c      if ( ieotri .ne. 0 ) then
c        write(*,*) ' Serious error in processing tracking information',
c     &    'run/event = ', idrun, idevt,' -- skip event!'
c        goto 999
c      end if

      CALL TIMEX(TEND)
c      TIME_TKHDET=TIME_TKHDET+(TEND-TSTART)
      TSTART     =TEND

C     adjust the size of the hits storage

c      CALL TKPUSH(1)

C.....skip if no hits found
c      IF (NTKHIT.EQ.0) THEN
*        WRITE(6,*) 'TKTREV: no hits in detector. skipping event.'
c        GOTO 999
c      END IF

      CALL TIMEX(TSTART)

C.....all silicon pattern recognition
c      IF (IPATRC(1).EQ.2) THEN
c ---   dummy patrec
c        CALL LP_SI
c      ELSE
c ---   RJH's silicon detector pattern recognition
c        CALL SIPATR

c ---   also add hit TEs for all silicon detectors

C        CALL TKADDH

c        CALL HCDIR('//PAWC',' ')
c      END IF

      CALL TIMEX(TEND)
c      TIME_LPVTX=TIME_LPVTX+(TEND-TSTART)
      TSTART     =TEND

C.....TPC pattern recognition
c      IF (IPATRC(2).EQ.2) THEN
c ---   dummy patrec
c        CALL LP_TPC
c      ELSE
c ---   GAB's ALEPH-style TPC pattern recognition
c      write(*,*) "---------------------"
c      write(*,*) "About to start TPCRUN"
c      write(*,*) "---------------------"


        CALL TPCRUN(IERROR)
        IF(IERROR.EQ.99) THEN
           GOTO 888
        ENDIF
        

C....  adjust size of hit banks in patrec
c        CALL TKPUSH_TPC(1)
c      END IF

c      CALL HCDIR('//PAWC',' ')

      CALL TIMEX(TEND)
c      TIME_LPTPC=TIME_LPTPC+(TEND-TSTART)
      TSTART     =TEND

C.....FCH pattern recognition
c      IF (IPATRC(3).EQ.2) THEN
C ---    ideal pattern recognition
c         CALL LP_FCH
c      ELSE
C ---    Ivanka Bozovic's FCH pattern recognition
c         CALL FCHPAT
c      END IF

      CALL TIMEX(TEND)
c      TIME_LPFCH=TIME_LPFCH+(TEND-TSTART)
      TSTART     =TEND

C.....adjust size of TE storage banks
c      CALL TKPUSH(2)

C.....cancel tracking if not a single TE found
c      IF (NTE.EQ.0) THEN
c        WRITE(6,*) 'TKTREV: no TEs. skipping analysis.'
c        GOTO 999
c      END IF

C.....mark mutual exclusions of TEs
c      CALL TKMKEX('TE')

      CALL TIMEX(TEND)
c      TIME_TKMKEX =TIME_TKMKEX+(TEND-TSTART)
      TSTART     =TEND

C.....combine TEs to TSs (ambiguous full track candidates) using DELSEARCH.


      CALL FSFSTR

      CALL TIMEX(TEND)
c      TIME_FSFSTR=TIME_FSFSTR+(TEND-TSTART)
      TSTART     =TEND

C.....all TSs are now put into the ambiguity processor
      CALL TKAMBI

      CALL TIMEX(TEND)
c      TIME_TKAMBI=TIME_TKAMBI+(TEND-TSTART)
      TSTART=TEND

C.....was the exclusion list size insufficient?
c      IF (EXCLERR.GT.0) THEN
c        WRITE(6,*) 'TKTREV: ',EXCLERR,' exclusions could not be stored.'
c        WRITE(6,*) '        This should be considered FATAL!'
c        WRITE(6,*) '        Possible effect: Ambiguous track output!'
c        CALL TKVETO(1)
c      END IF

C.....performance analysis
      IF (ITKPERF.GT.0) CALL TKPERF

      CALL TIMEX(TEND)
c      TIME_TKPERF =TIME_TKPERF+(TEND-TSTART)
      TSTART     =TEND

CSJA Send Full TK tracks over to LEPTracking

      do k=1, NTK

         CALL TKREAD('TK',0,k,RDATA,IDATA,IERR)

         do l=1, 15
            cov(l) = (RDATA(l+22))
c      write(*,*) "TDATA(",l+22,") = " , RDATA(l+22)
c      write(*,*) "cov(",l,") = " , cov(l)
         enddo
        
         try = tkmktkcpp(RDATA(1),RDATA(2),RDATA(3),RDATA(4),
     +    RDATA(5),RDATA(6),RDATA(7),RDATA(8),RDATA(9),
     +    RDATA(10),RDATA(11),RDATA(12),RDATA(13),RDATA(14),
     +    RDATA(15),RDATA(16),RDATA(17),RDATA(18),RDATA(19),
     +    RDATA(20),RDATA(21),RDATA(22),cov)     

         try = writetkitkdatcpp(ITKDAT(1,k),1,k)
         try = writetkitkdatcpp(ITKDAT(2,k),2,k)
         try = writetkitkdatcpp(ITKDAT(3,k),3,k)

         do l=1,ITKDAT(2,k)
         try = addtetktkcpp(ITKTEL(ITKDAT(1,k)+l-1),k)      
         enddo

      enddo

C.....fill common block for calorimeter interface/energy flow
C  MVL        Comment this to work with ideal TPC
c      CALL TKCALI

C.....save histograms while running every 50 events, or after each event
C         if we are running with extremely high background level
c      IF ((MOD(ievent,50).EQ.0).OR.(TNNBX.GE.50)) THEN
c        IF (FSFHIST.GT.0.OR.IHIS.GT.0.OR.TKSTHST.GT.0) THEN
c          WRITE(6,*) 'TKTREV: saving histograms in save.hst'
c          CALL HCDIR('//PAWC/TRACKING',' ')
c          CALL HRPUT(0,'save.hst','T')
c          CALL HCDIR('//PAWC',' ')
c        END IF
c        WRITE(6,*) 'TKTREV: writing intermediate run summary.'
c        IF (ITKPERF.GT.0) CALL TKFPERF
c      ENDIF

C.....fill alignment checking histograms whenever in delfit histogramming mode
c      IF (IHIS.GT.0) CALL FSALICHK

C.....finally, some storage usage information for debugging purposes
      IF (TKSTDBG.GT.1) THEN
        WRITE(6,*) 'TKTREV: needed NEXCL=',NEXCL,' of', MXEXCL
      END IF

 888  TKTREV = IERROR

      call flush
 999  RETURN
      
      END


