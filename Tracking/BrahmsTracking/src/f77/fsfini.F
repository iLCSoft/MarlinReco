*  
* $Id: fsfini.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      SUBROUTINE FSFINI
*******************************************************************************
*
*
*     initialistion of track search
*
*     WRITTEN BY                                    Daniel Wicke
*     CHANGES BY                                    Kristian Harder, August 99
*******************************************************************************
      implicit none
#include "include/fsfclun.inc"
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfxerrs.inc"
#include "include/fsfarrays.inc"
#include "include/fsfpi.inc"
#include "include/fsfsim.inc"
#include "include/turnon.inc"
      integer i,j

* --- print version number from PATCHY title
      WRITE (IPRNT,777)
  777 FORMAT (/' Track search compiled from PAM file ',
     +20HDelsearch      1.10   ,
     +/)

      if(fsfdbg.gt.0)then
        print*,'FSFINI: searches starting at the following detectors:'
        write(6,'(1X,A,$)') 'FSFINI: '
        do i=1,maxpd
          write(6,'(A,$)') detname(i),' / '
        end do
        write(6,*)
      endif

* --- force MC mode for the time being (-KH-)
      delsim=.TRUE.

#if defined(FSTEST)
      call fstest('INIT')
#endif

      maxspaceused = 0
c
c --- Error scaling for TEs in trackfit to soften outlayer in search.
c
      fiterrscale = 2.
c      write(6,*) 'fsfini: warning: using large error scaling'
c      fiterrscale = 10.
c
c --- Error scaling for TPC TEs in trackfit to soften outlayer in search.
c
      tpcerrscale = 2.

c --- initialize all errors to be sure not to have missed one
      DO i=1,detmax
        DO j=1,detmax
         xcx     (i,j) = 1e30
         xcy     (i,j) = 1e30
         xcRPhiR (i,j) = 1e30
         xcz     (i,j) = 1e30
         xctheta (i,j) = 1e30
         xcphi   (i,j) = 1e30
         xcthest (i,j) = 1e30
         xcphist (i,j) = 1e30
         plx     (i,j) = 1e30
         ply     (i,j) = 1e30
         plRPhiR (i,j) = 1e30
         plz     (i,j) = 1e30
         pltheta (i,j) = 1e30
         plphi   (i,j) = 1e30
         plthest (i,j) = 1e30
         plphist (i,j) = 1e30
        ENDDO
      ENDDO

c -----------------------------------------------------------------------------
c     These are the errors, that are added to the errors reported
c     in the TE to account for the inaccurateness of the extrapolation
c     (e.g. due to energy loss, when using helix extrapolation)
c
c     Plphi will be used to compare detectors with the fitresult.
c -----------------------------------------------------------------------------

c -----------------------------------------------------------------------------
c     First for extrapolation from TPC to elsewhere. (24.11.99)
c     no polar inversion needed
c -----------------------------------------------------------------------------

c --- TPC->VTX
         xcRPhiR (i_TPC,i_VTX) = 0.05
         xcz     (i_TPC,i_VTX) = 1.0
         xctheta (i_TPC,i_VTX) = 0.04
         xcphi   (i_TPC,i_VTX) = 0.03

c --- TPC->SIT1 (inner SIT)
c         xcRPhiR (i_TPC,i_SIT1) = 0.05
c         xcz     (i_TPC,i_SIT1) = 1.0
c --- TPC->SIT2 (outer SIT)
c         xcRPhiR (i_TPC,i_SIT2) = 0.02
c         xcz     (i_TPC,i_SIT2) = 0.7

c --- TPC->FCH
         xcx     (i_TPC,i_FCH) = 1.0
         xcy     (i_TPC,i_FCH) = 1.0
         xctheta (i_TPC,i_FCH) = 0.015


c -----------------------------------------------------------------------------
C     for extrapolation from VTX to elsewhere. (24.11.99)
c -----------------------------------------------------------------------------

c --- VTX->SIT1
c         xcRPhiR (i_VTX,i_SIT1) = 0.015
c         xcz     (i_VTX,i_SIT1) = 0.5

c --- VTX->SIT2
c         xcRPhiR (i_VTX,i_SIT2) = 0.02
c         xcz     (i_VTX,i_SIT2) = 0.5

c --- VTX->TPC
         xcRPhiR (i_VTX,i_TPC) = 0.015
         xcz     (i_VTX,i_TPC) = 0.4
         xctheta (i_VTX,i_TPC) = 0.008
         xcphi   (i_VTX,i_TPC) = 0.02

c --- VTX->FCH
         xcx (i_VTX,i_FCH)     = 0.5
         xcy (i_VTX,i_FCH)     = 0.6
         xctheta (i_VTX,i_FCH) = 0.016
         xcthest (i_VTX,i_FCH) = 0.01
         xcphist (i_VTX,i_FCH) = 0.05

c --- VTX->FTD
c         xcx     (i_VTX,i_FTD) = 0.2
c         xcy     (i_VTX,i_FTD) = 0.2
c         xctheta (i_VTX,i_FTD) = 0.005
c         xcphi   (i_VTX,i_FTD) = 0.01

c -----------------------------------------------------------------------------
c     extrapolation from FTD to elsewhere. (24.11.99)
c -----------------------------------------------------------------------------

c --- FTD->VTX
c         xcRPhiR (i_FTD,i_VTX)  = 0.02
c         xcz     (i_FTD,i_VTX)  = 0.2

c --- FTD->TPC
c         xcRPhiR (i_FTD,i_TPC) = 0.015
c         xcz     (i_FTD,i_TPC) = 1.0
c         xctheta (i_FTD,i_TPC) = 0.01
c         xcphi   (i_FTD,i_TPC) = 0.02

c -----------------------------------------------------------------------------
c     Extrapolation of strings (1) to elsewhere. (17.02.00)
c     no polar inversion needed
c -----------------------------------------------------------------------------

c --- SCAN1->SIT1 (inner SIT)
         xcRPhiR (i_SCAN1,i_SIT1) = 0.02
C(KH)    xcz     (i_SCAN1,i_SIT1) = 0.5
         xcz     (i_SCAN1,i_SIT1) = 1.0
c --- SCAN1->SIT2 (outer SIT)
         xcRPhiR (i_SCAN1,i_SIT2) = 0.02
C(KH)    xcz     (i_SCAN1,i_SIT2) = 0.5
         xcz     (i_SCAN1,i_SIT2) = 1.0

c --- SCAN1->TPC
         xcRPhiR (i_SCAN1,i_TPC) = 0.05
         xcz     (i_SCAN1,i_TPC) = 2.
         xctheta (i_SCAN1,i_TPC) = 0.05
         xcphi   (i_SCAN1,i_TPC) = 0.15

c --- SCAN1->FCH
         xcx (i_SCAN1,i_FCH)     = 1.0
         xcy (i_SCAN1,i_FCH)     = 1.0
         xctheta (i_SCAN1,i_FCH) = 0.04


c -----------------------------------------------------------------------------
c     Extrapolation of strings (2) to elsewhere. (25.11.99)
c     no polar inversion needed
c -----------------------------------------------------------------------------

c --- SCAN2->VTX1
c         xcRPhiR (i_SCAN2,i_VTX1) = 0.025
c         xcz     (i_SCAN2,i_VTX1) = 0.25
c --- SCAN2->VTX2
c         xcRPhiR (i_SCAN2,i_VTX2) = 0.02
c         xcz     (i_SCAN2,i_VTX2) = 0.25
c --- SCAN2->VTX3
         xcRPhiR (i_SCAN2,i_VTX3) = 0.006
         xcz     (i_SCAN2,i_VTX3) = 0.15
c --- SCAN2->VTX4
         xcRPhiR (i_SCAN2,i_VTX4) = 0.006
         xcz     (i_SCAN2,i_VTX4) = 0.15
c --- SCAN2->VTX5
         xcRPhiR (i_SCAN2,i_VTX5) = 0.006
         xcz     (i_SCAN2,i_VTX5) = 0.15

c --- SCAN2->FTD1
         IF (LVXDC) THEN
           xcx    (i_SCAN2,i_FTD1) = 0.15
           xcy    (i_SCAN2,i_FTD1) = 0.15
           xcthest(i_SCAN2,i_FTD1) = 0.015
           xcphist(i_SCAN2,i_FTD1) = 0.05
         END IF
c --- SCAN2->FTD2
         xcx    (i_SCAN2,i_FTD2) = 0.15
         xcy    (i_SCAN2,i_FTD2) = 0.15
         xcthest(i_SCAN2,i_FTD2) = 0.015
         xcphist(i_SCAN2,i_FTD2) = 0.05
c --- SCAN2->FTD3
         xcx (i_SCAN2,i_FTD3)    = 0.2
         xcy (i_SCAN2,i_FTD3)    = 0.2
         xcthest(i_SCAN2,i_FTD3) = 0.015
         xcphist(i_SCAN2,i_FTD3) = 0.08
c --- SCAN2->FTD4
         xcx (i_SCAN2,i_FTD4) = 0.25
         xcy (i_SCAN2,i_FTD4) = 0.25
         xcthest(i_SCAN2,i_FTD4) = 0.02
         xcphist(i_SCAN2,i_FTD4) = 0.1
c --- SCAN2->FTD5
         xcx (i_SCAN2,i_FTD5) = 0.5
         xcy (i_SCAN2,i_FTD5) = 0.5
         xcthest(i_SCAN2,i_FTD5) = 0.03
         xcphist(i_SCAN2,i_FTD5) = 0.15
c --- SCAN2->FTD6
         xcx (i_SCAN2,i_FTD6) = 0.7
         xcy (i_SCAN2,i_FTD6) = 0.7
         xcthest(i_SCAN2,i_FTD6) = 0.05
         xcphist(i_SCAN2,i_FTD6) = 0.2
c --- SCAN2->FTD7
         xcx (i_SCAN2,i_FTD7) = 0.7
         xcy (i_SCAN2,i_FTD7) = 0.7
         xcthest(i_SCAN2,i_FTD6) = 0.06
         xcphist(i_SCAN2,i_FTD6) = 0.2

c -----------------------------------------------------------------------------
c     Extrapolation of strings (3) to elsewhere. (06.02.01)
c     no polar inversion needed
c -----------------------------------------------------------------------------

c --- SCAN3->VTX2
         IF (LVXDP) THEN
           xcRPhiR (i_SCAN3,i_VTX2) = 0.02
           xcz     (i_SCAN3,i_VTX2) = 0.25
         ELSE
           xcRPhiR (i_SCAN3,i_VTX2) = 0.01
           xcz     (i_SCAN3,i_VTX2) = 0.1
         END IF

c -----------------------------------------------------------------------------
c     Extrapolation of strings (4) to elsewhere. (06.02.01)
c     no polar inversion needed
c -----------------------------------------------------------------------------

c --- SCAN4->VTX1
         IF (LVXDP) THEN
           xcRPhiR (i_SCAN4,i_VTX1) = 0.025
           xcz     (i_SCAN4,i_VTX1) = 0.25
         ELSE
           xcRPhiR (i_SCAN4,i_VTX1) = 0.005  ! was 0.016 for the LC note
           xcz     (i_SCAN4,i_VTX1) = 0.1
         END IF


      END

