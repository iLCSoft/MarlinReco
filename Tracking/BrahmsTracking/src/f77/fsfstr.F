*  
* $Id: fsfstr.F,v 1.5 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.3  2005/08/04 12:54:51  aplin
* *** empty log message ***
*
* Revision 1.2  2005/08/02 10:19:51  aplin
* corrected debug output
*
* Revision 1.1  2005/07/21 16:37:00  aplin
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
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      SUBROUTINE FSFSTR
*******************************************************************************
*
*
*               steering of track search
*
*                                                                  Daniel Wicke
*******************************************************************************
      IMPLICIT NONE
C +cde,UXCEVT.
#include "include/fsfclun.inc"
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfreject.inc"
#include "include/fsfflags.inc"
#include "include/turnon.inc"
#include "include/fsfmc.inc"

      integer i,imethod
      integer nid,none,isave,data_labl
      integer method(detmax)
      real tstart,tend

      parameter (none=4711) !Must be different from all possible methods for
                            !pairing

      if(fsfdbg.gt.0)then
        write(6,*) '****************************************'
        print*,'FSFSTR:   starting new event'
        write(6,*) '****************************************'
      endif

c --- copy beam spot into primary vertex for this event
      call fsfprivtx

c --- Read event data from TKBANK to local arrays
      call fsfinp

c --- Analyse event :
c

c --- list some information about the TE input  -KH-
      if (fsfdbg.gt.1) then
        write(6,*) '-----------------------------'
        do i=1,detmax
          do nid=1,teindexused(i)
            write(6,*) 'FSFSTR: TE ',nid,' of ',detname(i),': label ',
     >                 telabl(nid,i)
          end do
        end do
        write(6,*) '-----------------------------'
      end if

c --- set method. 0=no matching, 1=real matching, -1=matching using MC info
c --- no matching is default
      call vzero(method,detmax)

      if (.not.mcmerge) then
        if(fsfdbg.gt.0)then
          print*,'FSFSTR: normal extrapolations for making pairs'
        endif
        imethod=1
        data_labl = FLAG_DATA
      else
        if(fsfdbg.gt.0)then
          print*,'FSFSTR: used labels for making pairs'
        endif
        imethod=-1
        data_labl = FLAG_LABL
      endif

c --- never extrapolate onto strings
      method(i_SCAN1)=0
      method(i_SCAN2)=0
      method(i_SCAN3)=0
      method(i_SCAN4)=0


c --- ---------------------------------------
c --- --- barrel search starting from VTX ---
c --- ---------------------------------------

      call vzero(method,detmax)
c --- set extrapolation methods
      method(i_FTD) =0
      method(i_ITC) =0
      method(i_SIT1)=0
      method(i_SIT2)=0
      method(i_TPC) =imethod
      method(i_FCH) =imethod

c --- init timing
      CALL TIMEX(TSTART)

c --- perform extrapolation
      if(fsfdbg.gt.0)then
        write(6,*) '****************************************'
        print*,'FSFSTR: Barrel Track Search from VTX'
        write(6,*) '****************************************'
      endif
      method(i_VTX) = none      ! this is the pivot
      call fsfana3 (i_VTX,method,
     1     data_labl+FLAG_GIVESINGLTE+FLAG_NOPLINV+FLAG_USEUSED+
     2     FLAG_CUT100)

c --- read timing
      CALL TIMEX(TEND)
      TIME(i_VTX) = TIME(i_VTX)+(TEND-TSTART)
      TSTART      = TEND

c --- remove duplicate and sub-TS
      call fsfnsub(FLAG_DBL+FLAG_SUB)

c --- list of TS as found up to now
      if (fsfdbg.gt.1) call fsftsdump

c --- ---------------------------------------
c --- --- barrel search starting from TPC ---
c --- ---------------------------------------

      call vzero(method,detmax)
c --- set extrapolation methods
      method(i_VTX) =imethod
      method(i_FTD) =0
      method(i_ITC) =0
      method(i_SIT1)=0
      method(i_SIT2)=0
      method(i_FCH)=imethod

c --- init timing
      CALL TIMEX(TSTART)

c --- perform extrapolation
      if(fsfdbg.gt.0)then
        write(6,*) '****************************************'
        print*,'FSFSTR: Barrel Track Search from TPC'
        write(6,*) '****************************************'
      endif
      method(i_TPC) = none      ! this is the pivot
      call fsfana3 (i_TPC,method,
CC     1     data_labl+FLAG_GIVESINGLTE+FLAG_NOPLINV+FLAG_USEUSED+
     1     data_labl+FLAG_GIVESINGLTE+FLAG_NOPLINV+
     2     FLAG_CUT100)

c --- read timing
      CALL TIMEX(TEND)
      TIME(i_TPC) = TIME(i_TPC)+(TEND-TSTART)
      TSTART      = TEND

c --- remove duplicate and sub-TS
      call fsfnsub(FLAG_DBL+FLAG_SUB)

c --- list of TS as found up to now
      if (fsfdbg.gt.1) call fsftsdump

c --- -------------------------------------------------------
c --- --- now cleanup and iterate searches in SIT,TPC,FCH ---
c --- -------------------------------------------------------

      call vzero(method,detmax)
c --- set extrapolation methods
      method(i_VTX) =0
      method(i_FTD) =0
      method(i_ITC) =0
      method(i_SIT1)=imethod
      method(i_SIT2)=imethod
      method(i_TPC) =imethod
      method(i_FCH) =imethod

c --- init timing
      CALL TIMEX(TSTART)

c --- now start new search to find all TEs not yet found before
      data_labl = FLAG_DATA
      if(fsfdbg.gt.0)then
        write(6,*) '****************************************'
        print*,'FSFSTR: STRING extrapolation to SIT/FCH '
        write(6,*) '****************************************'
      endif

c --- now create new pseudo TEs from all TS we have found.
      call fsfrest(i_SCAN1)

      method(i_SCAN1) = none  ! this is the pivot
      call fsfana3 (i_SCAN1,method,
     1                data_labl+FLAG_SCANDET+FLAG_NOPLINV+
     2                FLAG_CUT100)

c --- read timing
      CALL TIMEX(TEND)
      TIME(i_SCAN1) = TIME(i_SCAN1)+(TEND-TSTART)
      TSTART        = TEND

c --- remove double and sub-strings.
      if(fsfdbg.gt.0)then
        print*,'FSFSTR: remove substrings'
      endif
      call fsfnsub(FLAG_DBL+FLAG_SUB)

c --- ----------------------------------------------------------
c --- --- now cleanup and iterate searches in FTD, outer VTX ---
c --- ----------------------------------------------------------

      call vzero(method,detmax)
c --- set extrapolation methods
      method(i_VTX) =0
      method(i_VTX1)=0
      method(i_VTX2)=0
      method(i_VTX3)=imethod
      method(i_VTX4)=imethod
      method(i_VTX5)=imethod
      method(i_FTD) =0
      method(i_FTD1)=imethod
      method(i_FTD2)=imethod
      method(i_FTD3)=imethod
      method(i_FTD4)=imethod
      method(i_FTD5)=imethod
      method(i_FTD6)=imethod
      method(i_FTD7)=imethod
      method(i_ITC) =0
      method(i_SIT1)=0
      method(i_SIT2)=0
      method(i_TPC) =0
      method(i_FCH)=0

c --- init timing
      CALL TIMEX(TSTART)

c --- now start new search to find all TEs not yet found before
      data_labl = FLAG_DATA
      if(fsfdbg.gt.0)then
        write(6,*) '********************************************'
        print*,'FSFSTR: STRING extrapolation to FTD/outer VTX'
        write(6,*) '********************************************'
      endif

c --- now create new pseudo TEs from all TS we have found.
      call fsfrest(i_SCAN2)

      method(i_SCAN2) = none  ! this is the pivot
      call fsfana3 (i_SCAN2,method,
     1                data_labl+FLAG_SCANDET+FLAG_NOPLINV+
     2                FLAG_CUT100)
      method(i_SCAN2) = 0

c --- read timing
      CALL TIMEX(TEND)
      TIME(i_SCAN2) = TIME(i_SCAN2)+(TEND-TSTART)
      TSTART        = TEND

c --- remove double and sub-strings.
      if(fsfdbg.gt.0)then
        print*,'FSFSTR: remove substrings'
      endif
      call fsfnsub(FLAG_DBL+FLAG_SUB)

c --- -----------------------------------------------------
c --- --- now cleanup and iterate searches in inner VTX ---
c --- -----------------------------------------------------

      call vzero(method,detmax)
c --- set extrapolation methods
      method(i_VTX2)=imethod

c --- init timing
      CALL TIMEX(TSTART)

c --- now start new search to find all TEs not yet found before
      data_labl = FLAG_DATA
      if(fsfdbg.gt.0)then
        write(6,*) '*****************************************'
        print*,'FSFSTR: STRING extrapolation to inner VTX'
        write(6,*) '*****************************************'
      endif

c --- now create new pseudo TEs from all TS we have found.
      call fsfrest(i_SCAN3)

      method(i_SCAN3) = none  ! this is the pivot
      call fsfana3 (i_SCAN3,method,
     1                data_labl+FLAG_SCANDET+FLAG_NOPLINV+
     2                FLAG_CUT100)
      method(i_SCAN3) = 0

c --- read timing
      CALL TIMEX(TEND)
      TIME(i_SCAN3) = TIME(i_SCAN3)+(TEND-TSTART)
      TSTART        = TEND

c --- remove double and sub-strings.
      if(fsfdbg.gt.0)then
        print*,'FSFSTR: remove substrings'
      endif
      call fsfnsub(FLAG_DBL+FLAG_SUB)

c --- -----------------------------------------------------
c --- --- now cleanup and iterate searches in inner VTX ---
c --- -----------------------------------------------------

      call vzero(method,detmax)
c --- set extrapolation methods
      method(i_VTX1)=imethod

c --- init timing
      CALL TIMEX(TSTART)

c --- now start new search to find all TEs not yet found before
      data_labl = FLAG_DATA
      if(fsfdbg.gt.0)then
        write(6,*) '*****************************************'
        print*,'FSFSTR: STRING extrapolation to inner VTX'
        write(6,*) '*****************************************'
      endif

c --- now create new pseudo TEs from all TS we have found.
      call fsfrest(i_SCAN4)

      method(i_SCAN4) = none  ! this is the pivot
      call fsfana3 (i_SCAN4,method,
     1                data_labl+FLAG_SCANDET+FLAG_NOPLINV+
     2                FLAG_CUT100)
      method(i_SCAN4) = 0

c --- read timing
      CALL TIMEX(TEND)
      TIME(i_SCAN4) = TIME(i_SCAN4)+(TEND-TSTART)
      TSTART        = TEND

c --- remove double and sub-strings.
      if(fsfdbg.gt.0)then
        print*,'FSFSTR: remove substrings'
      endif
      call fsfnsub(FLAG_DBL+FLAG_SUB)

c --- ----------------------------------------
c --- --- refit without error scalings     ---
c --- ----------------------------------------

c --- final fit without any error scaling
      if(fsfdbg.gt.0)then
        write(6,*) '****************************************'
        print*,'FSFSTR: refit without error scaling to clean up'
        write(6,*) '****************************************'
      endif
      call fsfrefit

c --- now do the real test for funny TS
      if(fsfdbg.gt.0)then
        write(6,*) '****************************************'
        print*,'FSFSTR: remove double and test funny TS in list'
        write(6,*) '****************************************'
      endif
      call fsfnsub(FLAG_DBL+FLAG_SUB)

c --- ----------------------------------------
c --- --- now remove funny TS combinations ---
c --- ----------------------------------------

      do i = 1,tsindexused
        if(tsindex(i).ne.0)then
          call fsftts(tsindex(i),badcomb,badsingle)
          if(tsindex(i).eq.0)then
            if(fsfdbg.gt.0)then
              print*,'FSFSTR: funny TS not given !'
            endif
          endif
        endif
      enddo
c --- remove double and sub-strings.
      if(fsfdbg.gt.0)then
        print*,'FSFSTR: remove substrings'
      endif
      call fsfnsub(FLAG_DBL+FLAG_SUB)

c --- ----------------------------------------
c --- --- save final TS candidates         ---
c --- ----------------------------------------

c --- Save all the TS and make TK of what we have made.
      call fsfsts(FLAG_GIVESINGLTE)
      if(fsfdbg.gt.0) then
        write(6,*) '****************************************'
        print*,'FSFSTR: event done !'
        write(6,*) '****************************************'
      endif
      return
      END


