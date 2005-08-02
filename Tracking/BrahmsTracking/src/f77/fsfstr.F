*  
* $Id: fsfstr.F,v 1.2 2005-08-02 10:19:51 aplin Exp $
*  
* $Log: not supported by cvs2svn $
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
c      call fsfinp

c --- Analyse event :
CSJAc
CSJA
CSJAc --- list some information about the TE input  -KH-
CSJA      if (fsfdbg.gt.1) then
CSJA        write(6,*) '-----------------------------'
CSJA        do i=1,detmax
CSJA          do nid=1,teindexused(i)
CSJA            write(6,*) 'FSFSTR: TE ',nid,' of ',detname(i),': label ',
CSJA     >                 telabl(nid,i)
CSJA          end do
CSJA        end do
CSJA        write(6,*) '-----------------------------'
CSJA      end if
CSJA
CSJA --- set method. 0=no matching, 1=real matching, -1=matching using MC info
CSJA --- no matching is default
CSJA      call vzero(method,detmax)
CSJA
CSJA      if (.not.mcmerge) then
CSJA        if(fsfdbg.gt.0)then
CSJA          print*,'FSFSTR: normal extrapolations for making pairs'
CSJA        endif
CSJA        imethod=1
CSJA        data_labl = FLAG_DATA
CSJA      else
CSJA        if(fsfdbg.gt.0)then
CSJA          print*,'FSFSTR: used labels for making pairs'
CSJA        endif
CSJA        imethod=-1
CSJA        data_labl = FLAG_LABL
CSJA      endif
CSJA
CSJAc --- never extrapolate onto strings
CSJA      method(i_SCAN1)=0
CSJA      method(i_SCAN2)=0
CSJA      method(i_SCAN3)=0
CSJA      method(i_SCAN4)=0
CSJA
CSJA
CSJAc --- ---------------------------------------
CSJAc --- --- barrel search starting from VTX ---
CSJAc --- ---------------------------------------
CSJA
CSJA      call vzero(method,detmax)
CSJAc --- set extrapolation methods
CSJA      method(i_FTD) =0
CSJA      method(i_ITC) =0
CSJA      method(i_SIT1)=0
CSJA      method(i_SIT2)=0
CSJA      method(i_TPC) =imethod
CSJA      method(i_FCH) =imethod
CSJA
CSJAc --- init timing
CSJA      CALL TIMEX(TSTART)
CSJA
CSJAc --- perform extrapolation
CSJA      if(fsfdbg.gt.0)then
CSJA        write(6,*) '****************************************'
CSJA        print*,'FSFSTR: Barrel Track Search from VTX'
CSJA        write(6,*) '****************************************'
CSJA      endif
CSJA      method(i_VTX) = none      ! this is the pivot
CSJA      call fsfana3 (i_VTX,method,
CSJA     1     data_labl+FLAG_GIVESINGLTE+FLAG_NOPLINV+FLAG_USEUSED+
CSJA     2     FLAG_CUT100)
CSJA
CSJAc --- read timing
CSJA      CALL TIMEX(TEND)
CSJA      TIME(i_VTX) = TIME(i_VTX)+(TEND-TSTART)
CSJA      TSTART      = TEND
CSJA
CSJAc --- remove duplicate and sub-TS
CSJA      call fsfnsub(FLAG_DBL+FLAG_SUB)
CSJA
CSJAc --- list of TS as found up to now
CSJA      if (fsfdbg.gt.1) call fsftsdump
CSJA
CSJAc --- ---------------------------------------
CSJAc --- --- barrel search starting from TPC ---
CSJAc --- ---------------------------------------
CSJA
CSJA      call vzero(method,detmax)
CSJAc --- set extrapolation methods
CSJA      method(i_VTX) =imethod
CSJA      method(i_FTD) =0
CSJA      method(i_ITC) =0
CSJA      method(i_SIT1)=0
CSJA      method(i_SIT2)=0
CSJA      method(i_FCH)=imethod
CSJA
CSJAc --- init timing
CSJA      CALL TIMEX(TSTART)
CSJA
CSJAc --- perform extrapolation
CSJA      if(fsfdbg.gt.0)then
CSJA        write(6,*) '****************************************'
CSJA        print*,'FSFSTR: Barrel Track Search from TPC'
CSJA        write(6,*) '****************************************'
CSJA      endif
CSJA      method(i_TPC) = none      ! this is the pivot
CSJA      call fsfana3 (i_TPC,method,
CSJACC     1     data_labl+FLAG_GIVESINGLTE+FLAG_NOPLINV+FLAG_USEUSED+
CSJA     1     data_labl+FLAG_GIVESINGLTE+FLAG_NOPLINV+
CSJA     2     FLAG_CUT100)
CSJA
CSJAc --- read timing
CSJA      CALL TIMEX(TEND)
CSJA      TIME(i_TPC) = TIME(i_TPC)+(TEND-TSTART)
CSJA      TSTART      = TEND
CSJA
CSJAc --- remove duplicate and sub-TS
CSJA      call fsfnsub(FLAG_DBL+FLAG_SUB)
CSJA
CSJAc --- list of TS as found up to now
CSJA      if (fsfdbg.gt.1) call fsftsdump
CSJA
CSJAc --- -------------------------------------------------------
CSJAc --- --- now cleanup and iterate searches in SIT,TPC,FCH ---
CSJAc --- -------------------------------------------------------
CSJA
CSJA      call vzero(method,detmax)
CSJAc --- set extrapolation methods
CSJA      method(i_VTX) =0
CSJA      method(i_FTD) =0
CSJA      method(i_ITC) =0
CSJA      method(i_SIT1)=imethod
CSJA      method(i_SIT2)=imethod
CSJA      method(i_TPC) =imethod
CSJA      method(i_FCH) =imethod
CSJA
CSJAc --- init timing
CSJA      CALL TIMEX(TSTART)
CSJA
CSJAc --- now start new search to find all TEs not yet found before
CSJA      data_labl = FLAG_DATA
CSJA      if(fsfdbg.gt.0)then
CSJA        write(6,*) '****************************************'
CSJA        print*,'FSFSTR: STRING extrapolation to SIT/FCH '
CSJA        write(6,*) '****************************************'
CSJA      endif
CSJA
CSJAc --- now create new pseudo TEs from all TS we have found.
CSJA      call fsfrest(i_SCAN1)
CSJA
CSJA      method(i_SCAN1) = none  ! this is the pivot
CSJA      call fsfana3 (i_SCAN1,method,
CSJA     1                data_labl+FLAG_SCANDET+FLAG_NOPLINV+
CSJA     2                FLAG_CUT100)
CSJA
CSJAc --- read timing
CSJA      CALL TIMEX(TEND)
CSJA      TIME(i_SCAN1) = TIME(i_SCAN1)+(TEND-TSTART)
CSJA      TSTART        = TEND
CSJA
CSJAc --- remove double and sub-strings.
CSJA      if(fsfdbg.gt.0)then
CSJA        print*,'FSFSTR: remove substrings'
CSJA      endif
CSJA      call fsfnsub(FLAG_DBL+FLAG_SUB)
CSJA
CSJAc --- ----------------------------------------------------------
CSJAc --- --- now cleanup and iterate searches in FTD, outer VTX ---
CSJAc --- ----------------------------------------------------------
CSJA
CSJA      call vzero(method,detmax)
CSJAc --- set extrapolation methods
CSJA      method(i_VTX) =0
CSJA      method(i_VTX1)=0
CSJA      method(i_VTX2)=0
CSJA      method(i_VTX3)=imethod
CSJA      method(i_VTX4)=imethod
CSJA      method(i_VTX5)=imethod
CSJA      method(i_FTD) =0
CSJA      method(i_FTD1)=imethod
CSJA      method(i_FTD2)=imethod
CSJA      method(i_FTD3)=imethod
CSJA      method(i_FTD4)=imethod
CSJA      method(i_FTD5)=imethod
CSJA      method(i_FTD6)=imethod
CSJA      method(i_FTD7)=imethod
CSJA      method(i_ITC) =0
CSJA      method(i_SIT1)=0
CSJA      method(i_SIT2)=0
CSJA      method(i_TPC) =0
CSJA      method(i_FCH)=0
CSJA
CSJAc --- init timing
CSJA      CALL TIMEX(TSTART)
CSJA
CSJAc --- now start new search to find all TEs not yet found before
CSJA      data_labl = FLAG_DATA
CSJA      if(fsfdbg.gt.0)then
CSJA        write(6,*) '********************************************'
CSJA        print*,'FSFSTR: STRING extrapolation to FTD/outer VTX'
CSJA        write(6,*) '********************************************'
CSJA      endif
CSJA
CSJAc --- now create new pseudo TEs from all TS we have found.
CSJA      call fsfrest(i_SCAN2)
CSJA
CSJA      method(i_SCAN2) = none  ! this is the pivot
CSJA      call fsfana3 (i_SCAN2,method,
CSJA     1                data_labl+FLAG_SCANDET+FLAG_NOPLINV+
CSJA     2                FLAG_CUT100)
CSJA      method(i_SCAN2) = 0
CSJA
CSJAc --- read timing
CSJA      CALL TIMEX(TEND)
CSJA      TIME(i_SCAN2) = TIME(i_SCAN2)+(TEND-TSTART)
CSJA      TSTART        = TEND
CSJA
CSJAc --- remove double and sub-strings.
CSJA      if(fsfdbg.gt.0)then
CSJA        print*,'FSFSTR: remove substrings'
CSJA      endif
CSJA      call fsfnsub(FLAG_DBL+FLAG_SUB)
CSJA
CSJAc --- -----------------------------------------------------
CSJAc --- --- now cleanup and iterate searches in inner VTX ---
CSJAc --- -----------------------------------------------------
CSJA
CSJA      call vzero(method,detmax)
CSJAc --- set extrapolation methods
CSJA      method(i_VTX2)=imethod
CSJA
CSJAc --- init timing
CSJA      CALL TIMEX(TSTART)
CSJA
CSJAc --- now start new search to find all TEs not yet found before
CSJA      data_labl = FLAG_DATA
CSJA      if(fsfdbg.gt.0)then
CSJA        write(6,*) '*****************************************'
CSJA        print*,'FSFSTR: STRING extrapolation to inner VTX'
CSJA        write(6,*) '*****************************************'
CSJA      endif
CSJA
CSJAc --- now create new pseudo TEs from all TS we have found.
CSJA      call fsfrest(i_SCAN3)
CSJA
CSJA      method(i_SCAN3) = none  ! this is the pivot
CSJA      call fsfana3 (i_SCAN3,method,
CSJA     1                data_labl+FLAG_SCANDET+FLAG_NOPLINV+
CSJA     2                FLAG_CUT100)
CSJA      method(i_SCAN3) = 0
CSJA
CSJAc --- read timing
CSJA      CALL TIMEX(TEND)
CSJA      TIME(i_SCAN3) = TIME(i_SCAN3)+(TEND-TSTART)
CSJA      TSTART        = TEND
CSJA
CSJAc --- remove double and sub-strings.
CSJA      if(fsfdbg.gt.0)then
CSJA        print*,'FSFSTR: remove substrings'
CSJA      endif
CSJA      call fsfnsub(FLAG_DBL+FLAG_SUB)
CSJA
CSJAc --- -----------------------------------------------------
CSJAc --- --- now cleanup and iterate searches in inner VTX ---
CSJAc --- -----------------------------------------------------
CSJA
CSJA      call vzero(method,detmax)
CSJAc --- set extrapolation methods
CSJA      method(i_VTX1)=imethod
CSJA
CSJAc --- init timing
CSJA      CALL TIMEX(TSTART)
CSJA
CSJAc --- now start new search to find all TEs not yet found before
CSJA      data_labl = FLAG_DATA
CSJA      if(fsfdbg.gt.0)then
CSJA        write(6,*) '*****************************************'
CSJA        print*,'FSFSTR: STRING extrapolation to inner VTX'
CSJA        write(6,*) '*****************************************'
CSJA      endif
CSJA
CSJAc --- now create new pseudo TEs from all TS we have found.
CSJA      call fsfrest(i_SCAN4)
CSJA
CSJA      method(i_SCAN4) = none  ! this is the pivot
CSJA      call fsfana3 (i_SCAN4,method,
CSJA     1                data_labl+FLAG_SCANDET+FLAG_NOPLINV+
CSJA     2                FLAG_CUT100)
CSJA      method(i_SCAN4) = 0
CSJA
CSJAc --- read timing
CSJA      CALL TIMEX(TEND)
CSJA      TIME(i_SCAN4) = TIME(i_SCAN4)+(TEND-TSTART)
CSJA      TSTART        = TEND
CSJA
CSJAc --- remove double and sub-strings.
CSJA      if(fsfdbg.gt.0)then
CSJA        print*,'FSFSTR: remove substrings'
CSJA      endif
CSJA      call fsfnsub(FLAG_DBL+FLAG_SUB)
CSJA
CSJAc --- ----------------------------------------
CSJAc --- --- refit without error scalings     ---
CSJAc --- ----------------------------------------
CSJA
CSJAc --- final fit without any error scaling
CSJA      if(fsfdbg.gt.0)then
CSJA        write(6,*) '****************************************'
CSJA        print*,'FSFSTR: refit without error scaling to clean up'
CSJA        write(6,*) '****************************************'
CSJA      endif
CSJA      call fsfrefit
CSJA
CSJAc --- now do the real test for funny TS
CSJA      if(fsfdbg.gt.0)then
CSJA        write(6,*) '****************************************'
CSJA        print*,'FSFSTR: remove double and test funny TS in list'
CSJA        write(6,*) '****************************************'
CSJA      endif
CSJA      call fsfnsub(FLAG_DBL+FLAG_SUB)
CSJA
CSJAc --- ----------------------------------------
CSJAc --- --- now remove funny TS combinations ---
CSJAc --- ----------------------------------------
CSJA
CSJA      do i = 1,tsindexused
CSJA        if(tsindex(i).ne.0)then
CSJA          call fsftts(tsindex(i),badcomb,badsingle)
CSJA          if(tsindex(i).eq.0)then
CSJA            if(fsfdbg.gt.0)then
CSJA              print*,'FSFSTR: funny TS not given !'
CSJA            endif
CSJA          endif
CSJA        endif
CSJA      enddo
CSJAc --- remove double and sub-strings.
CSJA      if(fsfdbg.gt.0)then
CSJA        print*,'FSFSTR: remove substrings'
CSJA      endif
CSJA      call fsfnsub(FLAG_DBL+FLAG_SUB)
CSJA
CSJAc --- ----------------------------------------
CSJAc --- --- save final TS candidates         ---
CSJAc --- ----------------------------------------
CSJA
CSJAc --- Save all the TS and make TK of what we have made.
CSJA      call fsfsts(FLAG_GIVESINGLTE)
CSJA
CSJA      if(fsfdbg.gt.0) then
CSJA        write(6,*) '****************************************'
CSJA        print*,'FSFSTR: event done !'
CSJA        write(6,*) '****************************************'
CSJA      endif
      return
      END


