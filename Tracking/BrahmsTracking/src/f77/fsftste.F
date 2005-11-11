*  
* $Id: fsftste.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      SUBROUTINE FSFTSTE(ts_index,te_index,ndat,ndat2,ishower)
********************************************************************
* Routine to transform a TS structure into a dummy TE structure.
*
* INPUT:  ts_index - index in space array to find the TS
*         te_index - index in space array to put TE, can be ts_index !
* OUTPUT: ndat2    - length of TE
*                                                      Markus Elsing
********************************************************************
      IMPLICIT NONE

#include "include/itcdims.inc"
#include "include/tpcdims.inc"
#include "include/fsfparams.inc"
#include "include/fsfflags.inc"
#include "include/brcntl.inc"
#include "include/fsfarrays.inc"
#include "include/fsfpi.inc"

      EXTERNAL numbit
      INTEGER  ts_index,te_index,ndat,ndat2,jbit,
     1         IERRSF,ier,numbit,ishower
      REAL     rstart,rend,r_surf
      REAL     TEXTRA(30)

      if(ndat.eq.0)then
        ndat2 = 0
        return
      else if (jbit(ispace(ts_index+3),1).eq.0)then
        print *,'WARNING in FSFTSTE: we need track parameters here !'
        ndat2 = 0
        return
      else

c --- find surface to extrapolate
        r_surf = 0. ! do not extrapolate !
c --- do we have more than one detector ?
        if (numbit(ispace(ts_index+2)).gt.1) then
c --- get middle of TS in R
          rstart = space(ts_index+11)**2+space(ts_index+12)**2
          rend   = space(ts_index+14)**2+space(ts_index+15)**2
          if (rstart.eq.0.or.rend.eq.0.)then
            if (jbit(ispace(ts_index+2),detid(i_TPC)+1).eq.1) then
              r_surf = (rtpcinn+rtpcout)/2.0   ! give TPC
            else if (jbit(ispace(ts_index+2),detid(i_ITC)+1).eq.1) then
              r_surf = (ritcinn+ritcout)/2.0   ! give ITC
            else
              if(fsfdbg.gt.0)then
                 print *,'WARNING in FSFTSTE: no R to extrapolate'
              endif
            endif
          else
            r_surf = (sqrt(rstart)+sqrt(rend))/2. ! give middle of R
          endif
        endif

c --- extrapolate TS
        if (r_surf.gt.0.) then
          call fsftsexx(ts_index,1,'CYLI',r_surf,-500.,500.,textra,
     >                  ierrsf,ier)
          if (.not.(ier+ierrsf.eq.0)) then
            if(fsfdbg.gt.0)then
              print *,'WARNING in FSFTSTE: extrapolation failed',
     1                ier,ierrsf
            endif
            r_surf = 0.
          endif
        endif

c --- get shower bits to the outside world, write into bits 17-19
        call sbyt(ispace(ts_index+4),ishower,17,3)

c --- now modify data in order to get TE out of it
c --- 1: module id -> keep TS info
        ispace(te_index+1) = ispace(ts_index+1)
c --- 2: submodule -> keep TS detector bitted word
        ispace(te_index+2) = ispace(ts_index+2)
c --- 5: pointer to mass id.= 0
        ispace(te_index+5) = 0
c --- 6: charge    -> keep TS info
        ispace(te_index+6) = ispace(ts_index+6)
c --- 7: NDF       -> from TS info
        ispace(te_index+7) = ispace(ts_index+8)
c --- 8: chi2      -> from TS info
        ispace(te_index+8) = ispace(ts_index+9)
c --- 9: length    -> from TS info
        ispace(te_index+9) = ispace(ts_index+10)
c --- 10-15: trackparameter -> from TS info
        CALL ucopy (space(ts_index+17),space(te_index+10),6)
c --- 16: beta ???? put 0
        ispace(te_index+16) = 0
c
c --- copy extrapolation result ?
c
        if (r_surf.gt.0.)then
c --- 4: measurement code -> exx is cylinder !
          ispace(te_index+4) = 1021 ! Rphi and 1/p
c --- 17-31: cov matrix     -> from TSX info
          call ucopy(textra(7) ,space(te_index+10),6)
          call ucopy(textra(13),space(te_index+17),15)
        else
c --- 4: measurement code -> decode TS, if plane or cylinder !
          if (jbit(ispace(ts_index+3),2).eq.1) then
           ispace(te_index+4) = 1021 ! Rphi and 1/p
          else
           ispace(te_index+4) = 1020 ! xyz  and 1/p
          endif
c --- 17-31: cov matrix     -> from TS info
          if (jbit(ispace(ts_index+3),3).eq.0.OR.ndat.lt.tssize) then
            call vzero (ispace(te_index+17),15)! do not use diagonal estimates
                                               ! this is the normal situation
                                               ! for TKSB since it gives no
                                               ! full TS info
          else
            call ucopy (space(ts_index+23),space(te_index+17),15) ! full matrix
          endif
        endif
c --- 3: extrapolation bank = 0
        ispace(te_index+3) = 0

        ndat2 = 31

        if(fsfdbg.gt.1)then
          print*,'FSFTSTE: TE-dump:'
          call fsfprintte(ispace(te_index+1),space(te_index+1),ndat2)
        endif

      endif
      end


