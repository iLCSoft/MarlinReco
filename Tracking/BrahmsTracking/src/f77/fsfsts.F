*  
* $Id: fsfsts.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfsts(flags)
*******************************************************************************
*     This routine shall save all TSs that are in the
*     tsindex array with a non-zero entry.
*
*     INPUT: flags - FLAG_LIFTTK to lift TK structure as well (IDVD mode)
*                    FLAG_GIVESINGLTE
*
*     FEATURES:
*     - It does handle exclusions by calling fsexcl for each TS.
*
*******************************************************************************
      implicit none
#include "include/brcntl.inc"
#include "include/fsfpi.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfclun.inc"
#include "include/fsfflags.inc"
#include "include/fsfsim.inc"
      integer i,j,n,idts,idtsr,ierr,itmp,pd,k
      integer labl,jbit
      integer idtkr,idtk,idatts(tssize)
      integer  ilabel
c      external EVTESE,EVFTES
      integer flags
      logical LIFTTK,GIVSNG

      LIFTTK = jbit(flags,FBIT_LIFTTK).eq.1
      GIVSNG = jbit(flags,FBIT_GIVESINGLTE).eq.1
      if(LIFTTK.and.fsfdbg.gt.0) then
        print*,'FSFSTS: write also TK structure'
      endif
      if(GIVSNG.and.fsfdbg.gt.0) then
        print*,'FSFSTS: return also TS with only a single TE'
      endif

      do i=1,tsindexused
c --- check that track parameters make any sense
        if(tsindex(i).ne.0)then
          if (space(tsindex(i)+20).lt.0. .or.
     1        space(tsindex(i)+20).gt.pi .or.
     2        space(tsindex(i)+21).lt.0. .or.
     3        space(tsindex(i)+21).gt.2*pi) then
            print *,'ERROR in FSFSTS: fit gave theta/phi off'
            print *,'TS dump: ',(ispace(tsindex(i)+k),k=1,8),
     1                           (space(tsindex(i)+k),k=9,23)
            print *,'drop it !'
            tsindex(i) = 0
          endif
        endif

c --- tsindex(i).eq.0 may be used to remove TSs.
        if(tsindex(i).ne.0)then
c --- minimum of two TEs in string !!!!
          if ((ispace(tsindex(i)+5).eq.1).and..not.GIVSNG)then
            if(fsfdbg.gt.0)then
              print*,'FSFSTS: single TE not given as TS !'
            endif
          else
c --- Save the TS Data.
            n=ispace(tsindex(i)+5)
            call tkstor('TS',idts,ispace(tsindex(i)+1),
     >                       n,ispace(tsindex(i)+tssize+1),0,ierr)
            idtsr=idts
            if(ierr.ne.0) then
              print*,'ERROR in FSFSTS: TKSTOR failed: ierr =',ierr
              CALL TKVETO(1)
              idts=0
            endif
            if (idts.ne.0) then
              if(fsfdbg.gt.1)then
                write(IPRNT,1000) idtsr
 1000 FORMAT(1X,'FSFSTS: written TSR has id:',I8)
              endif

c --- lift TK if needed (IDVD mode)
              if (LIFTTK) then
                call ucopy (space(tsindex(i)+1),idatts,tssize)
c --- set TKR(3) correct, it is different for TSR and TKR
                if (jbit(idatts(3),2).eq.1) then
                  idatts(3) = 1
                else
                  idatts(3) = 0
                endif
c --- save TKR
                write(*,*) "about to call tkstor"
                call tkstor('TK',idtk,idatts,0,0,0,ierr)
                idtkr=idtk
c                call TSAVF ('TK',1,idtsr,tssize,idatts,idtk,idtkr)
c                if (idtkr.eq.0) then
c                  call tquery(ierr)
                  if (ierr.ne.0) then
                    print*,'ERROR in FSFSTS: TKSTOR failed: ierr =',ierr
                    CALL TKVETO(1)
                  end if
c                else
                  if(fsfdbg.gt.1)then
                     write(IPRNT,1001) idtkr
 1001 FORMAT(1X,'FSFSTS: written TKR has id:',I8)
c                  endif
c --- give TK a label
c                  call extkla(idtkr)
                endif
              endif

c
c --- **** Some statistics about TSs. ****
c
              pd = ispace(tsindex(i)+1)-3000
              ntspd(pd)    = ntspd(pd)+1
              do j=1,detmax
                if(1.eq.jbit(ispace(tsindex(i)+2),detid(j)+1))then
                  detsum(pd,j) = detsum(pd,j)+1
                  detsum(0,j)  = detsum(0,j)+1
                endif
              enddo
              if(fsfdbg.gt.0)then
c --- receive label
                call tkglab('TS',idtsr,ilabel,ierr)
                if (ierr.ne.0) then
                  write(6,*) 'FSFSTS: TKGLAB failed.'
                  stop
                end if
                print*,'FSFSTS: Saved TS of module',
     1                 ispace(tsindex(i)+1),' label',ilabel
              endif
c
c --- ** Statistics about track charge
c
              if(space(tsindex(i)+22).eq.0)then
                nsign(3)=nsign(3)+1 !Unknown charge
              else if(space(tsindex(i)+22).gt.0)then
                nsign(2)=nsign(2)+1 !Negative charge
              else if(space(tsindex(i)+22).lt.0)then
                nsign(1)=nsign(1)+1 !Positive charge
              endif
c
c --- ** Statistics about track direction
c
              if(space(tsindex(i)+20).lt.Pi/2)then
                ndirect(1)=ndirect(1)+1 !Forward
              else if(space(tsindex(i)+20).gt.Pi/2)then
                ndirect(2)=ndirect(2)+1 !Backward
              endif
#if defined(HBOOK)
              call hcdir('//PAWC/TRACKING/DELSEARCH',' ')
              if (fsfhist.gt.0) then
                if (delsim) then
                  call tkglab('TS',idtsr,labl,ierr)
                  if(ierr.ne.0)then
                    print*,'ERROR in FSFSTS: TKGLAB failed: ierr =',
     1                     ierr
                  endif
                endif
c --- ** Statistics about the theta distribution
                IF (labl.ne.0.) THEN
                  call hf1(310110,space(tsindex(i)+20)*180/PI,1.)
                endif
              endif
#endif
            endif
          endif
        endif
      enddo
      END

