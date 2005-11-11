*  
* $Id: fsfplinv.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfplinv(idxv,pd,flags,accepted)
*******************************************************************************
*
*      Polar inversion of string candidate, cuts and then refit to build
*      TS out of good candidates.
*
*      INPUT:  idxv     - vector of TEs
*              pd       - primary detector
*              flags    - flags to steer methods
*      OUTPUT: accepted - good TS ?
*                                                                  Daniel Wicke
*******************************************************************************
      implicit none
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfflags.inc"
#include "include/fsfreject.inc"
#include "include/fsfpi.inc"
      integer i,j,pd,flags,idet,ntot
      integer idxv(detmax),itsnew
      integer jbit
      real    fsfpdiff,sum,wgtsum,the
      logical accepted
      logical uselabl,useused,givesnglte,averagetheta,
     1        noplinv

      real    a2,xp(detmax),yp(detmax),chi2,hin(3,detmax),hout(3)
      real    a(2)
      real    helixpar(6),helixerr(6)
      real    deltaphi,pinv
      integer ierr

      uselabl     =jbit(flags,FBIT_LABL).eq.1
      useused     =jbit(flags,FBIT_USEUSED).eq.1
      givesnglte  =jbit(flags,FBIT_GIVESINGLTE).eq.1
      averagetheta=jbit(flags,FBIT_AVERAGETHETA).eq.1
      noplinv     =jbit(flags,FBIT_NOPLINV).eq.1

      accepted = .FALSE.
      ntot    = 0

      do i=1,detmax
        if(idxv(i).ne.0)then
          ntot=ntot+1
c --- do we need polar inversion ?
          if (.not.noplinv) then
            if(jbit(datdescr(idxv(i),i),1).eq.0)then
               a2=1.0
               call fsplxy(a2,idxv(i),i,xp,yp,hin,ntot)
            else
               a2=1.0
               call fsplrp(a2,idxv(i),i,xp,yp,hin,ntot)
            endif
          endif
        endif
      enddo

c -----------------------------------------------------------------------------
c --- use primary detector to define parameters
c -----------------------------------------------------------------------------
      if((ntot.eq.1.and.givesnglte).or.noplinv)then
c ---  *Check if all coordinates we need are measured.
        if(1.eq.jbit(datdescr(idxv(pd),pd),3).and.
     1     1.eq.jbit(datdescr(idxv(pd),pd),4).and.
     2     1.eq.jbit(datdescr(idxv(pd),pd),5).and.
     3     1.eq.jbit(datdescr(idxv(pd),pd),8).and.
     4     1.eq.jbit(datdescr(idxv(pd),pd),9).and.
     5    (1.eq.jbit(datdescr(idxv(pd),pd),10).or.
     6     1.eq.jbit(datdescr(idxv(pd),pd),11) ) )then
c --- Give primary detector as TS.
          if(fsfdbg.gt.1) then
            print*,'FSFPLINV: no polar inversion !'
            if(ntot.eq.1) print*,'FSFPLINV: Created TS from single TE.'
          endif
          accepted=.TRUE.
c --- Crude estimations FK3TRK will fix this for us anyway.
          call fsfcpte(idxv(pd),pd,helixpar)
        endif
c -----------------------------------------------------------------------------
c --- now polar inversion
c -----------------------------------------------------------------------------
      else if(ntot.ge.2)then
        call fsflft(xp,yp,hin,a,hout,1,ntot,chi2)
c --- fallback on primary detector in case of error
        if(chi2.gt.1000.)then
          if(1.eq.jbit(datdescr(idxv(pd),pd),3).and.
     1       1.eq.jbit(datdescr(idxv(pd),pd),4).and.
     2       1.eq.jbit(datdescr(idxv(pd),pd),5).and.
     3       1.eq.jbit(datdescr(idxv(pd),pd),8).and.
     4       1.eq.jbit(datdescr(idxv(pd),pd),9).and.
     5      (1.eq.jbit(datdescr(idxv(pd),pd),10).or.
     6       1.eq.jbit(datdescr(idxv(pd),pd),11) ) )then
c --- Give primary detector as TS.
            if(fsfdbg.gt.1) print*,
     1        'FSFPLINV: Linefit failed, created TS using primary TE.'
            accepted=.TRUE.
c --- Crude estimations FK3TRK will fix this for us anyway.
            call fsfcpte(idxv(pd),pd,helixpar)
          else
            if(fsfdbg.gt.1)then
              print*,'FSFPLINV: Linefit failed, reject combination.'
            endif
          endif
        else
c -----------------------------------------------------------------------------
c --- use polar inversion to define parameters
c -----------------------------------------------------------------------------
c --- Helixpar from a and hout, starts at origine per definition
          helixpar(2)=0  !x
          helixpar(3)=0  !y
          helixpar(1)=0  !z

          if(averagetheta)then
c --- get theta from all TEs using the measured Z and origine
            sum    = 0.
            wgtsum = 0.
            do idet=1,detmax
c ---
              if (idxv(idet).ne.0) then
c --- do we have a Z from this TE ?
                if (jbit(datdescr(idxv(idet),idet),7).eq.1) then ! z measured ?
c --- we do have a z measurement and we take our R from local arrays
                  the = pi/2.
                  if (z(idxv(idet),idet).ne.0.) then
c --- take theta from dr/dz !!!
                    the  = atan(r(idxv(idet),idet)/z(idxv(idet),idet))
                    if (the.lt.0.) the = the+pi
                  endif
                  sum    = sum    + the
                  wgtsum = wgtsum + 1.0
                endif
              endif
            enddo
            if (wgtsum.eq.0.) then
              if (fsfdbg.gt.0)then
                print*,'FSFPLINV: theta average is zero - use primary'
              endif
              helixpar(4)=theta(idxv(pd),pd)
            else
              helixpar(4)=sum/wgtsum
            endif
          else
c --- use theta from primary detector
            helixpar(4)=theta(idxv(pd),pd)
          endif
#if defined(HBOOK)
          call hcdir('//PAWC/TRACKING/DELSEARCH',' ')
          if (fsfhist.gt.0) then
            if(simtheta(idxv(pd),pd).ne.0)then
              call hfill(310060,theta(idxv(pd),pd)
     1                          -simtheta(idxv(pd),pd),0,1.0)
              call hfill(310061,helixpar(4)
     1                          -simtheta(idxv(pd),pd),0,1.0)
            endif
          endif
#endif

c --- Phi shall be positiv
          helixpar(5)=atan(a(1)) !Phi
          if(helixpar(5).lt.0)  helixpar(5)=helixpar(5)+Pi
c --- Solve ambiguity phi <--> phi+pi
          if (R(idxv(pd),pd).gt.5.) then
c --- compare to Rphi/R of primary detector
            if(abs(fsfpdiff(helixpar(5),
     1             Rphi(idxv(pd),pd)/R(idxv(pd),pd))) .gt. pi/2)then
              helixpar(5)=helixpar(5)+Pi
            endif
          else
c --- compare to phi of primary detector
            if(abs(fsfpdiff(helixpar(5),
     1             phi(idxv(pd),pd))) .gt. pi/2)then
              helixpar(5)=helixpar(5)+Pi
            endif
          endif
c --- use primary detector to define sign of curvature
          deltaphi=fsfpdiff(phi(idxv(pd),pd),helixpar(5))

c --- Take the sign for R_T from deltaphi.
          helixpar(6)=sign(2*a(2)/sqrt(1+a(1)**2),deltaphi) !1/R_T
          if(helixpar(6).eq.0)then
            if(fsfdbg.gt.0) then
               print*,'WARNING in FSFPLINV: 1/R(fit)==0',
     1                ' using 1/R(pd) instead'
            endif
            helixpar(6)=overrt(idxv(pd),pd)
          endif

c --- Take the sign for R_T from primary detector.
c     deltaphi from pd is better (DW 15.1.96)
CC            helixpar(6)=sign(2*a(2)/sqrt(1+a(1)**2),overrt(idxv(pd),pd)) !1/R_T

#if defined(HBOOK)
          call hcdir('//PAWC/TRACKING/DELSEARCH',' ')
          if (fsfhist.gt.0) then
            if(simoverrt(idxv(pd),pd).ne.0)then
              call hfill(310070,overrt(idxv(pd),pd)
     1                          -simoverrt(idxv(pd),pd),0,1.0)
              call hfill(310071,helixpar(6)
     1                          -simoverrt(idxv(pd),pd),0,1.0)
            endif
          endif
#endif
c -----------------------------------------------------------------------------
c --- Now cut on phi:
c     This cut should be stronger, than the previous cuts as the fitted
c     phi should be much better.
c -----------------------------------------------------------------------------
c --- fsfplcut cannot cope with Theta==0!
          if(helixpar(4).eq.0)then
             print*,'WARNING in FSFPLINV: Theta==0!!!'
          else
            call fsfplcut(helixpar,idxv,accepted,pd)
            if(accepted)then
              if(fsfdbg.gt.1) then
                print *,'FSFPLINV: accepted after fsfplcut !'
              endif
            else
              if(fsfdbg.gt.0.and.labelsok)then
                print *,'WARNING in FSFPLINV:',
     1                  ' good TS not accepted after fsfplcut !'
              else if(fsfdbg.gt.1) then
                print *,'FSFPLINV: not accepted after fsfplcut !'
              endif
            endif
          endif
        endif
      endif
c -----------------------------------------------------------------------------
c --- Write TS ?
c -----------------------------------------------------------------------------
      if(accepted)then
        if(fsfdbg.gt.0)then
          print*,'FSFPLINV: fit TS for this combination:'
            do idet=1,detmax
              if (idxv(idet).ne.0)then
              print'(A,i4,A,A5,A,I8,A,i8)',
     1        '      TE # ',idxv(idet),
     2        ' from ',detname(idet),
     3        ' id #',teid(idxv(idet),idet),
     4        ' label ',telabl(idxv(idet),idet)
            endif
          enddo
        endif

        call fsfwts(idxv,helixpar,helixerr,pd,flags,ierr)

        if (ierr.ne.0) then
          itsnew = 0
        else
          itsnew = tsindex(tsindexused)
          pinv   = space(itsnew+22)
        endif
c -----------------------------------------------------------------------------
c --- now fit TS
c -----------------------------------------------------------------------------
        if (itsnew.ne.0)then
          call fsfits(itsnew,flags,pd)
          if(fsfdbg.gt.1)then
            if(itsnew.eq.0)then
              print*,'FSFPLINV: first trackfit failed !'
            else if(space(itsnew+22)*pinv.lt.0)then
              print*,'FSFPLINV: first fit changed sign of 1/p.'
            endif
          endif
        endif

        if (itsnew.ne.0)then
          if(ispace(itsnew+5).GT.1)then
c --- now fit with outlayer
            call fsfits(itsnew,flags+FLAG_OUTLAYER,pd)
            if(fsfdbg.gt.1.and.itsnew.eq.0)then
              print*,'FSFPLINV: trackfit with outlayer failed !'
            endif
c --- test if funny TS
            if(itsnew.ne.0)then
              call fsftts(itsnew,badcomb,0) ! single VFT for VD
              if(fsfdbg.gt.1.and.itsnew.eq.0)then
                print*,'FSFPLINV: TS after outlayer is funny !'
              endif
            endif
          endif
        endif
c -----------------------------------------------------------------------------
c --- now check that we are not left with only one "TE" from IDXV array
c -----------------------------------------------------------------------------
        if (itsnew.ne.0.and.ntot.ne.1) then
          if (ispace(itsnew+5).eq.1) then
            itsnew = 0
            if(fsfdbg.gt.0) then
              print *,'FSFPLINV: all but one TE outlayer, reject TS'
            endif
          else
c --- in case of multi-TEs in a pseudo TE I need at least one other TE
            do i = 1,detmax
              if (idxv(i).ne.0) then
c --- try to find this TE identifier in the TE list of the TS
                do j = 1,ispace(itsnew+5)
                  if (ispace(itsnew+tssize+j).eq.teid(idxv(i),i)) then
c --- fine, this TS contains non trivial information
                    GOTO 1000
                  endif
                enddo
              endif
            enddo
            itsnew = 0
            if(fsfdbg.gt.0) then
              print *,'FSFPLINV: TS after outlayer is multi-TE (PD),',
     1                ' reject TS'
            endif
          endif
        endif
 1000   CONTINUE

c --- TS rejected ?
        if (itsnew.eq.0)then
          if (tsindexused.gt.0) then
            spaceused            = tsindex(tsindexused)
            tsindex(tsindexused) = 0
            tsindexused          = tsindexused-1
          else
            write(6,*) 'FSFPLINV: warning - check if space is ok'
          endif
          accepted=.FALSE.
        endif
      endif

      if(fsfdbg.gt.0)then
        if (accepted) then
          print*,'FSFPLINV: TS is accepted !'
        else
          print*,'FSFPLINV: TS is not accepted !'
        endif
      endif
      end

