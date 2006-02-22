*  
* $Id: fsfits.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfits(its,flags,pd)
*******************************************************************************
*
*      Routine to fit TS and run beam constraint and outlayer...
*
*      INPUT:  its   - TS index in array
*              flags - FLAG_BEAMSPOT triggers beam spot constraint fit
*                      FLAG_OUTLAYER triggers outlayer removal
*                      FLAG_NOSCALING triggers no error scalings in fit
*                      FLAG_EXTRAPOL triggers full extrapolations for fit
*              pd    - primary detector which should not be removed in outlayer
*
*                                                                 Markus Elsing
*******************************************************************************
      implicit none
#include "include/brcntl.inc"
#include "include/fkdebug.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfflags.inc"
#include "include/fsfxerrs.inc"
#include "include/fsfpi.inc"
#include "include/fsfbits.inc"
#include "include/fsfsim.inc"
#include "include/fktkx.inc"

      integer anint
      integer its,flags
      integer n,jbit,k
      integer IDENT(NTEMX)
      integer IDATTE(MXTE,NTEMX)
      real    DATTE(MXTE,NTEMX)
      equivalence (idatte,datte)
      integer itk(MXTK)
      real     tk(MXTK),cov(15),r_ts
      equivalence (itk(1),tk(1))
      integer ier1,ier2,nout,idout(NTEMX),iflte(NTEMX),i,j,ierr,ikeep
      integer ier,ierx(NEXSMX),pd,ioutr,itefsf,index,num,ilabel,labelte
      logical usebeamspot,useoutlayer,usescaling,useextrapol
c
      real textra(30)
      integer   ierrsf

      usebeamspot=jbit(flags,FBIT_BEAMSPOT).eq.1
      useoutlayer=jbit(flags,FBIT_OUTLAYER).eq.1
      usescaling =jbit(flags,FBIT_NOSCALING).eq.0  !!!
      useextrapol=jbit(flags,FBIT_EXTRAPOL).eq.1

c --------------------------------------------------------------------------
c --- check that input track parameters make any sense
c --------------------------------------------------------------------------
      if (space(its+20).lt.0. .or. space(its+20).gt.pi .or.
     1    space(its+21).lt.0. .or. space(its+21).gt.2*pi) then
        if (fsfdbg.gt.0) then
          print *,'WARNING in FSFITS: entered fit with theta/phi off'
          print *,'TS dump:',(ispace(its+k),k=1,8),(space(its+k),k=9,23)
          print *,'drop it !'
        endif
        its = 0
        RETURN
      endif

c --------------------------------------------------------------------------
c --- load TE data from tracking banks direct
c --------------------------------------------------------------------------
      n = ispace(its+5)
      call vzero(IDENT,ntemx)
      call ucopy(space(its+tssize+1),IDENT,n)
c
      do i=1,n
        CALL TKREAD('TE',0,ident(i),datte(1,i),idatte(1,i),IERR)
        
        IF (IERR.NE.0) THEN
          WRITE(6,*) 'FSFITS: TKREAD FAILED.'
          STOP
        END IF
        CALL FSFCONVRT(idatte(1,i))
      enddo

c --------------------------------------------------------------------------
c --- scale error for trackfit
c --------------------------------------------------------------------------
      if(usescaling.and.(fiterrscale.ne.1..or.tpcerrscale.ne.1.))then
        do itefsf = 1,n
          index=0
          do i=1,6
            if(jbit(idatte(4,itefsf),5+i).eq.1) index=index+1
          enddo
          num=16+(index*(index+1))/2
          if (idatte(1,itefsf).eq.i_TPC) then
            do i=17,num
              datte(i,itefsf) = datte(i,itefsf) * tpcerrscale**2
            enddo
          else
            do i=17,num
              datte(i,itefsf) = datte(i,itefsf) * fiterrscale**2
            enddo
          endif
        enddo
      endif

c --------------------------------------------------------------------------
c --- Fake TE for Beamspot/Primary vertex:
c --------------------------------------------------------------------------
      if(usebeamspot)then
        if(.not.n.lt.maxtev)then
          print*,'ERROR in FSFITS: To many TEs in TS. ',
     1           'Could not fit with beamspot.'
          print*,'                 Increase parameter ',
     1           ' maxtev to cure this problem.'
          CALL TKVETO(1)
        else
c --- now create dummy TE based on beam spot, n=n+1 !!!
          if(fsfdbg.gt.0) then
            print*,'FSFITS: Fitting with beamspot constraint.'
          endif
          call FSFBEAMTE(datte,n)
c --- put TE identifer
          ident(n)=-maxpseudoTE   ! TE identifier
        endif
      endif

c --------------------------------------------------------------------------
c --- flag primary detector not to be the outlayer
c --------------------------------------------------------------------------
      if(useoutlayer)then
        if(fsfdbg.gt.0) then
          print*,'FSFITS: Fitting with outlayer.'
        endif
        if(pd.eq.i_TPC)then         ! TPC is primary detector
          ioutr = 2 ! do not throw out TPC
        else if(pd.eq.i_ITC)then    ! ITC is primary detector
          ioutr = 1 ! do not throw out ITC
        else if(pd.eq.i_FTD)then    ! FTD is primary detector
          ioutr = 8 ! do not throw out FTD
        else if(pd.eq.i_SCAN1)then  ! this is string primary
          ioutr = 3 ! anything can be thrown out
        else if(pd.eq.i_SCAN2)then  ! this is string primary
          ioutr = 3 ! anything can be thrown out
        else if(pd.eq.i_SCAN3)then  ! this is string primary
          ioutr = 3 ! anything can be thrown out
        else if(pd.eq.i_SCAN4)then  ! this is string primary
          ioutr = 3 ! anything can be thrown out
        else if(pd.eq.i_VTX)then     ! VTX is primary detector
          ioutr = 3 ! anything can be thrown out
c --- -ME- always crash if outlayer and primary detector is zero !
        else
          print*,'FATALERROR in FSFITS: Invalid primary detector given!'
          print*,'                      This *must* not happen.'
          stop
        endif
      else
        ioutr = 0
      endif

c --------------------------------------------------------------------------
c --- fit track
c --------------------------------------------------------------------------
      IF (useextrapol) THEN
         call fk3trk(ioutr,n,ident,idatte,space(its+1),
     1               tk,nout,idout,ier1,ier2,iflte,ierx)
      ELSE
         call fk3trk(1000+ioutr,n,ident,idatte,space(its+1),
     1               tk,nout,idout,ier1,ier2,iflte,ierx)
      END IF

c --------------------------------------------------------------------------
c --- kill this TS ?
c --------------------------------------------------------------------------
      if(ier1.ne.0)then
        if(fsfdbg.gt.0) then
          print*,'FSFITS: TS killed - fit not possible (ier1=',ier1,')'
        endif
        its = 0
        return
      endif

c --------------------------------------------------------------------------
c --- check that track parameters make any sense
c --------------------------------------------------------------------------
      if (tk(20).lt.0. .or. tk(20).gt.pi .or.
     1    tk(21).lt.0. .or. tk(21).gt.2*pi .or.
     1    tk(10).lt.0.) then
        if (fsfdbg.gt.0) then
          print *,'WARNING in FSFITS: fitted theta/phi/length is off'
          print *,'TS dump: ',(itk(k),k=1,8),(tk(k),k=9,22)
c --- debug fit
          print *,'Refit with debug and see:'
          ikeep      = ideb
          ideb = 2
CC          call fk3trk(ioutr,n,ident,idatte,space(its+1),
          call fk3trk(1000+ioutr,n,ident,idatte,space(its+1),
     1                tk,nout,idout,ier1,ier2,iflte,ierx)
          ideb = ikeep
          print *,'Drop it !'
        endif
        its = 0
        RETURN
      endif

c --------------------------------------------------------------------------
c --- remove all outlayers
c --------------------------------------------------------------------------
      if (nout.gt.0.and.useoutlayer)then
        if(fsfdbg.gt.0) then
          print*,'FSFITS: # of outlayers: ',nout,' of ',ispace(its+5)
        endif
        if(nout.eq.n)then
          print*,'FSFITS: TS killed - all TEs rejected'
c --- kill this TS
          its = 0
          return
        endif
        if (fsfdbg.gt.0.and.DELSIM)then
c --- get label for the TS
          call fsftlla (n,ident,ilabel)
c --- get label of TE
          if (ilabel.gt.0) then
            do i=1,nout
c --- protect against pseudo TEs being outlayer
              if (idout(i).gt.0) then
                CALL TKGLAB('TE',idout(i),labelte,IERR)
                if (labelte.eq.ilabel) then
                  print '(A,A,I8,A,I8)',
     1                  ' WARNING in FSFITS: good label outlayer,',
     1                  ' TE id ',idout(i),' label ',ilabel
c --- find TER data of this TE
                  do j = 1,n
                    if (ident(j).eq.idout(i)) then
                      print *,' TE module : ',idatte(1,j),
     1                        ' sub-module : ',idatte(2,j)
                    endif
                  enddo
                endif
              endif
            enddo
          endif
        endif
        do i=1,nout
c --- primary vertex can never be a out for IDVD fit !
          if (idout(i).eq.-maxpseudoTE) then
            if(fsfdbg.gt.0) then
              print*,'FSFITS: pseudo-TE fudged is outlayer'
            endif
          else
            if (idout(i).lt.0) then
              if(fsfdbg.gt.0) then
                print*,'FSFITS: pseudo-TE is outlayer'
              endif
            endif
c --- ok, normal job
            if(fsfdbg.gt.1) then
              print'(A,i8)',' FSFITS: TE rejected:',idout(i)
            endif
c --- limit to the original TEs
            j = ispace(its+5)
            do while (ident(j).ne.idout(i).and.j.gt.1)
              j = j-1
            enddo
            if (ident(j).eq.idout(i))then
              ispace(its+tssize+j)=ispace(its+tssize+ispace(its+5))
              ident(j)            =ident(ispace(its+5))   ! bugfix
              ispace(its+5)       =ispace(its+5)-1
            else
              print '(A,i8)',
     1         ' WARNING in FSFITS: rejected TE not found, id:',idout(i)
            endif
          endif
        enddo
      endif

c --------------------------------------------------------------------------
c --- copy fitted data back to TS; pay attention to the different
c     contents of the measurement code and to the number of TEs in
c     string.
c --------------------------------------------------------------------------
      ispace(its+1)=itk(1)
c --- check the detector bit code given back by FK3TRK   -KH-
      if (fsfdbg.gt.1) then
        write(6,*) 'FSFITS: detector bit code before FK3TRK: ',
     >             ispace(its+2)
        write(6,*) 'FSFITS: detector bit code after  FK3TRK: ',
     >             itk(2)
      end if
      ispace(its+2)=itk(2)
      if(1.eq.jbit(itk(3),1))then
        ispace(its+3)=7          ! cylinder parameters
      else
        ispace(its+3)=5          ! xyz
      endif
      ispace(its+4)=itk(4)
      call ucopy(tk(6),space(its+6),mxtk-5)

c --- reset primary vertex bit(s) if needed
      call sbit0(ispace(its+2),detid_pvtx  +1)
      call sbit0(ispace(its+2),detid_zfudge+1)
CC detid_shower is TPC !
CC      call sbit0(ispace(its+2),detid_shower+1)

c --------------------------------------------------------------------------
c --- if beam constraint, give track at R=5cm for internal consistency
c --------------------------------------------------------------------------
      if (ispace(its+3).eq.5) THEN
        R_TS = sqrt(space(its+17)**2+space(its+18)**2)
      else
        R_TS = space(its+17)
      endif
c     if (R_TS.LT.4.5) THEN   -KH- our beampipe is smaller!
      if (R_TS.LT.1.0) THEN
c --- ok, extrapolate back to surface of beam pipe
        if(fsfdbg.gt.1) then
          print*,'FSFITS: fit using beam spot/pvtx - give TS at R=1cm'
        endif
c --- extrapolate TS
        call fsftsexx(its,1,'CYLI',1.0,-1000.,1000.,textra,ierrsf,ier)
        if (ier+ierrsf.ne.0) then
          if (ierr.ne.0) then
            print*,'WARNING in FSFITS: exxtrk failed - drop TS !'
          endif
          its = 0
          return
        endif
c --- copy back result
        ispace(its+3)=7          ! cylinder parameters
        call ucopy (textra(7),space(its+17),21)
c --- check covariance matrix for possible corruption
        call fkmi5 (space(its+23),cov,ierr)
c --- reject TS if cov is crap
        if (ierr.ne.0) then
          if(fsfdbg.gt.0) then
            print*,
     1  'WARNING in FSFITS: TS with beam spot/pvtx killed - cov is crap'
          endif
c --- kill this TS
          its = 0
          return
        endif
      endif

c(KH) now transform TS parameters to inside beam pipe if not already done
c     (extrapolation was done automatically during track fit. just copy params)
      if (r_ts.gt.1.0.and.its.ne.0.and.iflx(1).eq.0.and.useextrapol)then
         do i=1,21
            space(its+16+i)=datx(i,1)
         enddo
         call sbit1(ispace(its+3),1)
      endif

      end

