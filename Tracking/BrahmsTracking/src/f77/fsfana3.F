*  
* $Id: fsfana3.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfana3(firstdet,method,flag)
*******************************************************************************
*
*      New steering routine to find strings.
*
*      INPUT:  firstdet  - detector index to start with
*              method    - type of pairsearch in FSFMKPAIR
*              flag      - flags to steer methods
*                                                                 Markus Elsing
*******************************************************************************
      implicit none
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfflags.inc"
#include "include/fsfreject.inc"
c ---
      integer firstdet,method(detmax),jbit,flag
      logical fsfmkpair,fsfftdm,fsfftds,fsfvtxm,fsfvtxs
c
      integer i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,idxv(detmax)
      integer ifirst,ibitted,ibad,sumhits,ite,ite2,idet,numdet
      integer det(detmax),idetused,iwant,ihit
      integer index(0:temax,detmax),nindex(0:detmax)
      logical uselabl,useused,usesingle
      logical accepted,hit(detmax),skipbad,usepconst,scandet,cut100
c
      integer ITEOUT(8),NOUT,ierr,found,k1,k2,prodhits

c --- decode FLAGs
      uselabl    =jbit(flag,FBIT_LABL).eq.1
      usepconst  =jbit(flag,FBIT_PCONST).eq.1
      useused    =jbit(flag,FBIT_USEUSED).eq.1
      usesingle  =jbit(flag,FBIT_GIVESINGLTE).eq.1
      scandet    =jbit(flag,FBIT_SCANDET).eq.1
      cut100     =jbit(flag,FBIT_CUT100).eq.1

c --- Copy TEs to arrays and adjust errors.
      call fsftec(firstdet,usepconst,method)

c --- now loop over TEs of primary detector and try to find strings
      do ifirst=1,teindexused(firstdet)

c --- Dont use this TE, if it is already used.(teused is set, if it is already
c --- used.) Furthermore require that the TE has some radius, either from the
c --- TE itself or from the beam constraint
        if((useused.or..NOT.teused(ifirst,firstdet)).and.
     1     (jbit(datdescr(ifirst,firstdet),10).eq.1 .or.
     2      jbit(datdescr(ifirst,firstdet),11).eq.1))then

c --- some debug output
          if(fsfdbg.gt.0) then
            print'(A,A5,A,I8,A,I8)',
     1           ' FSFANA3: start with ',detname(firstdet),
     2           ' TE id. ',teid(ifirst,firstdet),
     3           ' label ',telabl(ifirst,firstdet)
          endif
c ------------------------------------------------------------------------
c --- preselect all TEs which give a pair with primary TE
c ------------------------------------------------------------------------
c --- reset index arrays for TEs paired; set fixed value for the firstdet !
          call vzero (nindex,detmax+1)
          nindex(firstdet)  = 1      ! number of TEs in det(1)==firstdet is 1
          index(1,firstdet) = ifirst ! thats the index in TE arrays
c --- now create list of TEs in detectors which give a pair with primary TE
          do idet = 1,detmax
c --- do we try pairs with this detector ?
            if (idet.ne.firstdet.and.method(idet).ne.0) then
c --- do we have to scan for TEs in pseudo-TE from strings ?
              ite = 0
              if (scandet) THEN
c --- do not allow for others than simple detectors, no pseudo TEs
                if (detid(idet).eq.0) then
                  ite = 1
c --- submodule is detector bitted word for TS in TE !
                else if (jbit(ispace(teindex(ifirst,firstdet)+2),
     1                        detid(idet)+1).eq.1) then
                  ite = 1
                endif
              endif
              if (ite.eq.0) then
c --- loop over all TEs
                do ite = 1,teindexused(idet)
c --- if we have a pair, store index of this TE in array
                  if(fsfmkpair(method(idet),ifirst,firstdet,
     1                         ite,idet))then
                    nindex(idet)             = nindex(idet)+1
                    index(nindex(idet),idet) = ite
                  endif
                enddo
              endif
            endif
          enddo
c ------------------------------------------------------------------------
c --- now do bookkeeping
c ------------------------------------------------------------------------
 100      CONTINUE              ! dirty jump
          call vzero(det,detmax)
          idetused = 1
          det(1)   = firstdet
          ibitted  = 0
          sumhits  = 0
          prodhits = 1
          if(fsfdbg.gt.0)print *,'---------------------------'
          do idet = 1,detmax
            sumhits  = sumhits + nindex(idet)
            if(nindex(idet).gt.0)then
              prodhits = prodhits * nindex(idet)
c --- write debug output
              if(fsfdbg.gt.0)then
                print*,'FSFANA3: number of TEs in ',detname(idet),
     1                 ' successfully paired is ',nindex(idet)
              endif
c --- Renumber the detectors WITH HITS to have firstdet==det(1)
c     and det(2..detmax) containing each other detector once, THAT HAS HITS
              if (idet.ne.firstdet) then
                idetused      = idetused+1
                det(idetused) = idet
              endif
c --- set up bitted word for test of possible good TS combination
              if (detid(idet).gt.0)then
                call sbit1(ibitted,detid(idet)+1)
              else if(idet.eq.firstdet) then
c --- use copy of TS measurement code
                ibitted = ispace(teindex(ifirst,firstdet)+2)
              endif
            endif
          enddo
c --- look if any combination possible is not in BADCOMB
          skipbad = .FALSE.
          if (.not.usesingle) then
            ibad    = 1
            do while(badcomb(ibad).ne.0)
              if(ibitted.eq.badcomb(ibad))then
                skipbad = .TRUE.
              endif
              ibad=ibad+1
            enddo
          endif
c
c --- now decide if we should try to find a string for this TE
c
CC          print *,'FSFANA: prodhits :', prodhits

          if (sumhits.eq.1.and..not.usesingle)then
            if(fsfdbg.gt.0) print*,'FSFANA3: reject TE - no hits'
          else if (prodhits.gt.100.and.cut100) then
            if (fsfdbg.gt.0) then
              print*,'FSFANA3: more than 100 combinations, cut down'
            endif
            do idet = 1,detmax
              if (nindex(idet).gt.1) then
                if (fsfdbg.gt.0) then
                  print*,'FSFANA3: => drop ',detname(idet)
                endif
                nindex(idet) = 0
              endif
            enddo
            goto 100
          else if (skipbad) then
            if(fsfdbg.gt.0) then
              print*,'FSFANA3: reject TE - all combinations in badcomb'
            endif
          else
c ------------------------------------------------------------------------
c --- Now try polar inversion for each combination of hits
c ------------------------------------------------------------------------
c --- We always have a hit in 'firstdet==det(1)'. That's where we start.
c     hit(det)=.FALSE. shall indicate: this detector has a TE in a accepted TS
C(KH)       call vzero(hit,detmax)
            do iwant=1,detmax
              hit(iwant)=.false.
            enddo
            ihit          = 1
            hit(firstdet) = .TRUE.
c --- start with number of detectors wanted = number of detectors paired
            iwant = idetused
c --- look for iwant detectors in string
            do while (iwant.ge.ihit.and.(iwant.gt.1.or.usesingle))
              if(fsfdbg.gt.0)then
                print*,'FSFANA3: want ',iwant,' detectors, ihit=',ihit
                print*,'======================='
              endif

c ------------------------------------------------------------------------
c --- loops over all index lists, innermost index is firstdet to save CPU !
c ------------------------------------------------------------------------
c --- set up idxv array to give to subroutines
              call vzero(idxv,detmax)
              idxv(det(1))  = ifirst             ! this is always the case
c --- now the loops !
              do 1012 i12  = nindex(det(12)),0,-1
              if(det(12).ne.0)then
                if(i12.eq.0.and.hit(det(12))) GOTO 1012
                idxv(det(12)) = index(i12,det(12))
              endif
              do 1011 i11  = nindex(det(11)),0,-1
              if(det(11).ne.0)then
                if(i11.eq.0.and.hit(det(11))) GOTO 1011
                idxv(det(11)) = index(i11,det(11))
              endif
              do 1010 i10  = nindex(det(10)),0,-1
              if(det(10).ne.0)then
                if(i10.eq.0.and.hit(det(10))) GOTO 1010
                idxv(det(10)) = index(i10,det(10))
              endif
              do 1009 i9  = nindex(det(9)),0,-1
              if(det(9).ne.0)then
                if(i9.eq.0.and.hit(det(9))) GOTO 1009
                idxv(det(9)) = index(i9,det(9))
              endif
              do 1008 i8  = nindex(det(8)),0,-1
              if(det(8).ne.0)then
                if(i8.eq.0.and.hit(det(8))) GOTO 1008
                idxv(det(8)) = index(i8,det(8))
              endif
              do 1007 i7  = nindex(det(7)),0,-1
              if(det(7).ne.0)then
                if(i7.eq.0.and.hit(det(7))) GOTO 1007
                idxv(det(7)) = index(i7,det(7))
              endif
              do 1006 i6  = nindex(det(6)),0,-1
              if(det(6).ne.0)then
                if(i6.eq.0.and.hit(det(6))) GOTO 1006
                idxv(det(6)) = index(i6,det(6))
              endif
              do 1005 i5  = nindex(det(5)),0,-1
              if(det(5).ne.0)then
                if(i5.eq.0.and.hit(det(5))) GOTO 1005
                idxv(det(5)) = index(i5,det(5))
              endif
              do 1004 i4  = nindex(det(4)),0,-1
              if(det(4).ne.0)then
                if(i4.eq.0.and.hit(det(4))) GOTO 1004
                idxv(det(4)) = index(i4,det(4))
              endif
              do 1003 i3  = nindex(det(3)),0,-1
              if(det(3).ne.0)then
                if(i3.eq.0.and.hit(det(3))) GOTO 1003
                idxv(det(3)) = index(i3,det(3))
              endif
              do 1002 i2  = nindex(det(2)),0,-1
              if(det(2).ne.0)then
                if(i2.eq.0.and.hit(det(2))) GOTO 1002
                idxv(det(2)) = index(i2,det(2))
              endif

c ------------------------------------------------------------------------
c --- loop and count number of detectors we have
c ------------------------------------------------------------------------
                numdet  = 0
                do idet=1,detmax
                  if (idxv(idet).ne.0) numdet = numdet + 1
                enddo
c --- do we have correct number of detectors ?
                if (numdet.ne.iwant) goto 1002 ! ---> next combination
c --- is this combination possible ?
                if (numdet.gt.1) then
                  call fsftte(idxv,badcomb,0,skipbad)
                  if (skipbad) goto 1002 ! ---> next combination
                endif

c --- what about labels ?
                labelsok= .TRUE.
                label   = 0
                do idet=1,detmax
                  if(idxv(idet).ne.0)then
                    if(label.eq.0)then
                      label = abs(telabl(idxv(idet),idet))
                      if(label.eq.0) labelsok=.FALSE.
                    else if(label.ne.abs(telabl(idxv(idet),idet)))then
                      labelsok= .FALSE.
                      label   = 0
                    endif
                  endif
                enddo
                if (label.eq.0.or.numdet.le.1.or.
     1              telabl(idxv(det(1)),det(1)).lt.0) then ! label > 0
                  labelsok=.FALSE.
                endif
c --- debug output to let me know
                if(fsfdbg.gt.1.or.
     1             (fsfdbg.gt.0.AND.labelsok))then
                  print*,'FSFANA3: try this combination:'
                  do idet=1,detmax
                    if (idxv(idet).ne.0)then
                      print'(A,i4,A,A5,A,I8,A,i8)',
     1                '      TE # ',idxv(idet),
     2                ' from ',detname(idet),
     3                ' id #',teid(idxv(idet),idet),
     4                ' label ',telabl(idxv(idet),idet)
                    endif
                  enddo
                  if (labelsok)then
                    print*,'FSFANA3: good combination of labels'
                  endif
                endif
c ------------------------------------------------------------------------
c --- now call the string finding engine
c ------------------------------------------------------------------------
                call fsfplinv(idxv,det(1),flag,accepted)

                if(.not.accepted)then
                  if(fsfdbg.gt.0)then
                    if(labelsok)then
                      print*,'WARNING in FSFANA3: ',
     1                       'good labels not accepted'
                    else
                      print*,'FSFANA3: was not accepted'
                    endif
                    print*,'======================='
                  endif
                  goto 1002 ! ---> next combination
                endif

c --- If this combination was accepted:
c     Mark the involved detectors to have a hit.
                if(fsfdbg.gt.0)then
                  print*,'FSFANA3: combination was accepted.'
                  print*,'======================='
                endif

c ------------------------------------------------------------------------
c --- we have to do again some bookkeeping, does not find VD
c ------------------------------------------------------------------------
c
c --- loop over all TEs in TS and try to find it in local arrays
                do 1000 ite = 1,ispace(tsindex(tsindexused)+5)
c --- now loop over all detectors
                  do idet = 1,detmax
c --- check bit to speed up
                    if (jbit(ispace(tsindex(tsindexused)+2),
     1                       detid(idet)+1).eq.1) then
c --- loop over all TEs and look for this identifier
                      do ite2 = 1,teindexused(idet)
                        if (ispace(tsindex(tsindexused)+tssize+ite).eq.
     1                      teid(ite2,idet)) then
c --- mark TE as used
                          teused(ite2,idet)=.TRUE.
c --- special case for VTX/FTD, bookkeeping for multipoints
                          if(FSFFTDS(ispace(teindex(ite2,idet)+1)))then
c --- check, if we where looping on a FTD multipoint
                            if(idxv(i_FTD).ne.0.and..not.hit(i_FTD))then
                              teused(idxv(i_FTD),i_FTD)=.TRUE.
c --- set hit=.TRUE. for this detector
                              ihit      = ihit + 1
                              hit(i_FTD) = .TRUE.
                              if(fsfdbg.gt.0)then
                                print*,'FSFANA3: Setting hit(det)=',
     1                               '.TRUE. det=',i_FTD,
     2                               '(',detname(i_FTD),'), ihit=',ihit
                              endif
                            endif
                          endif
                          if(FSFVTXS(ispace(teindex(ite2,idet)+1)))then
c --- check, if we where looping on a VTX multipoint
                            if(idxv(i_VTX).ne.0.and..not.hit(i_VTX))then
                              teused(idxv(i_VTX),i_VTX)=.TRUE.
c --- set hit=.TRUE. for this detector
                              ihit      = ihit + 1
                              hit(i_VTX) = .TRUE.
                              if(fsfdbg.gt.0)then
                                print*,'FSFANA3: Setting hit(det)=',
     1                               '.TRUE. det=',i_VTX,
     2                               '(',detname(i_VTX),'), ihit=',ihit
                              endif
                            endif
                          endif
c --- do we have this detector in idxv and no hit ?
                          if (idxv(idet).ne.0.and..not.hit(idet))then
c --- set hit=.TRUE. for this detector
                            ihit      = ihit + 1
                            hit(idet) = .TRUE.
                            if(fsfdbg.gt.0)then
                              print*,'FSFANA3: Setting hit(det)=',
     1                               '.TRUE. det=',idet,
     2                               '(',detname(idet),'), ihit=',ihit
                            endif
                          endif
c --- next TE
                          goto 1000 ! next TE from TS
                        endif
                      enddo
                    endif
                  enddo
 1000           continue
c ------------------------------------------------------------------------
c --- end loops on detectors !
c ------------------------------------------------------------------------
 1002         continue ! i2
 1003         continue ! i3
 1004         continue ! i4
 1005         continue ! i5
 1006         continue ! i6
 1007         continue ! i7
 1008         continue ! i8
 1009         continue ! i9
 1010         continue ! i10
 1011         continue ! i11
 1012         continue ! i12

c ------------------------------------------------------------------------
c --- next time one detector less !
c ------------------------------------------------------------------------
              iwant = iwant - 1
              if(fsfdbg.gt.0.and.iwant.lt.ihit)then
                print*,'FSFANA3: ihit is ',ihit,
     1                 ' --> dont look for less hits !'
                print*,'======================='
              endif
            enddo ! iwant
          endif ! hits found with primary TE
        endif ! good primary TE
      enddo ! primary TEs

      return
      end


