*  
* $Id: fsfwts.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfwts(idxv,parami,errori,pd,flag,ierr)
*******************************************************************************
*      This subroutine will create a TS like data structure in the
*      space array and fill tsindex() accordingly.
*
*      INPUT
*        idxv    - list of TEs
*        parami  - track parameters
*        errori  - diagonal error matrix elements
*        pd      - index of primary detector (used for HISTOs only)
*        flag    - bitted word to flage methods used to create TS
*      OUTPUT:
*        ierr    - set to 1 in case of a problem
*
*                                                                  Daniel Wicke
*******************************************************************************
      implicit none
#include "include/fsfb.inc"
#include "include/brcntl.inc"
#include "include/fsfpi.inc"
#include "include/fsfparams.inc"
#include "include/fsfflags.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
c#include "include/gcflag.inc"
      integer pd,flag,k
      integer idxv(detmax),itmp,id,od
      real    parami(6),errori(6),paramf(6),dat(MXTE)
      integer exxerr,JBYT,idat(MXTE),deldet

      equivalence (dat,idat)

      integer i,j,n,ierr,jbit,ntes
      logical uselabl,useused,givesnglte,fsfgetspace,fsfftdm
      integer iteout(8),nout,inum(DETMAX),delid

      ierr = 0

c --- catch theta or phi out of range
      if (parami(4).lt.0.or.parami(4).gt.pi.or.
     1    parami(5).lt.0.or.parami(5).gt.2*pi) then
        print*,'WARNING in FSFWTS: theta/phi off'
        print*,'parameters :',(parami(i),i=1,6)
        print*,'drop it !'
        ierr = 1
        return
      endif

      id = pd
      od = pd

c --- Catch fatal errors
      if (.not.fsfgetspace('FSFWTS  ',tssize+maxtev)) THEN
        ierr = 1
        return
      endif
      if(tsindexused.ge.tsmax)then
        print*,'ERROR in FSFWTS: ',
     1               '**** TSMAX to small to store this TS ****'
        print*,'            Enlarge parameter tsmax recompile'
        print*,'            and rerun to cure this problem!!'
        ierr = 1
        CALL TKVETO(1)
        return
      endif
      if(parami(4).eq.0)then
        print*,'ERROR in FSFWTS: Theta == 0. Track not written as TS!'
        ierr = 1
        return
      endif

      uselabl    =jbit(flag,FBIT_LABL).eq.1
      useused    =jbit(flag,FBIT_USEUSED).eq.1
      givesnglte =jbit(flag,FBIT_GIVESINGLTE).eq.1

c --- zero array to be sure that we have no junk
      call vzero(ispace(spaceused+1),tssize)
c --- Now start filling the TS structure
      ispace(spaceused+1)=3000+pd !The module identifier
      ispace(spaceused+2)=0   ! The used detectors will be marked beyond.
c --- do not give start/end in the bits !
      ispace(spaceused+3)=1+8 ! Crude estimation, plane reference, no start/end
      ispace(spaceused+4)=0   ! we do later

c --- Count number of found TEs,
c     set corresponding bits in TS measurement code
c     as well as in itmp (index for our own statistics).
c     and fill in teid vector
      n    = 0
      itmp = 0
      do i=1,detmax
        if(idxv(i).ne.0)then

* --- write V0/conv/hadr.shower bits to bits 1-3 here !
          if (showerte(idxv(i),i).ne.0)
     1      ispace(spaceused+4)=JBYT(showerte(idxv(i),i),17,3)

c --- test if inner and outer detector
          if (R(idxv(i),i).lt.R(idxv(id),id))then
            id = i
          else if (R(idxv(i),i).gt.R(idxv(od),od))then
            od = i
          endif
c
c --- new code for strings in dummy TEs
c
          if(ispace(teindex(idxv(i),i)+1).gt.3000)then
c --- get number of TEs in this TS given as TE
            ntes = tesize(idxv(i),i)-31 ! 31 is default here
            if(fsfdbg.gt.1)then
              print '(A,A,A,i5,A,i8,A,i8)',' FSFWTS: dummy TE of ',
     1              detname(i), ' #',idxv(i),' id.',teid(idxv(i),i),
     2              ' # of TEs ',ntes
            endif
c --- copy TE identifiers
            if (ntes.le.0)then
              print*,
     1         'WARNING in FSFWTS: dummy TE for string corrupted ?'
            else
              do j=1,ntes
                n=n+1
                ispace(spaceused+tssize+n) =
     1            ispace(teindex(idxv(i),i)+31+j)
                if(fsfdbg.gt.0)then
                  print '(A,i8)',' FSFWTS: --> TE id. ',
     1                  ispace(teindex(idxv(i),i)+31+j)
                endif
              enddo
c --- detector bitted word is still in TE !
              ispace(spaceused+2) = ior(ispace(spaceused+2),
     1                              ispace(teindex(idxv(i),i)+2))

              if (fsfdbg.gt.1) then
                write(6,*) 'FSFWTS: added bits',
     >          (JBIT(ispace(teindex(idxv(i),i)+2),k),k=32,1,-1),
     >                     ' to TS submod.'
                write(6,*) '        this comes from TE ',idxv(i),
     >                     'of subdet ',detname(i)
                write(6,*) '        TE ID   : ',teid(idxv(i),i)
                write(6,*) '        TE index: ',teindex(idxv(i),i)
              endif

            endif
c --- replace FTD and VTX multipoints
          else if (i.eq.i_FTD.or.i.eq.i_VTX) then
c -- replace multipoint by exclusion list
            call tklhte (teid(idxv(i),i),10,ntes,inum,ierr)
            if (ntes.eq.0) then
              write(6,*) 'FSFWTS: no hit TEs for VTX/FTD found.'
              write(6,*) '        this is most likely caused by'
              write(6,*) '        exclusion list overflow. this'
              write(6,*) '        severely affects track reconstruction'
            endif
c --- loop over new TEs and set bits correct
            do j=1,ntes
              if(fsfdbg.gt.0)then
                print '(A,i8)',' FSFWTS: --> TE id. ',inum(j)
              endif
              call TKREAD('TE',0,inum(j),dat,idat,ierr)
              if (ierr.ne.0) then
                stop 'FSFWTS: TKREAD error !'
              endif
              if(delid(idat(1)).eq.i_VTX.or.delid(idat(1)).eq.i_FTD)then
                if(fsfdbg.gt.0)then
                  print *,' => MULTIPUNKT replacement, drop !!!'
                endif
              else
                n=n+1
                ispace(spaceused+tssize+n) = inum(j)
                call sbit1(ispace(spaceused+2),deldet(idat(1))+1)
                if (fsfdbg.gt.1) then
                 write(6,*) 'FSFWTS: added detid ',deldet(idat(1))
                endif
              endif
            enddo
c --- old code for normal TEs
          else
            n=n+1
            ispace(spaceused+tssize+n)=teid(idxv(i),i)
            call sbit1(ispace(spaceused+2),
     1           ispace(teindex(idxv(i),i)+1)+1) !!! BUG !!!
            if(fsfdbg.gt.1)then
              print '(A,A,A,i5,A,i8,A,i8)',' FSFWTS: TE of ',
     1             detname(i),' #',idxv(i),' id.',teid(idxv(i),i),
     2             ' sub.',ispace(teindex(idxv(i),i)+2)
              write(6,*) '   TS submodule bit ',
     >             ispace(teindex(idxv(i),i)+1)+1,' set.'
            endif
          endif
        endif
      enddo

      if(fsfdbg.gt.1)then
        print*,'FSFWTS: Tracktype created:',
     1          ispace(spaceused+4),'; number of TEs: ',n
        write(6,*) '        Detector code    : ',ispace(spaceused+2)
      endif

      ispace(spaceused+5)=n    !number of TEs in this string

      ispace(spaceused+6)=3    !unknown charge
                               !(charge may change later on, when the fit is
                               ! called anyway)

      ispace(spaceused+7)=0    ! mass code
      ispace(spaceused+8)=0    ! NDF
       space(spaceused+9)=0    ! chi squared

      call ucopy (parami,paramf,6)
c --- First extrapolate to innermost detector:
      if (n.gt.1) then
        call fsfxts(parami,paramf,idxv(id),id,exxerr)
        if(exxerr.ne.0)then
          if(fsfdbg.gt.0)then
            print*,'WARNING in FSFWTS: Reference point not given at',
     1                                 ' innermost detector.'
            print*,'WARNING in FSFWTS: also startpoint of TS incorrect.'
          endif
        endif
      endif

c --- Reference point
      space(spaceused+17)=paramf(2) !x
      space(spaceused+18)=paramf(3) !y
      space(spaceused+19)=paramf(1) !z
c --- Track segment parameter
      space(spaceused+20)=paramf(4) !theta
      space(spaceused+21)=paramf(5) !phi
      space(spaceused+22)=paramf(6)*sin(paramf(4))/bfieldz
                        !parami(6)==1/R_T; space(spaceused+22) should be 1/p .
c --- Diagonal Error Matrix
      call vzero(ispace(spaceused+23),15)
      space(spaceused+23)=errori(2)
      space(spaceused+24)=errori(3)
      space(spaceused+25)=errori(4)
      space(spaceused+26)=errori(5)
      space(spaceused+27)=errori(6)
c --- no Startpoint (same as reference point)
      space(spaceused+11)=0.
      space(spaceused+12)=0.
      space(spaceused+13)=0.
c --- no Endpoint
      space(spaceused+14)=0.
      space(spaceused+15)=0.
      space(spaceused+16)=0.
c --- no Track length
      space(spaceused+10)=0.

c --- update pointers

      tsindexused         =tsindexused+1

c(KH) ---check range
      if (tsindexused.gt.tsmax) then
        write(6,*) 'FSWTS: TSMAX too small.'
        tsindexused=tsmax
        CALL TKVETO(1)
      end if

      tsindex(tsindexused)=spaceused
      spaceused           =spaceused+tssize+maxtev

      if (fsfdbg.gt.1) then
        write(6,*) 'FSFWTS: tsindexused increased to ',tsindexused
      end if

      end

