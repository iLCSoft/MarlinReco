/*  
* $Id: fsfinp.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*/  
      subroutine fsfinp
*******************************************************************************
*
*     Read all TEs into local arrays to prepare search.
*
*
*     WRITTEN BY                                     Daniel Wicke
*     CHANGES BY                                     Kristian Harder, August 99
*
*******************************************************************************
      implicit none
#include "include/fkfild.inc"
#include "include/fsfclun.inc"
#include "include/brcntl.inc"
#include "include/fsfb.inc"
#include "include/fsfpi.inc"
#include "include/fsfsim.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfxerrs.inc"
#include "include/detnums.inc"

c +SEQ, TEBITS.   ! user bits tanagra
      integer i,j,tknumb
      integer ierr,nid,detindex
      logical fsfgetspace

c --- for the statistics
      integer n_unknwn,n_singleFTD,n_multiFTD,n_singleVTX,n_multiVTX

c --- reset statistics
      n_unknwn  = 0
      n_singleFTD = 0
      n_multiFTD  = 0
      n_singleVTX = 0
      n_multiVTX  = 0

c --- Get magnetic field at origin from trackfit
      bfieldz=consb
      if(fsfdbg.gt.0)then
        print*,'FSFINP: BField used by DELSEARCH is',bfieldz
      endif

c --- Reset FSFARRAYS:
      spaceused   = 0
      tsindexused = 0
      do i=1,detmax
        teindexused(i)=0
        do j=1,temax
          teused(j,i)=.FALSE.
        enddo
      enddo

c --- Fill FSFARRAYS:
      NID=TKNUMB('TE',0)
      if(fsfdbg.gt.1)then
        print *,'FSFINP: number of TEs in event: ', nid
      endif
      do i=1,nid

c --- do we have the space needed ?
        if(.not.fsfgetspace('FSFINP  ',MXTE)) RETURN

c --- copy TE from BRAHMS array into FSF array         -KH-

        CALL TKREAD('TE',0,I,SPACE(SPACEUSED+1),ISPACE(SPACEUSED+1),
     +       IERR)

        IF (IERR.NE.0) THEN
          WRITE(6,*) 'FSFINP: TKREAD FAILED.'
          WRITE(6,*) 'TE number ',I,' of ',NID
          STOP
        END IF

c
c --- get DELSEARCH detector index for this detector
c     once we have distinct IDs for all VTX/FTD layers, this is easily
c     done by DELDET in TKSTEER.
c
        if(ispace(spaceused+1).eq.ID_FTD)then
          detindex=i_FTD
          n_multiFTD  = n_multiFTD + 1
        else if(ispace(spaceused+1)/100.eq.ID_FTD/100)then
          n_singleFTD  = n_singleFTD + 1
          if(ispace(spaceused+1).eq.ID_FTD+1)then
            detindex=i_FTD1
          else if(ispace(spaceused+1).eq.ID_FTD+2)then
            detindex=i_FTD2
          else if(ispace(spaceused+1).eq.ID_FTD+3)then
            detindex=i_FTD3
          else if(ispace(spaceused+1).eq.ID_FTD+4)then
            detindex=i_FTD4
          else if(ispace(spaceused+1).eq.ID_FTD+5)then
            detindex=i_FTD5
          else if(ispace(spaceused+1).eq.ID_FTD+6)then
            detindex=i_FTD6
          else if(ispace(spaceused+1).eq.ID_FTD+7)then
            detindex=i_FTD7
          else
            detindex=0
            n_unknwn=n_unknwn+1
          endif
        else if(ispace(spaceused+1).eq.ID_VTX)then
          detindex=i_VTX
          n_multiVTX  = n_multiVTX + 1
        else if(ispace(spaceused+1)/100.eq.ID_VTX/100)then
          n_singleVTX  = n_singleVTX + 1
          if(ispace(spaceused+1).eq.ID_VTX+1)then
            detindex=i_VTX1
          else if(ispace(spaceused+1).eq.ID_VTX+2)then
            detindex=i_VTX2
          else if(ispace(spaceused+1).eq.ID_VTX+3)then
            detindex=i_VTX3
          else if(ispace(spaceused+1).eq.ID_VTX+4)then
            detindex=i_VTX4
          else if(ispace(spaceused+1).eq.ID_VTX+5)then
            detindex=i_VTX5
          else
            detindex=0
            n_unknwn=n_unknwn+1
          endif
        else if(ispace(spaceused+1).eq.ID_VTXP)then
          detindex=i_VTX
          n_multiVTX  = n_multiVTX + 1
        else if(ispace(spaceused+1).eq.ID_ITC)then
          detindex=i_ITC
        else if(ispace(spaceused+1).eq.ID_SIT+1)then
          detindex=i_SIT1
        else if(ispace(spaceused+1).eq.ID_SIT+2)then
          detindex=i_SIT2
        else if(ispace(spaceused+1).eq.ID_TPC)then
          detindex=i_TPC
        else if(ispace(spaceused+1).eq.ID_FCH)then
          detindex=i_FCH
        else
          detindex=0
          n_unknwn=n_unknwn+1
        endif

c --- convert detector ID from BRAHMS convention into DELSEARCH convention
        call fsfconvrt(ispace(spaceused+1))

c
c --- copy into arrays if TE is to be processed
c
        if(detindex.ne.0)then
          if (fsfdbg.gt.1) then
            write(6,*) 'FSFINP: saving TE ',i,' of length ',
     >                 ispace(spaceused+5)
          end if
          call fsfwte ('FSFINP',detindex,i,ispace(spaceused+5),ierr)
          if (ierr.ne.0) RETURN
        endif
      enddo

c --- try to get the labels for simulation
      if (delsim) call fsfgetlabls

c --- print out summary of input
      if(fsfdbg.gt.0)then
        print*,'FSFINP: summary of TE input:'
        print*,'FSFINP:   # of TEs found in event      :',nid
        print*,'FSFINP:   # of unknown TEs (why ??)    :',n_unknwn
        print*,'FSFINP:   # of single VTX TEs found    :',n_singleVTX
        print*,'FSFINP:   # of multi-TE VTX TEs found  :',n_multiVTX
        print*,'FSFINP:   # of single FTD TEs found    :',n_singleFTD
        print*,'FSFINP:   # of multi-TE FTD TEs found  :',n_multiFTD
        do i = 1,detmax
          if (teindexused(detorder(i)).gt.0)then
            print*,'FSFINP:   # of TEs from ',detname(detorder(i)),
     1             ' is:',teindexused(detorder(i))
          endif
        enddo
      endif
      end


