*  
* $Id: fsfprintte.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfprintte(idat,dat,ndat)
*******************************************************************************
*
*     Print contents of a TE.
*     idat and dat *MUST* be equivalent arrays to produce useful results!
*
*      INPUT:  idat - equivalenced array
*              dat  - same as real array
*              ndat - length of array
*                                                                  Daniel Wicke
*******************************************************************************
      implicit none
#include "include/fsfclun.inc"
#include "include/detnums.inc"
      real dat(*)
      integer ndat,idat(*),jbit,textptr,j
      character*20 textbuf

         if(idat(1).eq.ID_TPC)then
            textbuf='TPC'
         else if(idat(1).eq.ID_ITC)then
            textbuf='ITC'
         else if(idat(1).eq.ID_SIT+1)then
            textbuf='SIT1'
         else if(idat(1).eq.ID_SIT+2)then
            textbuf='SIT2'
         else if(idat(1).eq.ID_FTD)then
            textbuf='FTD multipoint'
         else if(idat(1).gt.ID_FTD.and.idat(1).lt.ID_FTD+100)then
            textbuf='FTD hit'
         else if(idat(1).eq.ID_VTX)then
            textbuf='CCD VTX multipoint'
         else if(idat(1).gt.ID_VTX.and.idat(1).lt.ID_VTX+100)then
            textbuf='CCD VTX hit'
         else if(idat(1).eq.ID_VTXP)then
            textbuf='APS VTX multipoint'
         else if(idat(1).gt.ID_VTXP.and.idat(1).lt.ID_VTXP+100)then
            textbuf='APS VTX hit'
         else if(idat(1).eq.ID_FCH)then
            textbuf='FCH'
         else if(idat(1).gt.3000)then
            textbuf='TS as dummy TE'
         else
            textbuf='unknown'
         endif
         print*,' ',idat(1),' ',textbuf
         print*,' ',idat(2)
         print*,' ',idat(3)

         if(0.eq.jbit(idat(4),1))then
            textbuf='x,y,z'
            textptr=6
         else
            textbuf='R,RPhi,z'
            textptr=9
         endif
         if(1.eq.jbit(idat(4),2))then
            textbuf(textptr:)='(non-std)'
            textptr=textptr+9
         else
            textbuf(textptr:)='(std)'
            textptr=textptr+5
         endif
         print*,' ',idat(4),' ',textbuf,jbit(idat(4),1),
     1          jbit(idat(4),2)
         if(0.eq.jbit(idat(4),3))then
            print*,' First coordinate *NOT* measured.'
         endif
         if(0.eq.jbit(idat(4),4))then
            print*,' Second coordinate *NOT* measured.'
         endif
         if(0.eq.jbit(idat(4),5))then
            print*,' Third coordinate *NOT* measured.'
         endif
         if(0.eq.jbit(idat(4),6))then
            print*,' *NO* error for t1 (x or RPHi).'
         endif
         if(0.eq.jbit(idat(4),6))then
            print*,' *NO* error for t2 (y or z).'
         endif
         if(0.eq.jbit(idat(4),8))then
            print*,' Theta direction *NOT* measured.'
         endif
         if(0.eq.jbit(idat(4),9))then
            print*,' Phi direction *NOT* measured.'
         endif
         if(1.eq.jbit(idat(4),10))then
            print*,' idat(15) is 1/p.'
         endif
         if(1.eq.jbit(idat(4),11))then
            print*,' idat(15) is 1/p_T.'
         endif
         if(1.eq.jbit(idat(4),12))then
            print*,' idat(15) is E.'
         endif
         if(1.eq.jbit(idat(4),13))then
            print*,' idat(16) is beta.'
         endif
         if(1.eq.jbit(idat(4),14))then
            print*,' idat(16) is dE/dx.'
         endif
         if(1.eq.jbit(idat(4),15))then
            print*,' idat(16) is time.'
         endif

         if(idat(6).eq.3)then
            textbuf=' unknown charge'
         else if(idat(6).eq.2)then
            textbuf=' negative charge'
         else if(idat(6).eq.1)then
            textbuf=' positive charge'
         else if(idat(6).eq.0)then
            textbuf=' neutral particle'
         endif
         print*,' ',idat(6),textbuf
         print*,' ',idat(7),' Degrees of freedom'
         do j=8,16
            print*,' dat(',j,')=',dat(j)
         enddo
         print*,' no error matrix printout'
      END

