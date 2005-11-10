*  
* $Id: tk2fit.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:52  hvogt
* Brahms V308
*
*  
      SUBROUTINE TK2FIT(NTE,DTE,IOUTR,DTK,ierror,nout,idou,iflte,ierx)
C****************************************************************************
c* Version number: 0.0      Date of last modification  09.07.99             *
c*------------------------------------------------------------------------- *
c* List of modifications since last version:                                *
c*   KH   inhibit fkref calls with nte<=2 that sometimes caused crashes     *
c*------------------------------------------------------------------------- *
C* Main Authors:                                                            *
c* Name          Institute                        Contact email             *
C* K. Moenig     DESY-Zeuthen                     moenig@ifh.de             *
c*                                                                          *
c****************************************************************************
c input:
c   NTE   = number of TEs
C   DTE   = track elements to be fitted
c   ioutr = flag to handle outliers/extrapolations
c Output:
C  dtk    = fit result in tk format
c  ierror = error flag
c  nout   = number of outliers
c  idout  = identfier of outliers
c  iflte  = flag for each input TER , Not 0 means TER dropped from fit
c  ierx   = flag for each extrapolation surface
C ---------------------------------------------------------------------------
      implicit none
      save
#include "include/coildims.inc"
#include "include/fkparm.inc"
c#include "include/gconst.inc"
      integer jbyt, jbit
      integer NTE,IOUTR,ierror
      real    DTE(mxte,NTE),DTK(*)
      integer i, imeas
      integer idte(ntemx),idou(*),iflte(*)
      integer ier1,ier2, nout, ierx(*)
      real ref(mxts)
      integer irf(mxts)
      equivalence (ref,irf)
      real perrmn, terrmn
      logical lref
      integer imts
      parameter (imts = 1+8)
      integer imtsi
      parameter (imtsi = 1+2+4)

*
* if one TE contains already full information use as reference
*

      call vzero(ref,mxts)
      call vzero(idte,ntemx)
      if (nte.ge.3) then
         perrmn = 0.1**2
         terrmn = 0.001**2
      else
         perrmn = 1.e10
         terrmn = 1.e10
      endif
      lref = .false.
      do i=1,nte
         idte(i) = i
         imeas = jbyt(dte(4,i),6,6)
         if (imeas.eq.31 .or. imeas.eq.47) then
            if (dte(22,i).lt.terrmn.and.dte(31,i).lt.perrmn) then
* take this as reference
*
               lref  =.true.
               irf(1)=0
               irf(2)=0
               irf(3)=imts
               if (jbit(dte(4,i),1).eq.1) call sbit1(irf(3),2)
               call ucopy(dte(10,i),ref(17),22)
               if (jbit(dte(4,i),11).eq.1) ref(22)=sin(ref(20))/ref(22)
               perrmn = dte(31,i)
            endif
         endif
      enddo
*
*  if no complete TE get reference from a simple polar inversion
*
      if (.not.lref) then
         if (nte.ge.3) then
           call fkref(nte,dte,ref,irf)
         else
           ierror=-1
           return
         endif
      endif
      ierror = 0
c      call fk3trk(ioutr+1000,nte,idte,dte,ref,dtk,
      call fk3trk(ioutr,nte,idte,dte,ref,dtk,
     $            nout,idou,ier1,ier2,iflte,ierx)
      if (ier1.ne.0 .and.lref) then
         if (nte.ge.3) then
           call fkref(nte,dte,ref,irf)
         else
           ierror=-2
           return
         endif
         call fk3trk(ioutr,nte,idte,dte,ref,dtk,
     $               nout,idou,ier1,ier2,iflte,ierx)
      endif
      if (ier1.ne.0.or.ier2.ne.0) then

*
         ierror=ier1+ier2*10
c         print *,' ===> FK3TRK fails with ierror:',ierror
c        print*, ' => first fit',nout, abs(1./dtk(22)), sqrt(dtk(37)),
c     $          ' was:',abs(1./ref(22)), sqrt(ref(37))
         return
      endif
*
* ME copy output and iterate fit for testing
*
      irf(1) = 0
      irf(2) = 0
      irf(3) = imtsi
      if (jbit(dtk(3),1).eq.0) call sbit0(irf(3),2)
      call ucopy(dtk(4),ref(4),34)
      call fk3trk(10+ioutr,nte,idte,dte,ref,dtk,
     $            nout,idou,ier1,ier2,iflte,ierx)
      if (ier1.ne.0.or.ier2.ne.0) then
         ierror=ier1+ier2*10
c         print *,' ===> FK3TRK fails with ierror:',ierror
         return
      else
c        print*, ' => second fit',nout, abs(1./dtk(22)), sqrt(dtk(37)),
c     $          ' was:',abs(1./ref(22)), sqrt(ref(37))
      endif
*


      end
