*  
* $Id: fxrefit.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fxrefit(tkr,newtkr,nte,teidlist,ier1,refitdebug)
      implicit none
#include "include/fkparm.inc"
      integer IDATTE(MXTE,NTEMX)
      real    DATTE(MXTE,NTEMX)
      equivalence (idatte,datte)
      integer tkr(MXTK),newtkr(MXTK),nte,teidlist(*),refitdebug
      integer nout,idout(NTEMX),ier1,ier2,iflte(NTEMX),ierx(NEXSMX)
      integer i,ierr

      do i=1,nte
        call tkread('TE',0,teidlist(i),datte(1,i),IDATTE(1,i),ierr)
        if(ierr.ne.0)then
          print*,'ERROR in FXREFIT: TKREAD failed with',ierr
        endif
      enddo

      ier1=0
CC      call fk3trk(3,nte,teidlist,idatte,tkr,
      call fk3trk(1000+3,nte,teidlist,idatte,tkr,
     1     newtkr,nout,idout,ier1,ier2,iflte,ierx)

      if(ier1.ne.0.and.refitdebug.gt.1)then
        print*,'FXREFIT: SubTK could not be fitted!!'
      endif
      end

