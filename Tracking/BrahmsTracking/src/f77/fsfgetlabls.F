*  
* $Id: fsfgetlabls.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      SUBROUTINE FSFGETLABLS
*******************************************************************************
*
*      Routine to get labels for simulation to debug.
*
*      WRITTEN BY                                Daniel Wicke
*      CHANGES BY                                Kristian Harder, August 99
*
*******************************************************************************
      implicit none
C +CDE, UXCOM.
#include "include/brcntl.inc"
#include "include/fsfb.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfpi.inc"
#include "include/fsfsim.inc"
      integer te,det,nte,ierr,tknumb
      integer ilabdet(detmax),nolabdet(detmax),badlabdet(detmax),
     1        idmirror,nidmirror,nfudged,inegdet(detmax)
      real    trdat(11)

      if (.not.DELSIM) RETURN

      call vzero(ilabdet,detmax)
      call vzero(inegdet,detmax)
      call vzero(nolabdet,detmax)
      call vzero(badlabdet,detmax)
      idmirror = 0
      nidmirror= 0
      nfudged  = 0

      nte=tknumb('TE',0)

      do det=1,detmax
        do te=1,teindexused(det)
          if(teid(te,det).gt.nte)then
            write(6,*) 'FSFGETLABLS: fatal error - TE non-existent'
            stop
          endif
          CALL TKGLAB('TE',teid(te,det),telabl(te,det),IERR)
          IF (IERR.NE.0) THEN
            WRITE(6,*) 'FSFGETLABLS: TKGLAB FAILED.'
            STOP
          END IF
          IF (telabl(te,det).ne.0) THEN
c            CALL TKREAD('MCT',0,iabs(telabl(te,det)),trdat,IERR)
            IF (trdat(1).eq.0.0.and.trdat(2).eq.0.0) IERR=999
          ELSE
            IERR=999
          END IF
          IF (IERR.NE.0) THEN
            simtheta(te,det)=0.0
            simphi(te,det)=0.0
            simoverrt(te,det)=0.0
          ELSE
            simtheta(te,det)=atan2(sqrt(trdat(1)**2+trdat(2)**2),
     >                       trdat(3))
            simphi(te,det)=atan2(trdat(2),trdat(1))
            simoverrt(te,det)=bfieldz/sqrt(trdat(1)**2+trdat(2)**2)
          END IF
          if (telabl(te,det).ne.0.and.abs(telabl(te,det)).lt.10**7)then
            if (telabl(te,det).gt.0) then
              ilabdet(det)  = ilabdet(det)+1
            else
              inegdet(det)  = inegdet(det)+1
            endif
          else if (telabl(te,det).eq.0) then
            nolabdet(det) = nolabdet(det)+1
          else
            badlabdet(det)= badlabdet(det)+1
          endif
        enddo
      enddo
      if(fsfdbg.gt.0)then
        print*,'FSFGETLABLS: summary of label processing:'
        do det = 1,detmax
          if (ilabdet(detorder(det))+nolabdet(detorder(det))+
     1        badlabdet(detorder(det))+inegdet(detorder(det)).ne.0)then
            print'(A,A5,A,i5,A,i5,A,i5,A,i5)',
     1           ' FSFGETLABLS: ',detname(detorder(det)),
     2           ' labels ok:',ilabdet(detorder(det)),
     3           ' negative:',inegdet(detorder(det)),
     4           ' zero:',nolabdet(detorder(det)),
     5           ' bad:',badlabdet(detorder(det))
          endif
        enddo
      endif

      end


