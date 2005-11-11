*  
* $Id: fsfplcut.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsfplcut(helixpar,idxv,accepted,pd)
*******************************************************************************
*
*      Loop over TEs of combination, extrapolate and cut on dphi for
*      association. If all are within cut range, accept combination.
*
*      INPUT:  helixpar - parameters of helix from polar inversion
*              idxv     - list of TEs
*              pd       - primary detector, triggers cut value
*      OUTPUT: accepted - .TRUE. for good combination
*                                                                 Markus Elsing
*******************************************************************************
      IMPLICIT NONE
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfpi.inc"
#include "include/fsfxerrs.inc"

      real    helixpar(6),paramf(6)
      integer idxv(detmax),pd
      logical accepted
      integer i,ierr,ij

      accepted = .FALSE.

      do i=1,detmax
        if(idxv(i).ne.0)then
c --- extrapolate
          if(helixpar(4).eq.0.or.helixpar(6).eq.0)then
            ierr=-1
          else
            CALL fsfxts (helixpar,paramf,idxv(i),i,ierr)
          endif
          if(ierr.ne.0)then
            if(fsfdbg.gt.0)then
              print*,'WARNING in FSFPLCUT: failed for ',detname(pd),
     1               ' ->',detname(i)
            endif
            accepted=.FALSE.
            return
          else

#if defined(HBOOK)
            if (fsfhist.gt.0) then
c --- do some HISTOs here !
              ij=10000*pd+50*i
c --- labels ok ?
              if(labelsok.and.
     >           telabl(idxv(pd),pd).eq.abs(telabl(idxv(i),i)))then
                ij = ij + 30
              else
                ij = ij + 40
              endif
              call fsfhst (paramf,idxv(pd),pd,idxv(i),i,ij)
            endif
#endif
c --- now call test routine to compare extrapolation and TEs
            CALL fsfcomp (paramf,idxv(pd),pd,idxv(i),i,
     1                    plx,ply,plRPhiR,plz,
     2                    pltheta,plphi,plthest,plphist,
     3                    accepted)
c --- return, if failed !
            if (.not. accepted) then
              if(fsfdbg.gt.1.or.
     1          (fsfdbg.gt.0.and.labelsok))then
                print'(A,A5,A,I8,A,I8)',
     1               '   target detector off was: ',detname(i),
     2               ' id # ',teid(idxv(i),i),
     3               ' label ',telabl(idxv(i),i)
              endif
              return
            endif
          endif
        endif
      enddo
      end


