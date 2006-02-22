*  
* $Id: fsfmkpair.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      logical function fsfmkpair(method,te1,det1,te2,det2)
*******************************************************************************
*    INPUT:
*      method   - This describes the way how the two TEs shall be
*                 proved to be a pair. The following methods are
*                 allowed:
*                 -1  - Use TE labels to decide, wether they belong together.
*                  0  - Define them NOT to be a pair
*                  1  - use helixextrapolation
*      te1,det1 - The TE and detector indizes of the starting detector.
*      te2,det2 - The TE and detector indizes of the second detector.
*
*    OUTPUT
*      Guess what :-)
*                                                                  Daniel Wicke
*******************************************************************************
      implicit none
#include "include/brcntl.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fsfdebug.inc"
#include "include/fsfpi.inc"
      integer method,te1,det1,te2,det2,jbit
      integer labli,lablj

      if (telabl(te1,det1).ne.0.and.
     1    telabl(te1,det1).eq.abs(telabl(te2,det2)).and. ! te1 must have +labl
     2    (te1.ne.te2.or.det1.ne.det2)) then
        labelsok = .TRUE.
      else
        labelsok = .FALSE.
      endif
c
      if(te1.eq.0.or.te2.eq.0)then
        fsfmkpair=.FALSE. !Index zero is no hit and thus wont give a pair
c
c ---   check for that theta of TE1 points to z of TE2 (if <> 0)
c       check only explicit these detectors
      else if((det2.eq.i_VTX.or.det2.eq.i_TPC.or.det2.eq.i_ITC).and.
     1        ((theta(te1,det1).lt.pi/2.and.z(te2,det2).lt.0).or.
     1         (theta(te1,det1).gt.pi/2.and.z(te2,det2).gt.0)))then
        fsfmkpair=.FALSE.
        if(fsfdbg.gt.1.and.z(te1,det1)*z(te2,det2).gt.0)then
          print*,'FSFMKPAIR: likesign z values, but theta1 opposite ?'
        endif
c
c --- protect shower TEs to extrapolate inwards
c
      else if(showerte(te1,det1).ne.0 .and.
     1        R(te1,det1).lt.R(te2,det2))then
        fsfmkpair=.FALSE.
        if(fsfdbg.gt.1)then
          print*,'FSFMKPAIR: Did not extrapolate shower TE of ',
     1            detname(det1),' to ',detname(det2)
        endif
c
c --- do not associate shower TEs
c
      else if(showerte(te2,det2).ne.0)then
        fsfmkpair=.FALSE.
        if(fsfdbg.gt.1)then
          print*,'FSFMKPAIR: Did not extrapolate to shower TE of ',
     1           detname(det2)
        endif
      else
        if(method.eq.0) then
          fsfmkpair=.FALSE.
        else if(method.eq.1) then
          call fsfhelix(fsfmkpair,te1,det1,te2,det2)
        else if(method.eq.-1) then
c ---     Check, if first te has enough parameters
c         measured to make an extrapolation.
c         Otherwise we might get trouble in the fit
          if(0.eq.jbit(datdescr(te1,det1),3) .or.
     1       0.eq.jbit(datdescr(te1,det1),4) .or.
     2       0.eq.jbit(datdescr(te1,det1),5) .or.
     3       0.eq.jbit(datdescr(te1,det1),8) .or.
     4       0.eq.jbit(datdescr(te1,det1),9) .or.
     5      (0.eq.jbit(datdescr(te1,det1),10)
     6        .and.0.eq.jbit(datdescr(te1,det1),11))
     7      )then
             fsfmkpair=.FALSE.
          else
             lablj=telabl(te1,det1)
             labli=telabl(te2,det2)
             fsfmkpair= abs(labli).eq.abs(lablj).and.lablj.ne.0
          endif
        else
          print*,'ERROR in FSFMKPAIR: Unknown method',method
          fsfmkpair=.FALSE.
        endif
      endif

      if(fsfdbg.gt.0)then
        if(.not.fsfmkpair.and.labelsok)then
          print*,'WARNING in FSFMKPAIR:',
     1           ' good TE combination (labels) not paired'
          print'(A,A5,A,I8,A,I8)',
     1           '   starting detector: ',detname(det1),
     2                           ' id # ',teid(te1,det1),
     3                          ' label ',telabl(te1,det1)
          print'(A,A5,A,I8,A,I8)',
     1           '     target detector: ',detname(det2),
     2                           ' id # ',teid(te2,det2),
     3                          ' label ',telabl(te2,det2)
        endif
      endif

      end

