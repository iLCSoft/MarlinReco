*  
* $Id: fsfshwr.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      function FSFSHWR(teid)
*******************************************************************************
*     Check if a given TE is used in a TS marked as hadronic shower TS.
*
*     INPUT: teid - Tanagra identifier of the TE to be checked.
*                                                                  Daniel Wicke
*******************************************************************************
      implicit none
      integer FSFSHWR
      integer teid,i
c +cde, EXTVCH.

      i=1
      FSFSHWR=0
c      do while(i.le.nntvch)
c         if(idtvch(i).eq.teid)then
c           call sbyt(mctvch(i),FSFSHWR,17,3)
c           i=nntvch+1 !**Quit Loop
c         else
c           i=i+1
c         endif
c      enddo

      end

