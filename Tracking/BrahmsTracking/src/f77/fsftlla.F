*  
* $Id: fsftlla.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsftlla (n,list,label)
*******************************************************************************
*
*      get simple label for list of TEs:
*
*      INPUT: n        - number of TEs
*             list(*)  - list of identifiers
*      OUTPUT:label    - label (set to TE label if all TEs
*                               in the list come from same track,
*                               set to 0 otherwise)
*
* WRITTEN BY                                                    Markus Elsing
* CHANGES BY                                                    Kristian Harder
********************************************************************************
      IMPLICIT NONE

#include "include/fsfparams.inc"

      INTEGER n,list(*),label,label2,i,ierr

      if (n.lt.1) then
        write(6,*) 'ERROR IN FSFTLLA: no TE list entries'
        label=0
        return
      end if

      call tkglab('TE',list(1),label,ierr)
      if (ierr.ne.0) then
        write(6,*) 'ERROR IN FSTLLA: could not read label'
        label=0
        return
      end if

      do i=2,n
        call tkglab('TE',list(i),label2,ierr)
        if (ierr.ne.0) then
          write(6,*) 'ERROR IN FSTLLA: could not read label'
          label=0
          return
        end if
        if (label2.ne.label) label=0
      end do

      end


* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

* PATCH FSFALI:
*
*      FSALICHK    - steering routine for test hbooks
*
*      INIALIHIST  - books alignment histos to test detectors
*      FSCYLHIST   - books cylinder histos for RPHI,R,z detectors
*      FSPLHIST    - books plane histos for xyz detectors
*
*      ALICHECK    - fills histos to check alignment...
*      FSCYLFILL   - fills cylinder histos for tracks vs TE
*      FSPLFILL    - fills plane histos for tracks vs TE
*
*      FSGDAT      - get day of year from event date
*      FKPTKR      - fransform TKFB to tanagra structur
*      FKWTOC      - invert wiegth matrix of TE (even not complete)
*
* -----------------------------------------------------------------------------

* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

