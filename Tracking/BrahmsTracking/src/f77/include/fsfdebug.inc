*  
* $Id: fsfdebug.inc,v 1.1 2005-07-21 16:37:00 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:49  hvogt
* Brahms V308
*
*  
c --- simulated track parameters
      real    simtheta(temax,detmax),simoverrt(temax,detmax)
      real    simphi(temax,detmax)
      common /fsfsimdat/ simtheta,simoverrt,simphi
c --- for the label of the TEs
      integer telabl(temax,detmax)
      common /fsftelabl/ telabl
c --- for the label of the current combination tested
      logical labelsok
      integer label
      common /culabel/labelsok,label

