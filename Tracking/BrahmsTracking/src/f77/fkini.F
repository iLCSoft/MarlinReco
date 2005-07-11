*  
* $Id: fkini.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:18  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
************************************************************************
*                                                                      *
      SUBROUTINE FKINI
*                                                                      *
*   initialization of magnetic field                                   *
*   initialization of extrapolation surfaces from EXXSUR               *
*   printout of constants and tolerances used                          *
*   definitions of histograms                                          *
*                                                                      *
************************************************************************
*
#include "include/fksave.inc"
#include "include/fkdebug.inc"
#include "include/fkclun.inc"
#include "include/fkbdlo.inc"
#include "include/fkcons.inc"
#include "include/fkfild.inc"
#include "include/fkparm.inc"
#include "include/fkexts.inc"
#include "include/fkddes.inc"
*
*
* -ME- steer debugging and histograms
*
*      IDEB = 0
*      IHIS = 0
*
      WRITE(IPRNT,*) ' '
      WRITE(IPRNT,*) ' ************************************************'
      WRITE(IPRNT,*) ' **                                            **'
      WRITE(IPRNT,*) ' **   TRACK FIT IN BARREL AND FORWARD REGION   **'
      WRITE(IPRNT,*) ' **                                            **'
      WRITE(IPRNT,*) ' **      Constants and tolerances used :       **'
      WRITE(IPRNT,'(''  **      c*Bz ='',F10.6,T49,''**'')') CONSB
      WRITE(IPRNT,'(''  **      LHIST ='',I6,T49,''**'')') LHIST
      WRITE(IPRNT,'(''  **      Cut on chi2/ndeg ='',F6.2,T49,''**'')')
     &            CHI2MX
      WRITE(IPRNT,'(''  **      Cut on total chi2 prob. ='',F6.4,'//
     >            ' '' (fwd.)'',T49,''**'')')   PRTOT(0)
      WRITE(IPRNT,'(''  **      Cut on total chi2 prob. ='',F6.4,'//
     >            ' '' (brl.)'',T49,''**'')')   PRTOT(1)
      WRITE(IPRNT,'(''  **      Cut on smoo. chi2 prob. ='',F6.4,'//
     >            ' '' (fwd.)'',T49,''**'')')   PRSMO(0)
      WRITE(IPRNT,'(''  **      Cut on smoo. chi2 prob. ='',F6.4,'//
     >            ' '' (brl.)'',T49,''**'')')   PRSMO(1)
      WRITE(IPRNT,*) ' **                                            **'
      WRITE(IPRNT,*) ' ************************************************'
      WRITE(IPRNT,*) ' '
*
*   define histograms
*
c      IF(IHIS.NE.0) THEN
c        KHIST=LHIST
c        CALL HMDIR('//PAWC/TRACKING/DELFIT','S')
c        CALL HBOOK1(KHIST+11,'r of ref. point',100,0.,100.,0.)
c        CALL HBOOK1(KHIST+12,'Phi fit at ref. point',60,0.,360.,0.)
c        CALL HBOOK1(KHIST+13,'z fit at ref. point',100,-100.,100.,0.)
c        CALL HBOOK1(KHIST+14,'theta fit at ref. point',60,0.,180.,0.)
c        CALL HBOOK1(KHIST+15,'phi fit at ref. point',60,0.,360.,0.)
c        CALL HBOOK1(KHIST+16,'p fit at ref. point',100,-250.,250.,0.)
c        CALL HBOOK2(KHIST+17,'p vs theta',100,-250.,250.,30,0.,180.,0.)
c        CALL HBOOK2(KHIST+18,'p vs phi',100,-250.,250.,30,0.,360.,0.)
c        CALL HBOOK2(KHIST+19,'z vs theta',100,-100.,100.,30,0.,180.,0.)
c        CALL HBOOK1(KHIST+20,'impact parameter',100,-.1,.1,0.)
c        CALL HCDIR('//PAWC',' ')
*
c      endif
*
      END
