*  
* $Id: fkstkr.F,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:34  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:19  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:55  hvogt
* Brahms V308
*
*  
***********************************************************************
*                                                                     *
      SUBROUTINE FKSTKR(ISUREF,IDIR,IHIS,DAT,IERR)
*                                                                     *
*   computation of fitted track parameters at reference point, and    *
*   their covariance matrix                                           *
*                                                                     *
*   definition of reference point :  reference point in the first     *
*   or last detector seen, or intersection with an extrap. surface    *
*                                                                     *
*   Input  : ISUREF : serial number of reference surface              *
*            IDIR   : -1 for backward fit , 1 for forward fit         *
*            IHIS   : flag for histograms                             *
*                                                                     *
*   Output : DAT    : TKR data array (starting from word 17)          *
*            IERR   : error flag : 0 if OK                            *
*                                  2 if weight matrix not positive    *
*                                                                     *
***********************************************************************
*
#include "include/fksave.inc"
#include "include/fkcons.inc"
#include "include/fkparm.inc"
#include "include/fkcom.inc"
*
      DOUBLE PRECISION PAR,WGT,COV
      DIMENSION DAT(21),PAR(5),WGT(15),COV(15)
*
      KHIST=LHIST
*
*   fitted values of parameters and weight matrix at reference point
*
      INUM=INUMB(ISUREF)
*
*   reference point is a on a measurement surface (TE)
      IF(INUM.GT.0) THEN
        IF(IDIR.EQ.-1) THEN
          CALL FKDPVA(PARAMS(2,ISUREF),DQBCK2(1,INUM),PAR,5)
          CALL FKDPCO(WGBCK2(1,INUM),WGT,15)
        ELSE
          CALL FKDPVA(PARAMS(2,ISUREF),DQFRW2(1,INUM),PAR,5)
          CALL FKDPCO(WGFRW2(1,INUM),WGT,15)
        ENDIF
*
*   reference point is a on a extrapolation surface
      ELSE
        IF(IDIR.EQ.-1) THEN
          CALL FKDPVA(PARAMS(2,ISUREF),DQBCKE(1,-INUM),PAR,5)
          CALL FKDPCO(WGBCKE(1,-INUM),WGT,15)
        ELSE
          CALL FKDPVA(PARAMS(2,ISUREF),DQFRWE(1,-INUM),PAR,5)
          CALL FKDPCO(WGFRWE(1,-INUM),WGT,15)
        ENDIF
      ENDIF
      IF(ABS(PAR(3)).LT..001.OR.ABS(PAR(5)).LT.0.0000001)THEN
        IERR = 1
        RETURN
      ENDIF
*
*   invert weight matrix at reference point
*
      CALL FKMI5D(WGT,COV,IERR)
      IF(IERR.NE.0) RETURN
*
*   transform parameters and covariance matrix to TANAGRA parameters
*
      CALL FKSPAR(TYPSUR(ISUREF),PARAMS(1,ISUREF),PAR,COV,DAT)
*
*   histograms
*
ccc      IF(IHIS.NE.0) THEN
ccc        PMOM=1./DAT(6)
ccc        CALL HCDIR('//PAWC/TRACKING/DELFIT',' ')
ccc        CALL HFILL(KHIST+14,RADDEG*DAT(4),0,1.)
ccc        CALL HFILL(KHIST+15,RADDEG*DAT(5),0,1.)
ccc        CALL HFILL(KHIST+16,PMOM,0,1.)
ccc        CALL HFILL(KHIST+17,PMOM,RADDEG*DAT(4),1.)
ccc        CALL HFILL(KHIST+18,PMOM,RADDEG*DAT(5),1.)
ccc        IF(TYPSUR(ISUREF).EQ.'CYLI') THEN
ccc          CALL HFILL(KHIST+11,DAT(1),0,1.)
ccc          CALL HFILL(KHIST+12,RADDEG*DAT(2)/DAT(1),0,1.)
ccc          CALL HFILL(KHIST+13,DAT(3),0,1.)
ccc          CALL HFILL(KHIST+19,DAT(3),RADDEG*DAT(4),1.)
ccc*   compute impact parameter (approximation valid for small imp. par.)
ccc          SALPHA=.5*DAT(1)*PAR(5)
ccc          IF(ABS(SALPHA).LE.1.) THEN
ccc            ALPHA=ASIN(SALPHA)
ccc            PARIMP=DAT(1)*(PAR(4)-ALPHA)
ccc            CALL HFILL(KHIST+20,PARIMP,0.,1.)
ccc          ENDIF
ccc        ENDIF
ccc      ENDIF
*
      END
