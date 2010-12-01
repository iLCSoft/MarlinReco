*  
* $Id: fkcom.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:49  hvogt
* Brahms V308
*
*  
       integer  nte,idet,measur,ndter,ierte,iflag,nsurf,inumb,iords,
     $          ndforw,irffrw,irzfrw,ndfrw1,ndfrw2,
     $          ndback,irfbck,irzbck,ndbck1,ndbck2
       real*4   sinth,costh,pinv,chter,alrfte,
     $          chforw,chfrw1,chfrw2,
     $          chback,chbck1,chbck2
       real*4   alrphi,xrls,xeloss

       DOUBLE PRECISION PARREF,PARPLA,PARAM,paraft,DQTER,WGTER,
     & DQX,PARAMX,WGTX,COVX,PARAMS,DQFORW,WGFORW,DQFRW1,WGFRW1,
     & DQFRW2,WGFRW2,DQFRWE,WGFRWE,DQBACK,WGBACK,DQBCK1,WGBCK1,
     & DQBCK2,WGBCK2,DQBCKE,WGBCKE
*
       COMMON /FKSPXX/ SINTH,COSTH,PINV
     & ,NTE,IDET(NTEMX),MEASUR(NTEMX),NDTER(NTEMX),CHTER(NTEMX)
     & ,IERTE(NTEMX),IFLAG(NTEMX),ALRFTE(NTEMX)
     & ,NSURF,INUMB(NSRFMX),ALRPHI(NSRFMX),XRLS(NSRFMX),IORDS(NSRFMX)
     & ,XELOSS(NSRFMX)
     & ,NDFORW,CHFORW,IRFFRW(NTEMX),IRZFRW(NTEMX)
     & ,NDFRW1(NTEMX),CHFRW1(NTEMX),NDFRW2(NTEMX),CHFRW2(NTEMX)
     & ,NDBACK,CHBACK,IRFBCK(NTEMX),IRZBCK(NTEMX)
     & ,NDBCK1(NTEMX),CHBCK1(NTEMX),NDBCK2(NTEMX),CHBCK2(NTEMX)
*
       COMMON /FKDPXX/PARREF(6),PARPLA(6),PARAM(6,NTEMX),DQTER(5,NTEMX)
     & ,WGTER(15,NTEMX),DQX(5,NEXSMX),PARAMX(6,NEXSMX),WGTX(15,NEXSMX)
     & ,COVX(15,NEXSMX),PARAMS(6,NSRFMX),paraft(6,nsrfmx)
     & ,DQFORW(5),WGFORW(15),DQFRW1(5,NTEMX),WGFRW1(15,NTEMX)
     & ,DQFRW2(5,NTEMX),WGFRW2(15,NTEMX)
     & ,DQFRWE(5,NEXSMX),WGFRWE(15,NEXSMX)
     & ,DQBACK(5),WGBACK(15),DQBCK1(5,NTEMX),WGBCK1(15,NTEMX)
     & ,DQBCK2(5,NTEMX),WGBCK2(15,NTEMX)
     & ,DQBCKE(5,NEXSMX),WGBCKE(15,NEXSMX)
*
       COMMON /FKCOM/ TYPE(NTEMX),TYPSUR(NSRFMX)
       CHARACTER*4 TYPE,TYPSUR
*
