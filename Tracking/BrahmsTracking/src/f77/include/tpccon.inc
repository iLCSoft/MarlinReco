*  
* $Id: tpccon.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:36  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
      REAL TVDRFT,myBFIELD,BCFGEV,TPRWDT,TSIGMD,TOMTAU,TPRFN2,TVOFFS
      REAL TXDRFT,TYDRFT,TVOFS0,TPCCLK
      INTEGER JPEDST
      COMMON /TPCCON/ TVDRFT(2),myBFIELD,BCFGEV,TPRWDT,TSIGMD,TOMTAU,
     1                TPRFN2(LTPDRO),
     2                JPEDST(LTSECT,2),TVOFFS(LTSECT),
     3                TXDRFT(2),TYDRFT(2),TVOFS0(LTSECT),TPCCLK

