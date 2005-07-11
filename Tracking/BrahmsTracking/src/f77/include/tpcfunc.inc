*  
* $Id: tpcfunc.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:36  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:49  hvogt
* Brahms V308
*
*  
*     BRTP1 bank (LBRT1)
*     (1)  number of TPC hits (not used)
*     n*(7 words) TPC hit information

*      TPCHITS(ITPCSET,ITPCOFF) =
*     &  Q( LBRTP1 + NTPC*(ITPCOFF-1)+ITPCSET+1 )

      TPCHITS(ATTRIBUTE,AHIT) = STFREADTPCHITS(ATTRIBUTE,AHIT)	
*     BRTP2 bank: hit link bank
*      HITLNK(ITPCOFF) = IQ(LBRTP2 + ITPCOFF)
* FIXME: Replaced zebra bank with array for now with maximum size NTPCOMX=50000. 
*        Should probably be done as stl vector in MARLIN

