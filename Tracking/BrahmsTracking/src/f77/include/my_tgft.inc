*  
* $Id: my_tgft.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
      INTEGER NTGFT
      INTEGER MY_JTGFTC,MY_JTGFTG,MY_LTGFTA,MY_JTGFDF
      REAL MY_JTGFIR,MY_JTGFTL,MY_JTGFP0
      REAL MY_JTGFD0,MY_JTGFZ0,MY_JTGFEM
      REAL MY_JTGFCD
      COMMON/MY_TGFT/NTGFT,
     & MY_JTGFIR(MAXCND),MY_JTGFTL(MAXCND),MY_JTGFP0(MAXCND),
     & MY_JTGFD0(MAXCND),MY_JTGFZ0(MAXCND),MY_JTGFEM(MAXCND,15),
     & MY_JTGFCD(MAXCND),MY_JTGFDF(MAXCND),MY_JTGFTC(MAXCND),
     & MY_JTGFTG(MAXCND),MY_LTGFTA(MAXCND)

