*  
* $Id: my_t1ft.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
      INTEGER NT1FT
      INTEGER MY_JT1FTC,MY_JT1FTG,MY_LT1FTA,MY_JT1FDF
      REAL MY_JT1FIR,MY_JT1FTL,MY_JT1FP0
      REAL MY_JT1FD0,MY_JT1FZ0,MY_JT1FEM
      REAL MY_JT1FCD
      COMMON/MY_T1FT/NT1FT,
     & MY_JT1FIR(MAXCND),MY_JT1FTL(MAXCND),MY_JT1FP0(MAXCND),
     & MY_JT1FD0(MAXCND),MY_JT1FZ0(MAXCND),MY_JT1FEM(MAXCND,15),
     & MY_JT1FCD(MAXCND),MY_JT1FDF(MAXCND),MY_JT1FTC(MAXCND),
     & MY_JT1FTG(MAXCND),MY_LT1FTA(MAXCND)

