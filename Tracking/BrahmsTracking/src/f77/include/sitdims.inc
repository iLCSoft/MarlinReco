*  
* $Id: sitdims.inc,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:48  hvogt
* Brahms V308
*
*  
      integer nsit
      parameter (nsit=2)
      REAL SITRAD(nsit), SITTHK, SITLEN(nsit), SITRPRES, SITZRES
      COMMON/SITDIMS/SITRAD, SITTHK, SITLEN, SITRPRES, SITZRES

