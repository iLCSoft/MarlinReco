*  
* $Id: fsfxerrs.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
c --- The following arrays shall contain cuts

      real xcx(detmax,detmax),xcy(detmax,detmax)
      real xcRPhiR(detmax,detmax),xcz(detmax,detmax)
      real xctheta(detmax,detmax),xcphi(detmax,detmax)
      real xcthest(detmax,detmax),xcphist(detmax,detmax)

c --- Extra error for phicut.
      real plx(detmax,detmax),ply(detmax,detmax)
      real plRPhiR(detmax,detmax),plz(detmax,detmax)
      real pltheta(detmax,detmax),plphi(detmax,detmax)
      real plthest(detmax,detmax),plphist(detmax,detmax)

      common /fsferrs/ xcx,xcy,xcRPhiR,xcz,
     1                 xctheta,xcphi,xcthest,xcphist,
     2                 plx,ply,plRPhiR,plz,
     3                 pltheta,plphi,plthest,plphist

c --- Error scaling for TEs in trackfit.
      real fiterrscale,tpcerrscale
      common/fsferrscale/fiterrscale,tpcerrscale

