*  
* $Id: fsfarrays.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
c --- This is where the TE data will be stored:
      real space(spacesize)
      integer ispace(spacesize)
      equivalence (ispace(1),space(1))
c --- Index of last used field in space/ispace
      integer spaceused,maxspaceused

c --- Now an Array containing the start indices for the TEs.
c     First index counts TEs, second one indicates the detector.
      integer teindex(temax,detmax),teindexused(detmax)
      integer tsindex(tsmax),tsindexused

c --- The corresponding TE ID number and size
      integer teid(temax,detmax), tesize(temax,detmax)
      logical teused(temax,detmax)

c --- The array where showertes are marked.
      integer  showerte(temax,detmax)

c --- the datadescriptor
      integer  datdescr(temax,detmax)

      common/fsfarrays/ tesize,teid,teindexused,teindex,
     1                  spaceused,space,teused,maxspaceused,
     2                  tsindex,tsindexused,datdescr,showerte

c --- I will copy the information from the TE's into the following
c     arrays. This will allow for easier access in the
c     fitroutines, as well for preprocessing some entries
c     (e.g. enlarge the errors to account for the shortcommings of
c      extrapolation procedures).

c --- First the coordinates:
      real x(temax,detmax),y(temax,detmax),z(temax,detmax)
      real R(temax,detmax),RPhi(temax,detmax)
      real theta(temax,detmax),phi(temax,detmax)
      real overrt(temax,detmax)

c --- Now the corresponding errors:
      real dx(temax,detmax),dy(temax,detmax)
      real dRPhi(temax,detmax),dz(temax,detmax)
      real dtheta(temax,detmax),dphi(temax,detmax)
      real doverrt(temax,detmax)


      common/fsftedata/x,y,R,Rphi,z,theta,phi,overrt,
     1                 dx,dy,dRphi,dz,dtheta,dphi,doverrt

      integer fsfnte(detmax),nsign(3),ndirect(2)
      integer ntspd(maxpd),detsum(0:maxpd,detmax)
      real    time(maxpd)
      common/fsfstat/fsfnte,ntspd,detsum,nsign,ndirect,time

