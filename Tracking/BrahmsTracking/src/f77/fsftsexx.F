*  
* $Id: fsftsexx.F,v 1.1 2005-11-03 15:16:14 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2003/05/23 13:17:53  hvogt
* Brahms V308
*
*  
      subroutine fsftsexx(its,n_surf,t_surf,r_surf,z_min,z_max,
     1                    textra,ierrsf,ier)
*******************************************************************************
*
*     Extrapolate internal TS structure to radius x.
*
*      INPUT:  its       - TS index in space array
*              n_surf    - number of surfaces
*              t_surf(*) - type of extr. surf (CYLI/PLAN)
*              r_surf(*) - array of radii of cylinders
*              z_min(*)  - minimal z
*              z_max(*)  - maximal z
*      OUTPUT: textra(30,*) - extrapolations at surfaces
*              ierrsf(*)    - error flag per surface
*              ier          - EXXTRK error flag
*                                                                  Daniel Wicke
*******************************************************************************
      IMPLICIT NONE

c +CDE, EXXPRM.
#include "include/exxgeo.inc"
#include "include/fsfpi.inc"
#include "include/fsfparams.inc"
#include "include/fsfarrays.inc"
#include "include/fkexts.inc"
#include "include/fktkx.inc"
#include "include/brcntl.inc"

      INTEGER i,its,n_surf,ierrsf(*),ier,jbit
      REAL    r_surf(*),z_min(*),z_max(*),textra(30,*)
      REAL    brsurf(nexsmx),bzmin(nexsmx),bzmax(nexsmx)
      CHARACTER*4 t_surf(*),btypx(nexsmx)
      INTEGER ibnexs

      REAL    tstart(mxte)
      INTEGER itstrt(mxte)
      EQUIVALENCE (tstart,itstrt)
      INTEGER nout, idout(ntemx), iflte(ntemx)
      REAL    dtk(mxtk)


c --- set exx parameters
      thxmin    =0.
      thxmax    =pi
      phymin    =0.
      phymax    =2.*pi

c --- reset output
      call vzero(textra,30*n_surf)
      call vzero(ierrsf,n_surf)
      ier    = 0

c --- check range
      if (n_surf.gt.nexsmx) then
        write(6,*) 'FSFTSEXX: too many extrapolation surfaces'
        stop
      end if

c --- save original extrapolation surfaces
      ibnexs=nexs
      do i=1,nexs
        btypx(i)=typx(i)
        brsurf(i)=rzsurf(i)
        bzmin(i)=zrmin(i)
        bzmax(i)=zrmax(i)
      end do

c --- copy surfaces
      nexs=n_surf
      do i = 1,n_surf
        typx(i)   = t_surf(i)
        rzsurf(i) = r_surf(i)
        zrmin(i)  = z_min(i)
        zrmax(i)  = z_max(i)
      enddo

c --- now set up tstart as pseudo-TE
      call vzero(itstrt,mxte)
      itstrt(1)=0
      itstrt(2)=0
      itstrt(3)=0
      itstrt(4)=jbit(ispace(its+ 3),2)+4+8+16+32+64+128+256+512
      itstrt(5)=32
      tstart(6)=space(its+ 6)
      call ucopy(space(its+ 8),tstart(7),3)
      call ucopy(space(its+17),tstart(10),6)
      tstart(16)=0.0
      call ucopy(space(its+23),tstart(17),15)
      tstart(32)=0.0

c --- extrapolate with error propagation
      call tk2fit(1,tstart,0,dtk,ier,nout,idout,iflte,ierrsf)

c --- convert extrapolation results into TANAGRA TiX format
c     (but do not fill words 1 to 6 - they are not used in BRAHMS)
      do i=1,nexs
        call ucopy(datx(1,i),textra(7,i),24)
      end do

c --- restore original extrapolation surfaces
      nexs=ibnexs
      do i=1,nexs
        rzsurf(i)=brsurf(i)
        zrmin(i)=bzmin(i)
        zrmax(i)=bzmax(i)
      end do

c --- check output
      if (ier.ne.0) then
        if (fsfdbg.gt.0)
     1    print*,'WARNING in FSFTSEXX: TK2FIT failed'
      else
        if (fsfdbg.gt.1) then
          do i = 1,n_surf
            print*,'FSFTSEXX: to surface #',i,' error (0=ok):',ierrsf(i)
          enddo
        endif
      endif

      end

