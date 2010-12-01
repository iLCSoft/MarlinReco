*  
* $Id: brcntl.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.2  2003/12/15 14:57:08  hvogt
* lcio I/O included
*
* Revision 1.1.1.1  2003/05/23 13:17:47  hvogt
* Brahms V308
*
*  
*     this common block contains variables which are
*     used to steer the main BRAHMS program, and its
*     dependencies. These contain input/output,
*     and, for some packages, histogramming
*
      INTEGER IFHITF, IHITF(20), IFWDST, IDSTF(20), IFFORM

      INTEGER LFITYP, NPFILE, NPAFIL, NAMELEN
      INTEGER IKIUN, BRSKIP, IKINSKIP
      INTEGER IECASEEN(10)
      INTEGER BRIOVERSION, BRBKGVERSION, BRWRVERSION
      real kinpcut
      LOGICAL BREVSKIP
      CHARACTER*10 KI_STATUS, KI_FORM
      CHARACTER*80 KI_FILE
      PARAMETER (NAMELEN=80,
     &  NPFILE=5+NAMELEN/4)
      COMMON / BRCNTL / LFITYP, IKIUN, BRSKIP, IECASEEN, IKINSKIP,
     &  BREVSKIP, IFHITF, IHITF, IFWDST, IDSTF, IFFORM,
     &  BRIOVERSION, BRBKGVERSION,
     &  BRWRVERSION,
     &  kinpcut
      COMMON / BRCNTLC / KI_STATUS, KI_FORM, KI_FILE
      SAVE /BRCNTL/
      SAVE /BRCNTLC/

*     logical to control the histogramming in the VTX patrec
      LOGICAL VXPHIS_HIST
*      PARAMETER (VXPHIS_HIST = .FALSE.)
      COMMON /VXPHIS/ VXPHIS_HIST
      SAVE /VXPHIS/

*     integer variable to control the track merging histogramming
      integer FSFDBG   ! 0,1,2
      integer FSFHIST  ! 0,1
      COMMON /FSFDBF/ FSFDBG,FSFHIST
      SAVE /FSFDBF/

#if defined(BRLCIO)
*     ponters needed for LCIO I/O in Brahms
#include "include/cbrlcio.inc"
#endif

