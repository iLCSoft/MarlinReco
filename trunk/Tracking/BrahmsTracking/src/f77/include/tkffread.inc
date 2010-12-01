*  
* $Id: tkffread.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
**************************************************************************
*
* steering flags from FFREAD cards
*
**************************************************************************

c     IPATRC:
c     steering of what is returned as TE from each subdetector.
c     for each subdetector (order: see TURNON sequence):
c       do local pattern recognition?
c       1=yes, realistic
c       2=yes, ideal
c
c     TKSTDBG: debug level for TKSTEER (0,1,2)
c     DELADBG: debug level for DELAMBI (0,1,2)
c     TKSTHST: fill tracking performance histograms
c     ITKPERF: switch on (1) or off (0, default) performance analysis
c
c     quality cuts a track has to pass in order to get counted in performance
c     analysis:
c     PFMINH  minimal number of hits (a value of at least 3 makes some sense)
c     PFPCUT  minimal momentum (3d)
c     PFCOST  maximal abs(cos(theta))
c     PFKINE  true to count only 'physics tracks', i.e. tracks that were read
c             from kinfile
c
c     SIMDET-like output steering:
c

      INTEGER IPATRC(3),TKSTDBG,TKSTHST,DELADBG
      INTEGER PFMINH,ITKPERF
      INTEGER IFBEST,IFHIST
      REAL    PFPCUT,PFCOST
      LOGICAL PFKINE
      COMMON /TKFFREAD/ IPATRC,TKSTDBG,TKSTHST,DELADBG,ITKPERF,
     >                  PFMINH,PFPCUT,PFCOST,PFKINE,
     >                  IFBEST,IFHIST
      SAVE   /TKFFREAD/


