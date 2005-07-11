*  
* $Id: chgtrak.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
*  
* $Log: not supported by cvs2svn $
* Revision 1.1.1.1  2005/06/18 10:05:35  gaede
* Initial local version
*
* Revision 1.1  2005/06/13 16:14:29  aplin
* *** empty log message ***
*
* Revision 1.1.1.1  2003/05/23 13:17:47  hvogt
* Brahms V308
*
*  
      INTEGER NMAXCHT,NMXCHIT
      PARAMETER(NMAXCHT=10000)
      PARAMETER(NMXCHIT=400)
      INTEGER NCHGTRK,ICHGTRK(NMAXCHT),NCHGHIT(NMAXCHT),NTKHITS(NMAXCHT)
      INTEGER ICHTRN(NMAXCHT)
      INTEGER IFITCODE(NMAXCHT)
*      REAL  CGTKHIT(NMAXCHT,NMXCHIT,7)
*      REAL CHTTMOM(NMAXCHT,4), CHTTVTX(NMAXCHT,3), CHTTQ(NMAXCHT)
*      REAL CHTRPAR(NMAXCHT,5), CHTRCOV(NMAXCHT,15), CHTRCHI(nmaxcht)

*      COMMON /CHGTRAK/NCHGTRK,ICHGTRK,NCHGHIT,CGTKHIT,NTKHITS
      COMMON /CHGTRAK/NCHGTRK,ICHGTRK,NCHGHIT,NTKHITS
C(KH) ---end

      SAVE /CHGTRAK/
**      SAVE /CHTTRAK/
*      SAVE /CHRTRAK/

**********************************************************************
*
*     common chgtrak
*
*     this set of common block is used to store
*     the hits generated and reconstructed and the
*     results from the track fit in the different
*     tracking detectors
*
*     The set contains three common blocks:
*
*     CHGTRAK   track parameters taken from GEANT
*               NCHGTRAK(ntrk)   number of tree level tracks
*               ICHGTRK(ntrk)    internal array linking to the
*                             number of GEANT tracks
*               NCHGHIT(ntrk)    number of hits on track
*               CHTKHIT(ntrk,nhit,5)    generated hits
*                    index 1   track number
*                    index 2   hit number
*                    index 3  1 x
*                             2 y
*                             3 z
*                             4 deposited charge, if appropriate, for hit
*                             5 code word for hit from DETNUMS
*                 NTKHITS
*
*     CHTTRAK:  true track parameters
*               CHTTMOM(ntrk,4) px,py,pz,E of track
*               CHTTVXT(ntrk,3) origin of track (x,y,z)
*               CHTTQ(ntrk)     charge of track
*
*     CHRTRAK   reconstructed track parameters
*               IFITCODE(ntrk) return code if the track fit
*               ICHTRN(ntrk)  number of hits used in the track fit
*               CHTRPAR(ntrk,5) parameter of the track after the fit:
*                               theta, phi, 1./p,
*               CHTRCOV(ntrk,15) upper diagonal elements of the
*                               covariance matrix
*               CHTRCHI(ntrk)   chi**2 of the track fit
*
*     Last update     V.L.Morgunov     04-Jule-2002                         *
**********************************************************************
