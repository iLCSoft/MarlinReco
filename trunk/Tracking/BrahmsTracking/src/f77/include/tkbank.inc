*  
* $Id: tkbank.inc,v 1.1.1.1 2005-07-11 15:15:48 aplin Exp $
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
* TRACKING SUBSYSTEM: PATTERN RECOGNITION AND TRACK MERGING
*
* TKBANK data structure
*
* Kristian Harder, 30 Sep, 1999
*
* many things taken from DELPHI TANAGRA manual by D. Bertrand, L.Pape
**************************************************************************

**************************************************************************
* a) MC track bank: true track parameters for all tracks in the detector
**************************************************************************
*
*    This bank is filled by the TK package before calling the local patrec
*    modules.
*
* --- total number of tracks that produced at least one hit

	INTEGER TKNTRK

* --- track information (equivalent arrays for integer and real access)
*     first index:
*       1 (real)    px
*       2 (real)    py
*       3 (real)    pz
*       4 (real)    E
*       5 (real)    x of track origin
*       6 (real)    y of track origin
*       7 (real)    z of track origin
*       8 (integer) GEANT particle code
*       9 (integer) GEANT track number
*      10 (integer) number of hits left in detector
*      11 (integer) HEPEVT track number (if any, zero otherwise)

***      REAL    TKMCTR(11,NTKMX)
***      INTEGER ITKMCT(11,NTKMX)
***      EQUIVALENCE (TKMCTR,ITKMCT)

**************************************************************************
* b) HIT bank: all hits collected in current event
**************************************************************************
*
*    All hits recorded in the current event are saved here for use
*    of the local pattern recognition modules. Filled by the TK package
*    before the local patrec modules are called.
*
*
* --- total number of hits

      INTEGER NTKHIT

* --- hit information (equivalent arrays for integer and real access)
*     first index:
*       1 (real)    x
*       2 (real)    y
*       3 (real)    z
*       4 (real)    energy
*       5 (integer) subdetector ID
*       6 (integer) track ID, 0 for real data (sic!)
*       7 (integer) pointer to first exclusion (0 if none)
*       8 (integer) number of exclusions
*       9 (integer) resolution code (see below)
*      10 (real)    resolution 1
*      11 (real)    resolution 2
*     position 7 and 8 are used for mutually exclusive hit candidates
*     (mirror hits etc.). they point to entries in the exclusion list
*     defined below.
*
*     resolution code: two bits are set. The lower value bit describes which
*                      resolution is described by resolution 1 (pos. 10), the
*                      higher value bit describes resolution 2 (pos.11).
*
*                      bit  0  - error on r*phi
*                      bit  1  - error on z
*                      bit  2  - error on r

****      REAL    RTKHIT(MXHT,NHMAX)
****      INTEGER ITKHIT(MXHT,NHMAX)
****      EQUIVALENCE (ITKHIT,RTKHIT)

* --- exclusion list (shared by all other banks)

      INTEGER NEXCL,IEXCL(MXEXCL),EXCLERR

* --- pointer array: IHPOINT(ID_DET) points to position of first hit
*     from subdetector ID_DET. IHNUMB(ID_DET) stores the number of hits
*     found in this detector

      INTEGER IHPOINT(3000),IHNUMB(3000)


**************************************************************************
* c) TE bank (subdetector output)
**************************************************************************
*
*
*    Following DELPHI nomenclature, the information delivered by a subdetector
*    is called track element (TE). A TE may contain fitted tracks (e.g. TPC),
*    but may also contain single hits (e.g. Presampler hits).
*    This bank must be filled by the local patrec modules.
*
* --- number of track elements (TEs).

      INTEGER NTE

* --- DELPHI TANAGRA TE structure for track elements/hits (equivalant arrays)
*             first index (m=length of packed covariance matrix):
*              1 (integer) subdetector ID
*              2 (integer) submodule: bit code for VTX/FTD layers used in TE
*              3 (integer) not used (yet) - MUST be set to zero
*              4 (integer) measurement code (see below)
*              5 (integer) pointer to end of TE (m+17)
*              6 (integer) charge (0=neutral, 1=positive, 2=negative, 3=unknown)
*              7 (integer) number of degrees of freedom
*              8 (real)    chi^2 of the fit
*              9 (real)    length of the track element
*             10 (real)    coordinate 1 of reference point: x or R
*             11 (real)    coordinate 2 of reference point: y or R*PHI
*             12 (real)    coordinate 3 of reference point: z
*             13 (real)    theta angle
*             14 (real)    phi angle
*             15 (real)    1/p at the reference point (in 1/(GeV/c))
*                          (or 1/p_t - see measurement code)
*             16 (real)    dE/dx (if measured - see measurement code)
*             17 (real)    covariance matrix for measured quantities
*             ...          (see below)
*             m+16 (real)  error on dE/dx if dE/dx is measured
*             m+17 (real)  must be 0.0
*
*             measurement code:
*             bit  0 = 0 for x,y,z coordinates
*                    = 1 for R,R*PHI,z coordinates
*             bit  1 = 0
*             bit  2 = 1 if coordinate 1 is known (index 10), 0 otherwise
*             bit  3 = 1 if coordinate 2 is known (index 11), 0 otherwise
*             bit  4 = 1 if coordinate 3 is known (index 12), 0 otherwise
*             ('known' means either defined by detector geometry or measured)
*             bit  5 = 1 if R*PHI (x) is measured for bit 0 = 1 (0)
*             bit  6 = 1 if z (y) is measured for bit 0 = 1 (0)
*             bit  7 = 1 if theta is measured
*             bit  8 = 1 if phi is measured
*             bit  9 = 1 if 1/p is measured
*             bit 10 = 1 if 1/p_t is measured
*             ('measured' really means measured, i.e. this quantity contributes
*              to the covariance matrix)
*             bit 11 = 0
*             bit 12 = 0
*             bit 13 = 1 if dE/dx is measured
*             bit 14-... = 0
*
*             covariance matrix:
*             every variable activated by a 1 in bit 5..10 is included
*             in the covariance matrix. for a 3d tracking detector like
*             a TPC, the measured quantities are in most cases
*
*             R*PHI and z (coordinates of the reference point at a reference
*                          radius R which therefore has no error itself)
*             theta,phi,1/p
*
*             In this example, the order of elements in the covariance matrix
*             is the following:
*
*             D(R*PHI,R*PHI),
*             C(R*PHI,z),     D(z,z),
*             C(R*PHI,theta), C(z,theta), D(theta,theta),
*             C(R*PHI,phi),   C(z,phi),   C(theta,phi),   D(phi,phi),
*             C(R*PHI,1/p),   C(z,1/p),   C(theta,1/p),   C(phi,1/p), D(1/p,1/p)

***      INTEGER ITE(MXTE,MXEVTE)
***      REAL    RTE(MXTE,MXEVTE)
***      EQUIVALENCE (ITE,RTE)

* --- additional information. first index:
*      1 (integer) position of first associated hit in hit list (see below)
*      2 (integer) number of hits
*      3 (integer) pointer to first exclusion (0 if none) in exclusion list
*      4 (integer) number of exclusions
*      5 (integer) track no. (positive if 95% of the hits belong to
*                             same track, negative if 75%, zero below.)
*     remark: the exclusion list is shared with the hit bank (see there).

***      INTEGER ITEDAT(5,MXEVTE)

* --- hit list: pointers to entries in the hit bank.

      INTEGER NHITTE,IHITTE(NHMAX)


**************************************************************************
* d) TS bank (ambiguous full track candidates)
**************************************************************************
*
*    Merging result as created by the TK package and DELSEARCH.
*    Ambiguous use of TEs is possible here; the unambiguous optimal
*    track reconstruction will be entered into the TK bank (see below)
*
* --- number of tracks segments (TS) found by DELSEARCH

      INTEGER NTS

* --- DELPHI TANAGRA TS structure
*             first index (m=length of packed covariance matrix):
*              1 (integer) module identifier (internal use only)
*              2 (integer) bit code for used detectors
*              3 (integer) measurement code (see below)
*              4 (integer) track type (see below)
*              5 (integer) number of TEs used in this TS
*              6 (integer) charge (0=neutral, 1=positive, 2=negative, 3=unknown)
*              7 (integer) not used (must be zero, though)
*              8 (integer) number of degrees of freedom of the track fit
*              9 (real)    chi^2 of the track fit
*             10 (real)    track segment length
*             11 (real)    x of TS start point
*             12 (real)    y of TS start point
*             13 (real)    z of TS start point
*             14 (real)    x of TS end point
*             15 (real)    y of TS end point
*             16 (real)    z of TS end point
*             17 (real)    coordinate 1 of reference point: x or R
*             18 (real)    coordinate 2 of reference point: y or R*PHI
*             19 (real)    coordinate 3 of reference point: z
*             20 (real)    theta angle
*             21 (real)    phi angle
*             22 (real)    1/p at the reference point (in 1/(GeV/c))
*                          (signed with geometric curvature of track
*                           measured from nside to outside of detector)
*             23 (real)    covariance matrix for measured quantities
*             ...          (as above)
*             37 (real)
*
*             measurement code:
*             bit  0 = 0 no estimate at all for the track segment parameters
*                    = 1 track segment parameters are given
*             bit  1 = 0 for x,y,z coordinates
*                    = 1 for R,R*PHI,z coordinates
*             bit  2 = 0 crude estimation of errors (diagonal cov. matrix)
*                    = 1 full covariance matrix
*             bit  3 = 0 all parameters are given if bit 0 is on
*                    = 1 no error matrix/starting point/end point
*
*             track type:
*             not yet documented. 1+8 (bits 0 and 3 set) means:
*             crude estimation, plane reference surface (cartesian coordinates),
*             no start/end points given

      INTEGER ITS(MXTS,NTSMX)
      REAL    RTS(MXTS,NTSMX)
      EQUIVALENCE (ITS,RTS)

* --- additional TS information. first index:
*      1 (integer) position of first associated TE in TE list (see below)
*      2 (integer) number of TEs
*      3 (integer) pointer to first exclusion (0 if none) in exclusion list
*      4 (integer) number of exclusions
*      5 (integer) track no. (positive if 95% of the hits belong to
*                             same track, negative if 75%, zero below.)
*     remark: the exclusion list is shared with the hit bank (see there).

      INTEGER ITSDAT(5,NTSMX)

* --- TE list. contains all TEs that are used by the TS. (referenced in ITSDAT)

      INTEGER NTSTEL,ITSTEL(MXEVTE)

**************************************************************************
* e) TK bank (the real tracking result - unambiguous track candidates)
**************************************************************************
*
*    Filled by TK package after resolving all TS ambiguities.
*
*
* --- number of tracks (TK)

      INTEGER NTK

* --- DELPHI TANAGRA TK structure
*             first index (m=length of packed covariance matrix):
*              1 (integer) module identifier (internal use only)
*              2 (integer) bit code for used detectors
*              3 (integer) measurement code (see below)
*              4 (integer) track type (as above)
*              5 (integer) number of TEs used in this TK
*              6 (integer) charge (0=neutral, 1=positive, 2=negative, 3=unknown)
*              7 (integer) not used (must be zero, though)
*              8 (integer) number of degrees of freedom of the track fit
*              9 (real)    chi^2 of the track fit
*             10 (real)    track length
*             11 (real)    x of TS start point
*             12 (real)    y of TS start point
*             13 (real)    z of TS start point
*             14 (real)    x of TS end point
*             15 (real)    y of TS end point
*             16 (real)    z of TS end point
*             17 (real)    coordinate 1 of reference point: x or R
*             18 (real)    coordinate 2 of reference point: y or R*PHI
*             19 (real)    coordinate 3 of reference point: z
*             20 (real)    theta angle
*             21 (real)    phi angle
*             22 (real)    1/p at the reference point (in 1/(GeV/c))
*                          (signed with geometric curvature of track
*                           measured from inside to outside of detector)
*             23 (real)    covariance matrix for measured quantities
*             ...          (as above)
*             37 (real)
*
*             measurement code:
*             bit  0 = 0 for x,y,z coordinates
*                    = 1 for R,R*PHI,z coordinates

      INTEGER ITK(MXTK,NTKMX)
      REAL    RTK(MXTK,NTKMX)
      EQUIVALENCE (ITK,RTK)

* --- additional TK information. first index:
*      1 (integer) position of first associated TE in TE list (see below)
*      2 (integer) number of TEs
*      3 (integer) track no. (positive if 95% of the hits belong to
*                             same track, negative if 75%, zero below.)

      INTEGER ITKDAT(3,NTKMX)

* --- TE list. contains all TEs that are used by the TK. (referenced in ITKDAT)

      INTEGER NTKTEL,ITKTEL(MXEVTE)



**************************************************************************

C --- the tracking common block
      COMMON /TKBANK/ TKNTRK,
*     &  TKMCTR,
     >                NTKHIT,
*     &  ITKHIT,
     >                NEXCL,IEXCL,EXCLERR,
     >                IHPOINT,IHNUMB,
     >                NTE,
*     ITE,ITEDAT,
     >                NHITTE,IHITTE,
     >                NTS,ITS,ITSDAT,NTSTEL,ITSTEL,
     >                NTK,ITK,ITKDAT,NTKTEL,ITKTEL
      SAVE   /TKBANK/

**************************************************************************

