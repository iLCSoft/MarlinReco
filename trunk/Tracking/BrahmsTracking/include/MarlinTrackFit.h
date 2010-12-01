#ifndef MARLINTRACKFIT_H
#define MARLINTRACKFIT_H 1
#include "HelixClass.h"

/**
 * Utility class which performs track fit using Kalman filter
 * Utility is based on the DELPHI code which performs track fit
 * Author : A. Raspereza, MPI Munich
 */

extern "C" {
  void tfithl_(int & NPT, double * XF,double * YF, float * RF, float *PF, double * WF,float * ZF,
	       float * WZF, int & IOPT, float * VV0,
	       float * EE0, float & CH2PH, float & CH2Z);
}

extern "C" {
  void trackfit_(int & nhits, int * idet, int * itype, 
		 float * x, float * y, float * z, float * phireso, float * zreso,
		 float * ref, int & ierr, float * rfit, float * rfite, float & chi2, int & ndf,
		 int & noutl, int * idoutl, int & fitcode);
  
}


class MarlinTrackFit {
 public:

  /**
   *  Constructor. 
   */ 
  MarlinTrackFit();

  /**
   *  Destructor. 
   */
  ~MarlinTrackFit();

  /**
   * Method performs track fitting taking into account energy loss and MS 
   * First simple helix fit is performed to define initial track parameters
   * at the PCA to primary IP. If simple helix fit converges and has qood quality  
   * set by variable chi2PrefitCut
   * then DELPHI fitting routine trkfit.F is envoked
   * Inputs :
   *          useExtraPoint : 1 - an additional extra point (P.C.A) is used in fit with
   *                               relatively large errors to improve D0 and Z0
   *                               determination @ P.C.A
   *                           0 - number of points in fit equals to the number of hits in track
   *          optFit : 0 - FORTRAN routine tfithl is used for 
   *                       prefit to determine initial track parameters
   *                   1 - ClusterShapes class is used for prefit
   *          nhits - number of tracker hits used in fit
   *          bField - magnetic field in Tesla
   *          idet - list of the hit detector identifiers 
   *          itype - list of the hit types ( 2 - planar disks, 
   *                                          3 - cyllindrical or laddered detector )
   *          chi2PrefitCut - cut on the chi2 of prefit. If chi2 of prefit > chi2
   *                          then prefit is regarded to fail and routine returns -1
   *          x,y,z - list of the hit Cartesian coordinates
   *          RPReso - list of the hit resolutions in R-Phi projection
   *          ZReso  - list of the hit resolutions along Z
   * Outpus :
   *          param  - list of the track parameters at the PCA
   *          param[0] - Omega (signed curvuture)
   *                     sign "+" corresponds to clock-wise direction
   *          param[1] - tan(Lambda) (tangence of the dip angle, lambda = pi/2 - Theta)
   *          param[2] - Azimuthal angle of the particle momentum @ PCA
   *          param[3] - D0 (IP in the R-Phi projection @ PCA )
   *          param[4] - Z0 (z coordinate displacement @ PCA )
   *          eparam[15] - covariance matrix of the track parameters
   *                       (lower-left part of symmetric matrix)
   *          RefPoint[3] - reference point of track fit
   *          chi2 - fit chi2
   *          ndf - number of degrees of freedom
   *          chi2rphi - prefit chi2 in the R-Phi (X-Y) plane
   *          chi2z    - prefit chi2 in the S-Z plane
   *          lhits - list of hit status after fit, 
   *                  lhits[i] = 0 means that hit is thrown away after fitting procedure
   *                  lhits[i] = 1 means that hit is kept in the fitting procedure
   *          returned integer indicate error flag
   *          error = 0 - no error, fit is successfull
   *          error = -1 - simple helix fit failed
   *          error = 1 - Delphi fit failed, but the simple helix fit is successfull 
   * 
   * Subroutine converts covariance matrix from the 
   * TANAGRA (DELPHI code) representation <br>
   * (R-Phi,R-Phi) <br>
   * (R-Phi,ZR)     (ZR,ZR) <br>
   * (R-Phi,Theta)  (ZR,Theta)  (Theta,Theta) <br>
   * (R-Phi,PhiR)   (ZR,PhiR)   (Theta,PhiR)   (PhiR,PhiR) <br>
   * (R-Phi,1/p)    (ZR,1/p)    (Theta,1/p)    (PhiR,1/p)  (1/p,1/p) <br>
   * into canonical representation <br>
   * (D0,D0) <br>
   * (Phi0,D0)   (Phi0,Phi0) <br>
   * (Omega,D0)  (Omega,Phi0)  (Omega,Omega) <br>
   * (Z0,D0)     (Z0,Phi0)     (Z0,Omega)    (Z0,Z0) <br>
   * (TanL,D0)   (TanL,Phi0)   (TanL,Omega)  (TanL,Z0)  (TanL,TanL) <br>
   *
   * TANAGRA   Vector Eta = (R-Phi,z,Theta,phi,1/p) <br>
   *            where R-Phi : R*Phi @ fixed Radius of the reference point <br>
   *                     ZR : z coordinate of the reference point <br>
   *                  Theta : polar angle of the momentum @ reference point <br>
   *                   PhiR : azimuthal angle of the momentum @ reference point <br>
   *                   1/p  : inverse of momentum magnitude (signed quantity) <br>
   *                          sign is "+" corresponds to anticlock-wise direction <br>
   *                          as opposed to definition of Omega (in canonical format) <br>
   * Canonical Vector X = (Omega,TanLambda,Phi,D0,Z0) <br>
   *
   * Covariance matrix is calculated according to <br>
   * UX[i,j] = VEta[l,m]*dXdEta[i,l]*dXdEta[j,m] <br>
   * where UX[i,j] is the cov matrix in the canonical representation <br>
   * Veta[l,m] is the cov matrix in the TANAGRA (DELPHI) representation <br>
   * and dXdEta[k,n] - is Jacobian dX_k/dEta_n <br>
   * <br>
   */  
  int DoFitting(int useExtraPoint, int fitOpt, // inputs
		int & nhits, float bField, int * idet, int * itype, // inputs 
		float chi2PrefitCut, // inputs 
		float * x, float * y, float * z, float * RPReso, float * ZReso, // inputs 
		float * param, float * eparam, float * RefPoint, float & chi2, int & ndf, // outputs
		float & chi2rphi, float & chi2z,int * lhits); // outputs

  /**
   * Method fills covariance matrix according to 
   * apriori specified resolutions 
   * Covariance matrix is assumed to be diagonal
   * Inputs :  reso[4]   - array of resolutions : <br>
   *                0    - momentum resolution delta(Pt)/Pt^2 <br>
   *                1    - angular resolution (radians) <br>
   *                2    - D0  resolution (mm) <br>
   *                3    - Z0  resolution (mm) <br>
   *           param[5]  - track parameters in the canonical representation <br>
   * Output : eparam[15] - covariance matrix <br>
   */ 
  void CrudeErrorEstimates(float BField, float * reso, float * param, float * eparam);
  void setParametersForIPErrors(float * par);

 private:    

  float _parIpReso[3];
    
  void CalculateChi2(HelixClass & helix, int & nhits, 
		     float * xhit, float * yhit, float * zhit,
		     int * iType, float * rphireso, float * zreso,
		     float * chi2hit, float & chi2rphi, float & chi2z);

  void ConvertLCtoTANAGRA(float * param, float & bField, float * ref);
  void ConvertTANAGRAtoLC(float * rfit, float & bField, int & fitCode, float * param);
  void ConvertCovMatrix(float * rfit, float * rfite, float * param, float & bField, int fitCode, float * eparam);

};

#endif
