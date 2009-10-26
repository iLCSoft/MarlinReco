#include "MarlinTrackFit.h"
#include "ClusterShapes.h"
#include <math.h>

MarlinTrackFit::MarlinTrackFit() {
    _parIpReso[0] = 0.0;
    _parIpReso[1] = 0.0;
    _parIpReso[2] = 1.0;
}

MarlinTrackFit::~MarlinTrackFit() {}


int MarlinTrackFit::DoFitting(int useExtraPoint, int fitOpt, // inputs
			      int & nhits, float bField, int * idet, int * itype, // inputs 
			      float chi2PrefitCut, // inputs 
			      float * x, float * y, float * z, float * RPReso, float * ZReso, // inputs 
			      float * param, float * eparam, float * RefPoint, float & chi2, int & ndf,
			      float & chi2rphi, float & chi2z, int * lhits) { // outputs

  // Create and fill up some intermediate arrays
  // needed to perform fits -->

  float * xhit = new float[nhits+1];
  float * yhit = new float[nhits+1];
  float * zhit = new float[nhits+1];
  int * iDet = new int[nhits+1];
  int * iTyp = new int[nhits+1];
  float * rphireso = new float[nhits+1];
  float * zreso = new float[nhits+1];
  float * phi = new float[nhits];
  double * xd = new double[nhits];
  double * yd = new double[nhits];
  float * r = new float[nhits];
  double * wfr = new double[nhits];
  float  * wfz = new float[nhits];
  float * ampl = new float[nhits];
  int * hitIndex = new int[nhits];
  int * idoutl = new int[nhits];
  int * idRep = new int[nhits];
  float * chi2hit = new float[nhits];
  float ref[6];
  float rfit[6];
  float rfite[15];

  float dzMin = 1.0e+20;
  float dzMax = -1.0e+20;
  float zBegin,zEnd;
  for (int i=0;i<nhits;++i) {
    // conversion from mm to cm
    // Delphi fit works with units of cm
    lhits[i] = 1;
    hitIndex[i] = i;
    xhit[i] = 0.1*x[i];
    yhit[i] = 0.1*y[i];
    zhit[i] = 0.1*z[i];
    ampl[i] = 1.0;
    xd[i] = xhit[i];
    yd[i] = yhit[i];
    phi[i] = atan2(yhit[i],xhit[i]);
    r[i] = sqrt(xhit[i]*xhit[i]+yhit[i]*yhit[i]);
    wfz[i] = 1./(0.01*ZReso[i]*ZReso[i]);
    wfr[i] = 1./(0.01*RPReso[i]*RPReso[i]);
    rphireso[i] = 0.01*RPReso[i]*RPReso[i];
    zreso[i] = 0.01*ZReso[i]*ZReso[i];    
    iDet[i] = idet[i];
    iTyp[i] = itype[i];
    float rhit = xhit[i]*xhit[i] +  yhit[i]*yhit[i];
    if (rhit>dzMax) {
      dzMax = rhit;
      zEnd = zhit[i];
    }
    if (rhit<dzMin) {
      dzMin = rhit;
      zBegin = zhit[i];
    }
  }
  
  HelixClass  helix;

  // do prefit using simple helix fit model ---->
  float chi2Sh;
  if (fitOpt == 0 || fitOpt >= 3) { // simple helix fit with subroutine tfithl
    int IOPT = 2;
    // fit with simple helix hypothesis using tfithl
    tfithl_(nhits,xd,yd,r,phi,wfr,zhit,
	    wfz,IOPT,param,eparam,chi2rphi,chi2z);    
    param[3] = param[3]*param[0]/fabs(param[0]);
    ConvertLCtoTANAGRA(param,bField,ref);
  }
  else if (fitOpt == 1) { // simple helix fit with ClusterShapes class
    // create instance of class ClusterShapes
    ClusterShapes * shapes = new ClusterShapes(nhits,ampl,xhit,yhit,zhit);
    float parSh[5];
    float dparSh[5];
    float distmax = 0;
    // do fitting
    shapes->FitHelix(500, 0, 1, parSh, dparSh, chi2Sh, distmax);
    float x0Sh = parSh[0];
    float y0Sh = parSh[1];
    float r0Sh = parSh[2];
    float bzSh = parSh[3];
    float phi0Sh = parSh[4];
    float signPz = 1;
    delete shapes;
    if (zEnd<zBegin)
      signPz = -1;
    helix.Initialize_BZ(x0Sh, y0Sh, r0Sh, 
			bzSh, phi0Sh, bField,signPz,
			zBegin);
    param[0] = helix.getOmega();
    param[1] = helix.getTanLambda();
    param[2] = helix.getPhi0();
    param[3] = helix.getD0();
    param[4] = helix.getZ0();
    ConvertLCtoTANAGRA(param,bField,ref);
  }
  else if (fitOpt == 2) { // Initial parameters are set by user
      param[0] = 10.*param[0];
      param[3] = 0.1*param[3];
      param[4] = 0.1*param[4];
      ConvertLCtoTANAGRA(param,bField,ref);
  }

  float omegaTmp = param[0];
  float tanLambdaTmp = param[1];
  float phi0Tmp = param[2];
  float d0Tmp = param[3];
  float z0Tmp = param[4];
  helix.Initialize_Canonical(phi0Tmp,d0Tmp,z0Tmp,omegaTmp,tanLambdaTmp,bField);
  CalculateChi2(helix, nhits, 
		xhit, yhit, zhit,
		iTyp, rphireso, zreso,
		chi2hit, chi2rphi, chi2z);  
   
  chi2Sh = chi2rphi + chi2z;
  int ndfSh = 2*nhits - 5;

  if (chi2Sh <=0 || ndfSh <=0) {
      std::cout << "MarlinTrackFit --> WARNING !" << std::endl;
      std::cout << "Chi2Sh=" << chi2Sh << std::endl;
      std::cout << "NdfSh=" << ndfSh << std::endl;
  }
  
  if (chi2Sh/float(ndfSh) > chi2PrefitCut) { // prefit condition not fulfilled
    delete[] xhit;
    delete[] yhit;
    delete[] zhit;
    delete[] rphireso;
    delete[] zreso;
    delete[] phi;
    delete[] xd;
    delete[] yd;
    delete[] r;
    delete[] wfr;
    delete[] wfz;
    delete[] iDet;
    delete[] iTyp;
    delete[] ampl;
    delete[] chi2hit;
    delete[] hitIndex;
    delete[] idoutl;
    delete[] idRep;
    chi2 = chi2Sh;
    ndf = ndfSh;
    // conversion from cm to mm
    param[0] = 0.1*param[0];
    param[3] = 10.0*param[3];
    param[4] = 10.0*param[4];
    return -1;
  }


  int nfit = nhits;

  float d0 = param[3];
  float phi0 = param[2];
  float z0 = param[4];

  RefPoint[0] = -d0*sin(phi0);
  RefPoint[1] = d0*cos(phi0);
  RefPoint[2] = z0;

  // Extra point @ PCA with large resolutions
  xhit[nhits] = RefPoint[0];
  yhit[nhits] = RefPoint[1];
  zhit[nhits] = RefPoint[2];
  rphireso[nhits] = 100.0;
  zreso[nhits] = 100.0;
  iDet[nhits] = 2;
  iTyp[nhits] = 2;  
  nfit = nhits ;
  // If extra point is used, effective number of
  // fit points is incremented 
  if (useExtraPoint > 0) 
      nfit = nhits+1 ;

  int ierr;
  int fitCode;
  int noutl = 0;

  if (fitOpt >= 3 )
      ref[0] = 0;  

  trackfit_(nfit,iDet,iTyp,xhit,yhit,zhit,rphireso,zreso,ref,ierr,
	    rfit,rfite,chi2,ndf,noutl,idoutl,fitCode);  
  
  if (ierr!=0) noutl = 0;

  int nh = nhits;
  if (noutl<nhits) 
    nh = noutl;

  for (int ioutl=0;ioutl<nh;++ioutl) {
    int id = idoutl[ioutl]-1;
    if (id>-1 && id<nhits)
      lhits[id] = 0;
  }
    

  // In case when fit fails, try to remove hits with large contributions to chi2  
  if ((ierr != 0) || (chi2 < 0) || (chi2/float(ndf) > 100.) || fabs(rfit[5]) < 1e-10 ) {
//       std::cout << "!!!!!!!!! Fit failure !!!!!!!!!" << std::endl;
//       std::cout << "Error = " << ierr 
// 		<< "  chi2 = " << chi2 
// 		<< "  NDF  = " << ndf
// 		<< "  fabs(rfit[5]) = " << fabs(rfit[5]) << std::endl;
    int NPT = 0;
    float chi2New = 0;
    for (int ihit=0;ihit<nhits;++ihit) {
      if (chi2hit[ihit] < 100.) {
	xhit[NPT] = xhit[ihit];
	yhit[NPT] = yhit[ihit];
	zhit[NPT] = zhit[ihit];
	iTyp[NPT] = iTyp[ihit];
	iDet[NPT] = iDet[ihit];
	rphireso[NPT] = rphireso[ihit];
	zreso[NPT] = zreso[ihit];
	chi2hit[NPT] = chi2hit[ihit];
	chi2New += chi2hit[NPT];
	idRep[NPT] = ihit;
	NPT++;
      }
      else {
	lhits[ihit] = 0;
      }
    }                
    nfit = NPT;
    if (nfit>=3) {
      trackfit_(nfit,iDet,iTyp,xhit,yhit,zhit,rphireso,zreso,ref,ierr,
		rfit,rfite,chi2,ndf,noutl,idoutl,fitCode);  
      if (ierr!=0) noutl = 0;
      int nh = nfit;
      if (noutl<nfit) 
	nh = noutl;

      for (int ioutl=0;ioutl<nh;++ioutl) {
	int id = idoutl[ioutl]-1;
	if (id>-1 && id<nfit) {
	  int idTrue = idRep[id];
	  if (idTrue>-1 && idTrue<nhits)
	    lhits[idTrue] = 0;	  
	}
      }      
      if ((ierr != 0) || (chi2 < 0) || (chi2/float(ndf) > 100.) || fabs(rfit[5]) < 1e-10 ) {
	ndfSh = nfit;
	chi2Sh = chi2New;
      }
      else {
	ierr = 1;
// 	std::cout << "!!!!!!! Repeated fit successfull !!! " << std::endl;
// 	std::cout << "Error = " << ierr 
// 	      << "  chi2 = " << chi2 
// 	      << "  NDF  = " << ndf
// 	      << "  fabs(rfit[5]) = " << fabs(rfit[5]) << std::endl;	
      }
    }
    else {
       for (int ihit=0;ihit<nhits;++ihit) {
	 lhits[ihit] = 1;
       }
    }
  }  

  // If fit didn't converge use the fitted quantities obtained with 
  // simple helix fit
  if ((ierr > 1) || (chi2 < 0) || (chi2/float(ndf) > 100.) || fabs(rfit[5]) < 1e-10 ) {
    chi2 = chi2Sh;
    ndf  = ndfSh;
    // conversion from cm to mm
    param[0] = 0.1*param[0];
    param[3] = 10.0*param[3];
    param[4] = 10.0*param[4];
    // Error estimates from tfithl and ClusterShapes are unreliable 
    // Hence we use some apriori values 
    // sigma(Omega)/Omega = sigma(pT)/pT = 1e-3*pT
    // angular resolution = 1mrad
    // D0 resolution = 1mm (as for TPC)
    // Z0 resolution = 1mm (as for TPC)
    float reso[4];
    reso[0] = 1e-3;
    reso[1] = 0.001;
    reso[2] = 1.0;
    reso[3] = 5.0;
    CrudeErrorEstimates(bField,reso,param,eparam);

//     std::cout << "!!!!!!! Repeated fit failure !!!!!!! " << std::endl;
//     std::cout << "Error = " << ierr 
// 	      << "  chi2 = " << chi2 
// 	      << "  NDF  = " << ndf
// 	      << "  fabs(rfit[5]) = " << fabs(rfit[5]) << std::endl;	
//     std::cout << "D0 Error = " << sqrt(eparam[0])
// 	      << "  Z0 Error = " << sqrt(eparam[9]) 
// 	      << "  Omg Error = " << sqrt(eparam[5]) 
// 	      << "  Phi Error = " << sqrt(eparam[2]) 
// 	      << "  TL Error = " << sqrt(eparam[14]) << std::endl;
    

    delete[] xhit;
    delete[] yhit;
    delete[] zhit;
    delete[] rphireso;
    delete[] zreso;
    delete[] phi;
    delete[] xd;
    delete[] yd;
    delete[] r;
    delete[] wfr;
    delete[] wfz;
    delete[] iDet;
    delete[] iTyp;
    delete[] ampl;    
    delete[] chi2hit;
    delete[] hitIndex;
    delete[] idoutl;
    delete[] idRep;
    if (ierr == 5)
	ierr = -1;
    return ierr;
  }

  ConvertTANAGRAtoLC(rfit,bField,fitCode,param);
  ConvertCovMatrix(rfit,rfite,param,bField,fitCode,eparam);

  // convert to mm
  param[0] = 0.1*param[0];
  param[3] = 10.*param[3];
  param[4] = 10.*param[4];

  d0 = param[3];
  phi0 = param[2];
  z0 = param[4];

  RefPoint[0] = -d0*sin(phi0);
  RefPoint[1] = d0*cos(phi0);
  RefPoint[2] = z0;
  

  if (noutl>50) {
    std::cout << "MarlinTrackFit ---> Too many outliers : " << noutl << std::endl;
  }

  // Check if the d0 or z0 errors have too low values
  float omega = param[0];
  float tanlambda = param[1];
  phi0 = param[2];
  d0 = param[3];
  z0 = param[4];
  helix.Initialize_Canonical(phi0, d0, z0, omega, tanlambda, bField);
  float momX = helix.getMomentum()[0];
  float momY = helix.getMomentum()[1];
  float momZ = helix.getMomentum()[2];
  float totMom = sqrt(momX*momX+momY*momY+momZ*momZ); 
  float TransMom = sqrt(momX*momX+momY*momY);
  float sinTheta = TransMom/totMom;
  float Parameter = totMom*pow(sinTheta,1.5);
  float MinIpReso = _parIpReso[0]+_parIpReso[1]/pow(Parameter,_parIpReso[2]);
  float ErrorD0 = sqrt(eparam[0]);
  float ErrorZ0 = sqrt(eparam[9]);

  //  SJA: This is far too dangerous. Removed. 26/10/09
//  if (ErrorD0<MinIpReso || ErrorZ0<MinIpReso) {
////       std::cout << "Par = " << _parIpReso[0]
//// 		<< "  " << _parIpReso[1]
//// 		<< "  " << _parIpReso[2]
//// 		<< "  " << _parIpReso[3]
//// 		<< " ; Min = " 
//// 		<< " ED0 = " << 1e+3*ErrorD0
//// 		<< " EZ0 = " << 1e+3*ErrorZ0
//// 		<< std::endl;
//      float errPhi = eparam[2];
//      float errOmg = eparam[5];
//      float errTL  = eparam[14];
//      for (int ip=0;ip<15;++ip) 
//	  eparam[ip] = 0;
//      float Error = MinIpReso + 0.002;
//      eparam[0] = Error*Error;
//      eparam[2] = errPhi;
//      eparam[5] = errOmg;
//      Error = MinIpReso + 0.002;
//      eparam[9] = Error*Error;
//      eparam[14]= errTL;
//  }



  delete[] xhit;
  delete[] yhit;
  delete[] zhit;
  delete[] rphireso;
  delete[] zreso;
  delete[] phi;
  delete[] xd;
  delete[] yd;
  delete[] r;
  delete[] wfr;
  delete[] wfz;
  delete[] iDet;
  delete[] iTyp;
  delete[] ampl;
  delete[] chi2hit;
  delete[] hitIndex;
  delete[] idoutl;
  return ierr;


}

void MarlinTrackFit::CrudeErrorEstimates(float BField, float * reso, float * param, float * eparam) {

  for (int iPar=0;iPar<15;++iPar) {
    eparam[iPar] = 0.0;
  }

  float omega = param[0];
  float tanlambda = param[1];
  float phi0 = param[2];
  float d0 = param[3];
  float z0 = param[4];
  float XTL = 1+tanlambda*tanlambda;
  HelixClass helix;
  helix.Initialize_Canonical(phi0, d0, z0, omega, tanlambda, BField);
  float pT = helix.getPXY();

  // (Omega,Omega)
  // delta(Omega) = Omega*delta(pT)/pT
  // reso[0] = delta(pT)/pT^2
  float sOmega = reso[0]*omega*pT;
  eparam[5] = sOmega*sOmega; 

  // (TanLambda,TanLambda)
  // TanLambda = cot(Theta) -> delta(Tanlambda) = delta(Theta)*(1+TanLambda^2)
  float sTL = reso[1]*XTL;
  eparam[14] = sTL*sTL;

  // (Phi0,Phi0)
  // delta(Phi0) = delta(Theta)/sin(Theta) = delta(Theta)*Sqrt(1+TanLambda^2)
  float sPhi = reso[1]*XTL;
  eparam[2] = sPhi*sPhi;

  // (D0,D0)
  eparam[0] = reso[2]*reso[2];

  // (Z0,Z0)
  eparam[9] = reso[3]*reso[3];


}

void MarlinTrackFit::CalculateChi2(HelixClass & helix, int & nhits, 
				   float * xhit, float * yhit, float * zhit,
				   int * iTyp, float * rphireso, float * zreso,
				   float * chi2hit, float & chi2rphi, float & chi2z) {
    float phi0L = helix.getPhi0();
    float z0L = helix.getZ0();
    float d0L = helix.getD0();
    chi2rphi = 0.0;
    chi2z = 0.0;
    float refp[3];
    refp[0] = -d0L*sin(phi0L);
    refp[1] = d0L*cos(phi0L);
    refp[2] = z0L;
    for (int i=0;i<nhits;++i) {
      if (iTyp[i] == 2) {
	float point[3];
	helix.getPointInZ(zhit[i],refp,point);
	float dX = xhit[i] - point[0];
	float dY = yhit[i] - point[1];
	float dR2 = dX*dX + dY*dY;
	chi2rphi = chi2rphi + dR2/rphireso[i];	
// 	if (dR2/rphireso[i]>1000.) rphireso[i] = 1e+20;
	chi2hit[i] = dR2/rphireso[i];
      }
      else {
	float point[6];
	float Radius = sqrt(xhit[i]*xhit[i]+yhit[i]*yhit[i]);
	helix.getPointOnCircle(Radius, refp, point);
	float dX1 = xhit[i] - point[0];
	float dY1 = yhit[i] - point[1];
	float dZ1 = zhit[i] - point[2];
	float dX2 = xhit[i] - point[3];
	float dY2 = yhit[i] - point[4];
	float dZ2 = zhit[i] - point[5];
	float dist1 = dX1*dX1+dY1*dY1+dZ1*dZ1;
	float dist2 = dX2*dX2+dY2*dY2+dZ2*dZ2;
	float dX,dY,dZ;
	if (dist1<dist2) {
	  dX = dX1;
	  dY = dY1;
	  dZ = dZ1;
	}
	else {
	  dX = dX2;
	  dY = dY2;
	  dZ = dZ2;
	} 
	float dRSquare = dX*dX + dY*dY;
	float dZSquare = dZ*dZ;
	chi2rphi = chi2rphi + dRSquare/rphireso[i];
	chi2z = chi2z + dZSquare/zreso[i];
//  	if (dRSquare/rphireso[i] > 1000. || dZSquare/zreso[i] > 1000.) {
//  	  rphireso[i] = 1e+20;
//  	  zreso[i] = 1e+20;
//  	}
	chi2hit[i] = dRSquare/rphireso[i] + dZSquare/zreso[i];
      }
    }

}

void MarlinTrackFit::ConvertLCtoTANAGRA(float * param, float & bField, float * ref) {

    float tanlambda = param[1];
    float omega = param[0];  
    float phi0 = param[2];
    float z0 = param[4];
    float d0 = param[3];
    float xInitial =  -d0*sin(phi0);
    float yInitial = d0*cos(phi0);
      
    // convert canonical helix parameters in the 
    // vector of the reference point in TANAGRA format
    ref[0] = xInitial;
    ref[1] = yInitial;
//   float RHO = sqrt(xInitial*xInitial+yInitial*yInitial);
//   float PHI = atan2(yInitial,xInitial);
//   if (PHI<0) PHI += 2*acos(-1.0);
//   ref[0] = RHO;
//   ref[1] = RHO*PHI;
    ref[2] = z0;
    ref[4] = phi0;
    ref[3] = 0.5*acos(-1.) - atan(tanlambda);
    ref[5] = -100.0*omega/(0.299792458*bField);
    ref[5] = ref[5]*fabs(sin(ref[3]));
    
}

void MarlinTrackFit::ConvertTANAGRAtoLC(float * rfit, float & bField, int & fitCode, float * param) {

    // Calculation of helix parameters 
    // in canonical form
    // measurement code (R,R-PHI,Z)
    float xx = rfit[0]*cos(rfit[1]/rfit[0]);
    float yy = rfit[0]*sin(rfit[1]/rfit[0]);
    float pos[3];
    // conversion to mm
    // measurement code (R,R-PHI,Z)
    pos[0] = 10*xx;
    pos[1] = 10*yy;
    pos[2] = 10*rfit[2];
    if (fitCode == 0) { // measurement code (X,Y,Z)
      pos[0] = 10*rfit[0];
      pos[1] = 10*rfit[1];
    }
    float mom[3];
    mom[0] = sin(rfit[3])*cos(rfit[4])/fabs(rfit[5]);
    mom[1] = sin(rfit[3])*sin(rfit[4])/fabs(rfit[5]);
    mom[2] = cos(rfit[3])/fabs(rfit[5]);
    float charge = -rfit[5]/fabs(rfit[5]);
    
    HelixClass helix;
    helix.Initialize_VP(pos,mom,charge,bField);
    
    // conversion to cm
    param[0] = 10.*helix.getOmega();
    param[1] = helix.getTanLambda();
    param[2] = helix.getPhi0();
    param[3] = 0.1*helix.getD0();
    param[4] = 0.1*helix.getZ0();    
    
}


void MarlinTrackFit::ConvertCovMatrix(float * rfit, float * rfite, float * param, float & bField, int fitCode, float * eparam) {
 // Calculation of cov matrix in canonical form
  // Subroutine converts covariance matrix from the 
  // TANAGRA format 
  // (R-Phi,R-Phi)
  // (R-Phi,ZR)     (ZR,ZR)
  // (R-Phi,Theta)  (ZR,Theta)  (Theta,Theta)
  // (R-Phi,phi)    (ZR,phi)    (Theta,phi)    (phi,phi)
  // (R-Phi,1/p)    (ZR,1/p)    (Theta,1/p)    (phi,1/p)  (1/p,1/p)
  // into canonical format
  // (Omega,Omega)
  // (Omega,TanLambda) (TanLambda,TanLambda)
  // (Omega,Phi0)      (TanLambda,Phi0)       (Phi0,Phi0)
  // (Omega,D0)        (TanLambda,D0)         (Phi0,D0)    (D0,D0)
  // (Omega,Z0)        (TanLambda,Z0)         (Phi0,Z0)    (D0,Z0)  (Z0,Z0)
  // TANAGRA   Vector Eta = (R-Phi,ZR,Theta,phi,1/p)
  // Canonical Vector X = (Omega,TanLambda,Phi,D0,Z0)
  // UX[i,j] = VEta[l,m]*dXdEta[i,l]*dXdEta[j,m]
  //
  // Output of the DELPHI fit routine
  // rfit[0] = fixed radius of the reference point (rho, eta0)
  // rfit[1] = R*Phi at the reference point 
  // rfit[2] = Theta angle of the momentum vector @ the reference point
  // rfit[3] = phi angle @ the reference point
  // declaration of cov matricies and Jacobian
  double VEta[5][5];
  double UX[5][5];
  double dXdEta[5][5];

   // Fill array of initial covariance matrix
  int counter = 0;
  for (int i=0;i<5;++i) {
    for (int j=0;j<=i;++j) {
      VEta[i][j] = rfite[counter];
      VEta[j][i] = VEta[i][j];
      counter++;
    }
  }

//   std::cout << "Counter = " << counter << std::endl;

//   for (int i=0;i<5;++i) {
//     printf("%19.15f %19.15f %19.15f %19.15f %19.15f\n",
//  	   VEta[i][0],VEta[i][1],VEta[i][2],VEta[i][3],VEta[i][4]);
//   }

  // Calculation of Jacobian 

  // some variables needed  for 
  // calculations of Jacobian matrix elements
  // Variables are set in cm !
  double cot_eta3 = 1.0/tan(rfit[3]);
  double sin_eta3 = sin(rfit[3]);
  double sin_eta4 = sin(rfit[4]);
  double cos_eta4 = cos(rfit[4]);
  double R = -100.0*sin_eta3/(0.299792458*bField*rfit[5]); 
  double rho = rfit[0];
  // Coordinates of the Circle central point (XC,YC)
  // XC = xRef + R*sin(phi)
  // YC = yRef - R*cos(phi)
  double XC = double(rho*cos(rfit[1]/rho) + R*sin_eta4);
  double YC = double(rho*sin(rfit[1]/rho) - R*cos_eta4);
  if (fitCode == 0) {
      XC = double(rfit[0] + R*sin_eta4);
      YC = double(rfit[1] - R*cos_eta4);
  }
  double RC2 = double(XC*XC+YC*YC);

  // X0 = Omega
  // Omega = -FCT*BField/(p*sin(Theta))
  dXdEta[0][0] = 0;
  dXdEta[0][1] = 0;
  dXdEta[0][2] = -cot_eta3/R;
  dXdEta[0][3] = 0;
  dXdEta[0][4] = 1./(R*rfit[5]);
    
  // X1 = TanLambda
  // TanLambda = cot(Theta)
  dXdEta[1][0] = 0;
  dXdEta[1][1] = 0;
  dXdEta[1][2] = -1/(sin_eta3*sin_eta3);
  dXdEta[1][3] = 0;
  dXdEta[1][4] = 0;
  
  // X2 = Phi0
  // Phi = atan2(-YC,-XC) - pi*q/2,
  // where q is the charge of particles
  double dXCdEta[5];
  double dYCdEta[5];
  
  dXCdEta[0] = -sin(rfit[1]/rho);
  dXCdEta[1] = double(0);
  dXCdEta[2] = R*cot_eta3*sin_eta4;
  dXCdEta[3] = R*cos_eta4;
  dXCdEta[4] = -R*sin_eta4/rfit[5];
  
  dYCdEta[0] = cos(rfit[1]/rho);
  dYCdEta[1] = double(0);
  dYCdEta[2] = -R*cot_eta3*cos_eta4;
  dYCdEta[3] = R*sin_eta4;
  dYCdEta[4] = R*cos_eta4/rfit[5];

  if (fitCode == 0) {
      dXCdEta[0] = double(1);
      dXCdEta[1] = double(0);
      dYCdEta[0] = double(0);
      dYCdEta[1] = double(1);      
  }

  for (int i=0;i<5;++i) 
    dXdEta[2][i] = (dYCdEta[i]*XC-dXCdEta[i]*YC)/RC2;
  
  
  // X3 = D0
  // D0 = R - XC/sin(Phi0) or D0 = R + YC/cos(Phi0)
  double dRdEta[5];
  dRdEta[0] = 0.0;
  dRdEta[1] = 0.0;
  dRdEta[2] = R*cot_eta3;
  dRdEta[3] = 0.0;
  dRdEta[4] = -R/rfit[5];
  
  // method using formulae XC = (R-d0)*sin(Phi0)
  //                       YC = (d0-R)*cos(Phi0) 


  float sinPhi0=sin(param[2]);
  float cosPhi0=cos(param[2]);
  float sin2Phi0=sinPhi0*sinPhi0;
  float cos2Phi0=cosPhi0*cosPhi0;
  
  if (fabs(cosPhi0)>fabs(sinPhi0)) {
    for (int i=0;i<5;++i) {
      dXdEta[3][i] = dRdEta[i] + (dYCdEta[i]*cosPhi0+YC*sinPhi0*dXdEta[2][i])/cos2Phi0;
    }
  }
  else {
    for (int i=0;i<5;++i) {
      dXdEta[3][i] = dRdEta[i] - (dXCdEta[i]*sinPhi0-XC*cosPhi0*dXdEta[2][i])/sin2Phi0;
    }    
  }
  

  // method using formulae : d0 = R-sqrt(XC**2+YC**2) , R > 0
  //                         d0 = R+sqrt(XC**2+YC**2) , R < 0
  // 
  //  double rXCRC = double(XC/RC);
  //  double rYCRC = double(YC/RC);
  //  double RC  = double(sqrt(RC2));

  //   for (int i=0;i<5;++i) {
  //     if (R > 0)
  //       dXdEta[3][i] =  dRdEta[i] - rXCRC*dXCdEta[i] - rYCRC*dYCdEta[i];
  //     else 
  //       dXdEta[3][i] =  dRdEta[i] + rXCRC*dXCdEta[i] + rYCRC*dYCdEta[i];
  //   }

  // X4 = Z0
  // Z0 = ZR + (Phi-PhiR)*sign(Omega)*TanLambda
  dXdEta[4][0] = -param[1]*R*dXdEta[2][0];
  dXdEta[4][1] = 1;
  dXdEta[4][2] = -dXdEta[1][2]*R*param[2] -param[1]*R*dXdEta[2][2] -param[1]*param[2]*dRdEta[2]
    +dXdEta[1][2]*R*rfit[4] +param[1]*dRdEta[2]*rfit[4];
  dXdEta[4][3] = -param[1]*R*dXdEta[2][3]+param[1]*R;
  dXdEta[4][4] = -param[1]*dRdEta[4]*param[2]-param[1]*R*dXdEta[2][4]+param[1]*dRdEta[4]*rfit[4];
  

  // Transformation of Cov matrix -->
  for (int i=0;i<5;++i) {
    for (int j=0;j<5;++j) {
      UX[i][j] = 0.0;
      for (int k=0;k<5;++k) {
	for (int l=0;l<5;++l)
	  UX[i][j] += VEta[k][l]*dXdEta[i][k]*dXdEta[j][l];
      }
    }
  }

//   for (int iC=0;iC<15;++iC) {
//     std::cout << eparam[iC] << " " ;
//   }
//   std::cout << std::endl;
  
//   for (int i=0;i<5;++i) {
//     printf("%19.15f %19.15f %19.15f %19.15f %19.15f\n",
//  	   UX[i][0],UX[i][1],UX[i][2],UX[i][3],UX[i][4]);
//   }

//   std::cout << std::endl;


  //.. (omega,tanl,phi,d0,z0) ->
  //.. (d0,phi0,omega,z0,tanl)

  // reordering errors according to the ILC Convention
  // T.Kraemer,  LC-DET-2006-004
  // Convertion from cm to mm           COV MATRIX
  //                                             ELEMENT          DIMENSION
//   eparam[0] = float(0.01*UX[0][0]);        // (Omega,Omega)           L-2
  
//   eparam[1] = float(0.1*UX[1][0]);         // (Omega,TanLambda)       L-1
//   eparam[2] = float(UX[1][1]);             // (TanLambda,TanLambda)   L0
    
//   eparam[3] = float(0.1*UX[2][0]);         // (Omega,Phi)             L-1
//   eparam[4] = float(UX[2][1]);             // (TanLambda,Phi)         L0
//   eparam[5] = float(UX[2][2]);             // (Phi,Phi)               L0
    
//   eparam[6] = float(UX[3][0]);             // (Omega,D0)              L0
//   eparam[7] = float(10*UX[3][1]);          // (TanLambda,D0)          L1
//   eparam[8] = float(10*UX[3][2]);          // (Phi,D0)                L1
//   eparam[9] = float(100*UX[3][3]);         // (D0,D0)                 L2
  
//   eparam[10] = float(UX[4][0]);            // (Omega,Z0)              L0
//   eparam[11] = float(10*UX[4][1]);         // (TanLambda,Z0)          L1
//   eparam[12] = float(10*UX[4][2]);         // (Phi,Z0)                L1
//   eparam[13] = float(100*UX[4][3]);        // (D0,Z0)                 L2
//   eparam[14] = float(100*UX[4][4]);        // (Z0,Z0)                 L2

  eparam[0] = float(100*UX[3][3]);         // (D0,D0)                 L2
  
  eparam[1] = float(10*UX[2][3]);          // (Phi,D0)                L1
  eparam[2] = float(UX[2][2]);             // (Phi,Phi)               L0
    
  eparam[3] = float(UX[0][3]);             // (Omega,D0)              L0
  eparam[4] = float(0.1*UX[0][2]);         // (Omega,Phi)             L-1
  eparam[5] = float(0.01*UX[0][0]);        // (Omega,Omega)           L-2
    
  eparam[6] = float(100*UX[4][3]);         // (Z0,D0)                 L2
  eparam[7] = float(10*UX[4][2]);          // (Z0,Phi)                L1
  eparam[8] = float(UX[4][0]);             // (Z0,Omega)              L0
  eparam[9] = float(100*UX[4][4]);         // (Z0,Z0)                 L2
  
  eparam[10] = float(10*UX[1][3]);         // (TanLambda,D0)          L1
  eparam[11] = float(UX[1][2]);            // (TanLambda,Phi)         L0
  eparam[12] = float(0.1*UX[1][0]);        // (TanLambda,Omega)       L-1
  eparam[13] = float(10*UX[1][4]);         // (TanLambda,Z0)          L1
  eparam[14] = float(UX[1][1]);            // (TanLambda,TanLambda)   L0

  //  std::cout << "sigma(D0)=" << sqrt(eparam[9]) << std::endl;
  //  std::cout << std::endl;

}

void MarlinTrackFit::setParametersForIPErrors(float * par) {

    for (int i=0;i<3;++i)
	_parIpReso[i] = par[i];
}
