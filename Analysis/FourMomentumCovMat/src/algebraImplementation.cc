#include "TMath.h"
#include "streamlog/streamlog.h"
#include "marlin/Global.h"
#include <gear/BField.h>
#include <gearimpl/Vector3D.h>

#include "algebraImplementation.h"

using namespace lcio;
using namespace marlin;


// ====================================================================
int getCovMatrixMomenta(ReconstructedParticle const *par,
                        TMatrixD &covMatrixMomenta)
// Obtain covariance matrix on (px,py,pz,E) from the
// covariance matrix on helix parameters.
//
// I define the jacobian as the 3x4 matrix:
// (W = omega)
//
//
//       Dpx/DtanL     Dpy/DtanL      Dpz/DtanL     DE/DtanL
//
// J =   Dpx/DW        Dpy/DW         Dpz/DW        DE/DW
//
//       Dpx/DPhi      Dpy/DPhi       Dpz/DPhi      DE/DPhi
//
//
//
//         0           0         PT        Pz*Pz / E / tanL
//
// J =    -Px/W       -Py/W      -Pz/W     -P*P / E / W
//
//        -Py          Px         0        0
//
//
// Order in the covariance matrix on helix parameters:
//
//            tanL tanL         tanL omega        tanL phi
//
// Cov =      omega tanL        omega omega       omega phi
//
//            phi tanL          phi omega         phi phi
//
//
//
{
  TrackVec tracks = par->getTracks();
  if ( tracks.size() == 0) {
      streamlog_out(ERROR) << " Warning: this pfo has no associated tracks. \n";
      return -999;
  }

  Track* mytrack = tracks[0];
  double charge  = par->getCharge();
  if ( TMath::Abs(charge) < 0.5 ) {
    streamlog_out(ERROR) << " Warning: null charge particle!"
                         << " This case is not supported." << std::endl;
      return -999;
  }

  const int rows      = 3; // n rows jacobian
  const int columns   = 4; // n columns jacobian

  double bField       = Global::GEAR->getBField().at(gear::Vector3D(0,0,0)).z() ;
  const double CT     = 2.99792458E-4;

  double track_omega  = mytrack->getOmega();
  double track_tanL   = mytrack->getTanLambda();
  double track_phi    = mytrack->getPhi();
  double track_pt     = (CT*bField) / fabs(track_omega);
  double track_pz     = (track_pt*track_tanL);
  double track_p      = sqrt( (track_pt * track_pt)+(track_pz * track_pz) );
  double track_px     = track_pt * TMath::Cos(track_phi);
  double track_py     = track_pt * TMath::Sin(track_phi);
  double energy       = par->getEnergy();


  // Define array with jacobian matrix elements by columns
  Double_t jacobian_by_columns[rows*columns] = {
    0,
    -track_px/track_omega,-track_py,
    0,
    -track_py/track_omega,track_px,
    track_pt,-track_pz/track_omega,
    0,
    track_pz*track_pz/energy/track_tanL,
    -track_p*track_p/energy/track_omega,
    0
  };
  // construct the Jacobian using previous array ("F" if filling by columns,
  //"C", if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
  TMatrixD jacobian(rows,columns, jacobian_by_columns, "F");


  // Get the elements of the helix covariance matrix that we need.
  // i keep notation remembering that this matrix is a submatrix of a 5x5
  // (skip d0, z0 part)
  Double_t cov11 = (Double_t) mytrack->getCovMatrix()[2];
  Double_t cov22 = (Double_t) mytrack->getCovMatrix()[5];
  Double_t cov44 = (Double_t) mytrack->getCovMatrix()[14];

  Double_t cov41 = (Double_t) mytrack->getCovMatrix()[11];
  Double_t cov42 = (Double_t) mytrack->getCovMatrix()[12];
  Double_t cov21 = (Double_t) mytrack->getCovMatrix()[4];

  // hex covariance matrix by columns
  Double_t helix_cov_matrix_by_columns[rows*rows] = {cov44, cov42, cov41,
                                                     cov42, cov22, cov21,
                                                     cov41, cov21, cov11};

  TMatrixD covMatrix_helix(rows,rows, helix_cov_matrix_by_columns, "F");

  covMatrixMomenta.Mult( TMatrixD( jacobian,
                                   TMatrixD::kTransposeMult,
                                   covMatrix_helix) ,
                         jacobian);


  return 0;
}



// ====================================================================
int getCovMatrixMomenta(ReconstructedParticle const *par,  FloatVec &covP)
// Partial calculations make with doubles
{

  const int kspace_time_dim = 4;

  TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
  getCovMatrixMomenta(par, covMatrixMomenta);
  //  (px,py,pz,E)
  covP.at(0) = static_cast<float>(covMatrixMomenta(0,0) ); // x-x
  covP.at(1) = static_cast<float>(covMatrixMomenta(1,0) ); // y-x
  covP.at(2) = static_cast<float>(covMatrixMomenta(1,1) ); // y-y
  covP.at(3) = static_cast<float>(covMatrixMomenta(2,0) ); // z-x
  covP.at(4) = static_cast<float>(covMatrixMomenta(2,1) ); // z-y

  covP.at(5) = static_cast<float>(covMatrixMomenta(2,2) ); // z-z
  covP.at(6) = static_cast<float>(covMatrixMomenta(3,1) ); // e-x
  covP.at(7) = static_cast<float>(covMatrixMomenta(3,2) ); // e-y
  covP.at(8) = static_cast<float>(covMatrixMomenta(3,3) ); // e-z
  covP.at(9) = static_cast<float>(covMatrixMomenta(3,0) ); // e-e

  return 0;
}
