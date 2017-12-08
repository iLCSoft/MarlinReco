/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

/**  Digitizer for Simulated Hits in the Vertex Detector. Track segment within
 a layer is approximated with a line.
 It is divided by n subsegments, where the length of the segments  can be
 specified by the user via external Processor parameter.
 For each subsegment the charge  is simulated according to Landau distribution
 as expected for Silicon. The charge transfer from the middle point of track
 subsegment (referred  hereafter to as ionisation point) to the collection plane
 is performed taking into account Lorentz effect int the depleted zone. It is assumed that on the collection plane the electron cloud from each ionisation point is spread according to a gaussian distribution divided by the lateral distance. This Formular is taken from Sinevs Code, where it is not justified. The width of the distribution depends on the distance between ionisationpoint and collection plane. 
 The charge on each fired pixel is then calculated as a sum of contributions from n distributions. The VTX ladder is assumed to have matrix of rectangular pixels. The output of the processor is the collection of Reconstructed Tracker Hits.There are to reconstruction methods which can be chosen by the user. Lorentz effect in the depleted zone is corrected.
Input collections and prerequisites:
Processor requires collection of simulated vertex tracker hits. If such a collection with the user specified name does not exist processor takes no action. <h4>Output</h4> Processor produces an output collection of the Tracker Hits.  Collection has name "VTXTrackerHits".

 author: Stefan Uebelacker, RAL, 2008
 the code is based on:
 - VTXDigitizer written A. Raspereza, MPI (Munich)
 - Simulation of CCD detection process by M. Sinev (http://source.freehep.or g/jcvsweb/ilc/LCSIM/list/lcsim/src/org/lcsim/mc/CCDSim)
*/
 
#ifndef CCDDigitizer_h
#define CCDDigitizer_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/SimTrackerHit.h"
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/SimTrackerHitImpl.h"
#include <string>
#include <vector>
#include "MyG4UniversalFluctuationForSi.h"
#include "EVENT/LCIO.h"
#include <IMPL/LCCollectionVec.h>



#ifdef CCD_diagnostics
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IAxis.h>
#include <AIDA/ITree.h>

#include <marlin/Global.h>
#include <marlin/VerbosityLevels.h>
#include <marlin/AIDAProcessor.h>
#endif

using namespace lcio ;
using namespace marlin ;

struct IonisationPoint {
  double x;
  double y;
  double z;
  double eloss;
};


typedef std::vector<IonisationPoint> IonisationPointVec;
typedef std::vector<TrackerHitImpl*> TrackerHitImplVec;
typedef std::vector<SimTrackerHitImpl*> SimTrackerHitImplVec;



#define maxpixx 9// grid size in local(ladder) x axis
#define maxpixy 9// grid size in local(ladder) y axis
//grid size in which diffusion is computed, the charge, which diffuses outside this grid, is lost;value may be redefined again after sigmacoefficient or thickness of active layer is changed

//we approach the continual distribution of the diffusion with discrete points:
#define Numstepx 10// Number of points at which amplitude of diffusion is calculated within one pixel in x direction
#define Numstepy 10// Number of points at which amplitude of diffusion is calculated within one pixel in y direction
//according to first tests increasing the number of steps effects the performance of the processor only slightly and non- systematically

//if using a table, numhitstep must be adjusted when changing Numstep
//table
//#define numhitstepx 6//must be =(Numstepx/2)+1
//#define numhitstepy 6//must be =(Numstepy/2)+1
//#define numsigstep 20
//table//


class CCDDigitizer : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CCDDigitizer ; }
  
  
  CCDDigitizer() ;
  
  /**
   * Initialisation member function
   */
  virtual void init() ;
  
  /**
   * Processing of run header
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Processing of one event
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  /** Produces check plots
   */
  virtual void check( LCEvent * evt ) ; 
  
  /** Termination member function
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */
  std::string _colName{};
  std::string _outputCollectionName{};
  std::string _colVTXRelation{};

  /** Run number
   */
  int _nRun{};
  
  /** Event number
   */
  int _nEvt{};
  /** tangent of Lorentz angle
   */
 
  double _cutOnDeltaRays{};
  /** Diffusion coefficient in mm for nominla layer thickness
   */
  double _diffusionCoefficient{}; 



  int _numberOfLayers{};
  double _pixelSizeX{};
  double _pixelSizeY{};
  double _electronsPerKeV{};
  double _segmentDepth{};
  double _currentTotalCharge{};

  std::vector<int> _laddersInLayer{};
  std::vector<float> _layerRadius{};
  std::vector<float> _layerThickness{};
  std::vector<float> _layerHalfThickness{};
  std::vector<float> _layerLadderLength{};
  std::vector<float> _layerLadderHalfWidth{};
  std::vector<float> _layerPhiOffset{};
  std::vector<float> _layerActiveSiOffset{};
  std::vector<float> _layerHalfPhi{};
  std::vector<float> _layerLadderGap{};
  std::vector<float> _bkgdHitsInLayer{};
  std::vector<float> _layerLadderWidth{};

  int _currentLayer{};
  int _currentModule{};
  int _generateBackground{};
  double _currentParticleMomentum{};
  double _currentParticleEnergy{};
  double _currentParticleMass{};
  double _currentPhi{};

  double PI{},TWOPI{},PI2{};

  int _produceFullPattern{};
  int _numberOfSegments{};
  int _debug{};
  int _PoissonSmearing{};
  int _electronicEffects{};
  int _useMCPMomentum{};
 
  int _recmethod{};
  double _threshold{};
  double _saturation{};
  int _framesize{};
  int maxnionpoint{};
  
  double _currentLocalPosition[3]{};
  double _electronicNoise{};
  double _segmentLength{};

  IonisationPointVec _ionisationPoints{};
 

  MyG4UniversalFluctuationForSi * _fluctuate{};

   
  //Finds coordinates    
  void FindLocalPosition(SimTrackerHit * hit,
                         double * localPosition,
                         double * localDirection);

  void TransformToLab(double * xLoc, double * xLab);
  void ProduceIonisationPoints( SimTrackerHit * hit);
  void diffusion(double xdif,double ydif, double sigma);

  void ProduceHits(SimTrackerHitImplVec & simTrkVec);
  void TransformXYToCellID(double x, double y, 
                           int & ix, 
                           int & iy,double & xdif,double & ydif);  
  void TransformCellIDToXY(int ix, int iy,
                           double & x, double & y);
  void PoissonSmearer( SimTrackerHitImplVec & simTrkVec );
  void GainSmearer( SimTrackerHitImplVec & simTrkVec );
  void PrintInfo( SimTrackerHit * simTrkHit, TrackerHitImpl * recoHit);
  TrackerHitImpl * ReconstructTrackerHit(SimTrackerHitImplVec & simTrkVec );
  void TrackerHitToLab( TrackerHitImpl * recoHit );

  void generateBackground(LCCollectionVec * col);
  void settanlorentzangle(double B, double E, double mu, double T);
  void settanlorentzangleb(double B, double E, double mu, double T);

//  double _currentLocalPosition[3];
//   double _currentEntryPoint[3];
//   double _currentExitPoint[3];
 
  double _energyLoss{};
  std::vector <SimTrackerHitImplVec> _hitsInLayer{};

  double depdep{};
  double undep{};
  double epitaxdep{};

  double pxl[maxpixx][maxpixy]{};
  int midpixx{};
  int midpixy{};
  double stepx{};
  double stepy{};
  double xobsoffset{};
  double yobsoffset{};  
 
  double sigmacoefficient{};
  double sigmin{};
  double _difcoef{};
  double _efield{};
  double _biasvolt{};
  double _bfield{};
  double _T{};
  double _mu{};
  double TanLorentzAngle{};


#ifdef CCD_diagnostics
  double dirraw[3]{};
  double posraw[3]{};
  AIDA::IHistogram1D* histdist{};
  AIDA::IHistogram1D* histcluster{};
  AIDA::IHistogram2D* histclustxy{};
  AIDA::IHistogram1D* histcharge{};
  AIDA::IHistogram2D* histdistxy{};
  AIDA::IHistogram1D* histNionpoint{};
  AIDA::IHistogram1D* histzcoord{};
  AIDA::IHistogram1D* histenergy{};
  AIDA::IHistogram1D* histsignal{};
  AIDA::IHistogram1D* histsignalframe{};
  AIDA::IHistogram1D* histenergycentre{};
  int Nionpoint{};
#endif

 // table
 //  double sigstep;
//   double table [numsigstep][numhitstepx][numhitstepy][maxpixx][maxpixy];
//   void diffusiontable(double xdif,double ydif, double sigma);
//   void settable();
// table//

};



#endif

