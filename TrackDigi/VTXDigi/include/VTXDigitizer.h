/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef VTXDigitizer_h
#define VTXDigitizer_h 1

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


using namespace lcio ;
using namespace marlin ;

struct IonisationPoint {
  double x{};
  double y{};
  double z{};
  double eloss{};
};

struct SignalPoint {
  double x{};
  double y{};
  double sigmaX{};
  double sigmaY{};
  double charge{};

};

typedef std::vector<IonisationPoint> IonisationPointVec;
typedef std::vector<SignalPoint> SignalPointVec;
typedef std::vector<TrackerHitImpl*> TrackerHitImplVec;
typedef std::vector<SimTrackerHitImpl*> SimTrackerHitImplVec;

/**  Digitizer for Simulated Hits in the Vertex Detector. <br>
 * Digitization follows the procedure adopted in the CMS software package. <br>
 * For each Simulated Tracker Hit the intersection points with <br>
 * the inner and outer boundaries of sensitive Si layer are calculated. <br>
 * Track segment within a layer is approximated with the line. <br>
 * It is divided by n subsegments, where n can be specified by a user <br>
 * via external Processor parameter. For each subsegment the charge  <br>
 * is simulated according to Landau distribution as expected for Silicon. <br>
 * The charge transfer from the middle point of track subsegment (referred <br> 
 * hereafter to as ionisation point) <br>
 * to the outer collection plane is performed taking <br>
 * into account Lorentz effect in the magnetic field <br>
 * It is assumed that on the collection plane the electron cloud from <br>
 * each ionisation point is spread according to the Gaussian distribution <br>
 * whose width  is proportional to the square-root of the drift distance. <br>
 * The charge on each fired pixel is then calculated as a sum of contributions <br>
 * from n Gaussians. The VTX ladder is assumed to have matrix of rectangular pixels <br>
 * (In the future I plan to implement possibility of variying pixel size with <br>
 * the z coordinate). <br>
 * The output of the processor is the collection of Reconstructed Tracker Hits. <br>
 * Each reconstructed hit is assigned the position of the center-of-gravity of <br>
 * the cluster of fired pixels. Lorentz effect is corrected for. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collection of simulated vertex tracker hits. <br>
 * If such a collection with the user specified name does not exist <br>
 * processor takes no action.
 * <h4>Output</h4>
 * Processor produces an output collection of the Tracker Hits. 
 * Collection has name "VTXTrackerHits".
 * @param CollectionName name of input SimTrackerHit collection <br>
 * (default parameter value : "vxd01_VXD", taken from Mokka)
 * @param TanLorentz tangent of the Lorentz angle <br>
 * (default parameter value : 0.8) <br>
 * @param CutOnDeltaRays cut on the energy of delta-electrons (in MeV) <br>
 * (default parameter value : 0.03) <br>
 * @param Diffusion diffusion coefficient for the nominal active layer thickness (in mm) <br>
 * (default parameter value : 0.002) <br>
 * @param LayerThickness thickness of the active Silicon layer (in mm) <br>
 * (default parameter value : 0.03744) <br>
 * @param PixelSizeX pixel size along direction perpendicular to beam axis (in mm) <br>
 * (default value : 0.025) <br>
 * @param PixelSizeY pixel size along beam axis (in mm) <br>
 * (default value : 0.025) <br>
 * @param ElectronsPerMeV number of electrons produced per MeV of deposited energy <br>
 * (default parameter value : 270.3) <br>
 * @param Threshold threshold on charge deposited on one pixel (in electons) <br>
 * (default parameter value : 200.0) <br>
 * @param LaddersInLayer vector of integers, numbers of phi-ladders in each layer <br>
 * (default parameter values : 8, 8, 12, 16, 20; taken from Mokka database for VXD model vxd01) <br>
 * @param LadderRadius vector of doubles, radii of layers (in mm) <br>
 * (default parameter values : 15.301, 26.301, 38.301, 49.301, 60.301; taken from Mokka database 
 * for VXD model vxd01) <br>
 * @param ActiveLadderOffset vector of doubles, offset of active Si layer along phi-angle (in mm) for each layer <br>
 * (default parameter values : 1.455, 1.39866, 2.57163, 3.59295, 4.42245) <br>
 * @param LadderHalfWidth vector of doubles, half-width of the ladders in each layer (in mm)<br>
 * (default parameter values : 6.5, 11.0, 11.0, 11.0, 11.0; taken from Mokka database for VXD model vxd01) <br>
 * @param LadderGaps vector of doubles, gaps between two symmetric ladders (+/-z) in mm<br>
 * (default parameter values : 0.0, 0.04, 0.04, 0.04, 0.04; taken from Mokka database) <br>
 * @param PhiOffset vector of doubles, offset in phi angle for starting ladder in each layer <br>
 * (default parameter values : 0, 0, 0, 0, 0; taken from Mokka database for VXD model vxd01) <br>
 * @param SegmentLength segment length along track path which is used to subdivide track into segments (in mm).
 * The number of track subsegments is calculated as int(TrackLengthWithinActiveLayer/SegmentLength)+1 <br>
 * (default parameter value : 0.005) <br>
 * @param WidthOfCluster defines width in Gaussian sigmas to perform charge integration for 
 * a given pixel <br>
 * (default parameter value : 3.0) <br>
 * @param ElectronicEffects flag to switch on gaussian smearing of signal (electronic noise) <br>
 * (default parameter value : 1) <br>
 * @param ElectronicNoise electronic noise in electrons <br>
 * (default parameter value : 100) <br>
 * @param GenerateBackground flag to switch on pseudo-generation of pair background hits
 * Background hits are uniformly generated in cosQ and Phi <br>
 * (default parameter value : 0)
 * @param BackgroundHitsPerLayer expected mean value of background hits accumulated in each layer 
 * over readout time. This number is calculated as Number of background hits per bunch crossing times
 * number of bunch crossings over integration time. <br>
 * (default values : 34400 23900 9600 5500 3100 corresponding to 100 overlaid bunch crossings) <br>
 * @param Debug boolean variable, if set to 1, printout is activated <br>
 * (default parameter value : 0) <br>
 * <br>
 * @author A. Raspereza, MPI Munich
 */
class VTXDigitizer : public Processor {
  
 public:
  VTXDigitizer(const VTXDigitizer&) = delete;
  VTXDigitizer& operator=(const VTXDigitizer&) = delete;

  virtual Processor*  newProcessor() { return new VTXDigitizer ; }
  
  
  VTXDigitizer() ;
  
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
  double _tanLorentzAngle{};
  /** cut in MeV on delta electrons
   * used in simulation of charge for each ionisation point
   */
  double _cutOnDeltaRays{};
  /** Diffusion coefficient in mm for nominla layer thickness
   */
  double _diffusionCoefficient{}; 
  /** layer thickness
   */
 //  double _layerThickness;
//   /** layer half-thickness
//    */
//   double _layerHalfThickness;


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
  double _widthOfCluster{};

  double PI{},TWOPI{},PI2{};
  double SCALING{};

  int _produceFullPattern{};
  int _numberOfSegments{};
  int _debug{};
  int _PoissonSmearing{};
  int _electronicEffects{};
  int _useMCPMomentum{};
  int _removeDrays{};

  double _threshold{};
  double _currentLocalPosition[3]{};
  double _currentEntryPoint[3]{};
  double _currentExitPoint[3]{};
  double _electronicNoise{};
  double _segmentLength{};

  IonisationPointVec _ionisationPoints{};
  SignalPointVec _signalPoints{};

  MyG4UniversalFluctuationForSi * _fluctuate{};

  /** 
   * Finds coordinates 
   */
  void FindLocalPosition(SimTrackerHit * hit,
                         double * localPosition,
                         double * localDirection);

  void TransformToLab(double * xLoc, double * xLab);
  void ProduceIonisationPoints( SimTrackerHit * hit);
  void ProduceSignalPoints( );
  void ProduceHits(SimTrackerHitImplVec & simTrkVec);
  void TransformXYToCellID(double x, double y, 
                           int & ix, 
                           int & iy);  
  void TransformCellIDToXY(int ix, int iy,
                           double & x, double & y);
  void PoissonSmearer( SimTrackerHitImplVec & simTrkVec );
  void GainSmearer( SimTrackerHitImplVec & simTrkVec );
  void PrintInfo( SimTrackerHit * simTrkHit, TrackerHitImpl * recoHit);
  TrackerHitImpl * ReconstructTrackerHit(SimTrackerHitImplVec & simTrkVec );
  void TrackerHitToLab( TrackerHitImpl * recoHit );
  void PositionWithinCell(double x, double y,
                          int & ix, int & iy,
                          double & xCell, double & yCell);
  void generateBackground(LCCollectionVec * col);

  double _xLayerReco{},_yLayerReco{},_zLayerReco{};
  double _xLayerSim{},_yLayerSim{},_zLayerSim{};
  int _iLayer{};

  int _nCoveredX{},_nCoveredY{},_nCells{};

  int _totEntries{};
  double _totMomentum{};
  double _momX{},_momY{},_momZ{};
  double _eDep{};

  double _amplX[20]{};
  double _amplY[20]{};
  double _amplC[100]{};
  double _ampl{};
  double _amplMax{};
  double _eSum{};
  double _energyLoss{};
  double _clusterWidthX{},_clusterWidthY{};
  double _ampl33{},_ampl55{},_ampl77{};
  int _ncell33{},_ncell55{},_ncell77{};
  int _storeNtuple{};
  double _xLocalRecoCOG{},_yLocalRecoCOG{};
  double _xLocalRecoEdge{},_yLocalRecoEdge{};
  double _xLocalSim{},_yLocalSim{};

  std::vector <SimTrackerHitImplVec> _hitsInLayer{};

  


} ;

#endif
