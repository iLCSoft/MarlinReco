#ifndef Compute_dEdxProcessor_hh
#define Compute_dEdxProcessor_hh 1

#include <EVENT/LCCollection.h>
#include <TH2.h>
#include <marlin/Processor.h>
#include <random>
#include <string>
#include <vector>

using namespace lcio;
using namespace marlin;

/** Compute dE/dx Processor <br>
 *  This processor calculates the dE/dx for every track.<br>
 * <h4>Input collections and prerequisites</h4>
 *  The processor requires a collection of tracks that contains the corresponding track hits.<br>
 *  Every track hit that lies within the boundaries of the TPC is used.<br>
 *  dE is the deposited energy of the hit.<br>
 *  dx is the distance between the hits and can be calculated in 3 different ways.<br>
 *  A truncation of the hits with the lowest and highest dE7Dx values is performed.<br>
 *  Then the mean is calculated (truncation-mean method).<br>
 *  <h4>Output</h4>
 *  The calculated dE/dx is attached to the track.<br>
 *  This is only possible if the the track collection allows write access.<br>
 *  Bethe-Bloch histograms (root TH2D) can be generated for every dx strategy.<br>
 *  Both outputs are optional.<br>
 *  @param _LDCTrackCollection - name of the input track collection <br>
 *  default: MarlinTrkTracks
 *  @param _writedEdx - flag indicating if calculated dE/dx should be attached to track<br>
 *  If fully reconstructed tracks are used as input collection this can be switched off to only generate histograms.<br>
 *  default: true<br>
 *  @param _energyLossErrorTPC - the dE/dx resolution<br>
 *  default: 0.054  (5.4%)<br>
 *  @param _lowerTrunFrac - lower truncation fraction for truncated-mean method<br>
 *  The hits with the lowest <_lowerTrunFrac> dE/dx values are rejected.<br>
 *  default: 0.08  (8%; ALEPH: 8%)<br>
 *  @param _upperTrunFrac - upper truncation fraction for truncated-mean method<br>
 *  The hits with the highest <_upperTrunFrac> dE/dx values are rejected.<br>
 *  default: 0.3  (30%; ALEPH: 40%)<br>
 *  @param _isSmearing - flag indicating if additional smearing should be applied<br>
 *  This compensates for a 'too good' processor outcome, compared to test beam results.<br>
 *  default: false<br>
 *  @param _smearingFactor - width of the Gaussian function used for the additional smearing<br>
 *  default: 0.035  (3.5%)<br>
 *  @param _ncorrpar - parameter for number-of-hits correction<br>
 *  default: 1.468<br>
 *  @param _acorrpar - parameters for angular correction<br>
 *  default: {0.635762, -0.0573237}<br>
 *  @param _errexp - scaling exponents of the dE/dx error for path length and number of hits, respectively<br>
 *  default: {-0.34, -0.45}<br>
 *  @param _dxStrategy - ID specifying which strategy for calculating dx should be used<br>
 *  Strategy 1: hit-to-hit distance<br>
 *  Strategy 2: hit-to-hit path length of projected hits  (do not use at the moment)<br>
 *  Strategy 3: path over hit row<br>
 *  default: 1<br>
 *  If none of the above is chosen, the processor defaults to 1.<br>
 *  @param _StratCompHist - flag indicating if Bethe-Bloch histograms for each dx strategy should created.<br>
 *  default: false<br>
 *  @param _StratCompHistWeight - flag indicating if Bethe-Bloch histograms (if chosen) should be filled with a
 * sqrt(number-of-track-hits) weighting.<br> default: false (-> weight for each track = 1)<br>
 *  @param _StratCompHistFiles - file names of the generated dx strategy comparison histograms (if chosen).<br>
 *  The respective strategy number and '.png' is added.<br>
 *  default: dEdx_Histo_Strategy  (-> "dEdx_Histo_Strategy1.png", etc.)<br>
 *  @author M. Kurata, KEK
 *  adapted by U. Einhaus, DESY
 *  @version $Id$
 */

class Compute_dEdxProcessor : public Processor {
public:
  virtual Processor* newProcessor() { return new Compute_dEdxProcessor; }
  Compute_dEdxProcessor();
  virtual void init();
  virtual void processRunHeader(LCRunHeader* run);
  virtual void processEvent(LCEvent* evt);
  virtual void check(LCEvent* evt);
  virtual void end();

private:
  Compute_dEdxProcessor(const Compute_dEdxProcessor&) = delete;
  Compute_dEdxProcessor& operator=(const Compute_dEdxProcessor&) = delete;

  std::pair<double, double> CalculateEnergyLoss(TrackerHitVec& hitVec, Track* trk);
  double getNormalization(double dedx, float hit, double trkcos);
  double getSmearing(double dEdx);

  std::string _description = "";
  std::string _LDCTrackCollection = "";
  LCCollection* _LDCCol = NULL;
  bool _writedEdx = true;

  float _energyLossErrorTPC = 0;
  float _lowerTrunFrac = 0;
  float _upperTrunFrac = 0;
  float _ncorrpar = 0;
  std::vector<float> _acorrpar = {};
  std::vector<float> _errexp = {};
  int _dxStrategy = 0;
  bool _StratCompHist = false;
  bool _StratCompHistWeight = false;
  std::string _StratCompHistFiles = "";

  // smearing
  std::random_device seed_gen{};
  std::default_random_engine* engine = NULL;
  std::uniform_real_distribution<> dist{};
  bool _isSmearing = 0;
  float _smearingFactor = 0;

  // geometry
  float _TPC_inner = 0;
  float _TPC_outer = 0;
  float _TPC_padHeight = 0;
  float _bField = 0;

  // root histograms for dx strategy comparison
  TH2* _BBHist_Strategy1{};
  TH2* _BBHist_Strategy2{};
  TH2* _BBHist_Strategy3{};
};

#endif
