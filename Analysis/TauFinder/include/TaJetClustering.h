#ifndef TaJetClustering_h
#define TaJetClustering_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

using namespace lcio ;
using namespace marlin ;

class TaJetClustering : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TaJetClustering ; }
  
  
  TaJetClustering() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */
	std::string _pfoCollectionName{};
	std::string _tauCollectionName{};
	std::string _remainCollectionName{};

	// parameters for tau clustering
	double _tauMass{};
	double _tauCosAngle{};

	// selection
	int _noSelection{};

	// primary cuts
	double _minimumTrackEnergy{};
	double _minimumJetEnergy{};
	double _minimumTrackEnergyAssoc{};
	int _acceptFlexibleLowEnergyTrack{};
	//int _maxTracks; // now fixed: allowed only 1 or 3

	// three prong cut
	double _maximumNeutralEnergyInThreeProng{};

	// cone cuts
	double _coneMinCosAngle{};
	double _coneMaxCosAngle{};
	double _coneMaxEnergyFrac{};

	// impact parameter cuts
	double _ipCutsMinimumTrackEnergy{};
	double _ipCutsSigmaOneProng{};			// and cut, d0 or z0
	double _ipCutsSigmaThreeProngNoNeutrals{};	// and cut, one of three, d0 or z0
	double _ipCutsSigmaThreeProngWithNeutrals{};	// or cut with cone, one of three, d0 or z0
	double _ipMaxAbsolute{}; // d0 or z0, max allowed absolute value in mm

	// lepton ID cuts: only applicable when ntr==1 && eneutral<1 GeV
	double _leptonCutsMinimumTrackEnergy{};
	double _muMaxFracEcal{};
	double _muMaxCalByTrack{};
	double _eMinFracEcal{};
	double _eMinCalByTrack{};
	double _eMaxCalByTrack{};

	//int _maxTaus; // throw away too many taus: energy order after all cuts
	//int _returnCutTaus; // put cut taus to remain collection or not

	int _nRun{};
	int _nEvt{};
} ;

#endif



