#ifndef CLICPFOSELECTOR_H
#define CLICPFOSELECTOR 1

#include "PfoUtilities.h"
#include "marlin/Processor.h"
#include <EVENT/ReconstructedParticle.h>
#include "lcio.h"
#include "TrackHitPair.h"
#include "HelixClass.h"
#include <string>
#include <map>
#include <set>
#include <algorithm>

#define FORMATTED_OUTPUT_TRACK_CLUSTER(out, N1, E1,E2,E3,N2,E4,N3,E5,E6,E7) \
    out <<                                                                                      \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthInt)      <<    N2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E4        <<                                   \
    std::right << std::setw(widthInt  )    <<    N3        <<                                   \
    std::right << std::setw(widthFloat)    <<    E5        <<                                   \
    std::right << std::setw(widthFloat)    <<    E6        <<                                   \
    std::right << std::setw(widthFloat)    <<    E7  << std::endl

using namespace lcio ;
using namespace marlin ;

/** === CLICPfoSelector Processor === <br>
 * Processor to select good pfos based on timing
 */

class CLICPfoSelector : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CLICPfoSelector ; }  
  CLICPfoSelector() ;  
  float  TimeAtEcal(const Track* pTrack, float &tof);
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

 protected:

  void CleanUp(); 

  int _nRun=-1;
  int _nEvt=-1;

  float _bField=0.0;

  std::string     m_inputPfoCollection{};                 ///< Input PFO collection name
  std::string     m_selectedPfoCollection{};              ///< Output PFO collection name
  int             m_monitoring=0;                         ///< Whether to display monitoring information
  int             m_displaySelectedPfos=0;                ///< Whether to display monitoring information concerning selected pfos
  int             m_displayRejectedPfos=0;                ///< Whether to display monitoring information concerning rejected pfos
  float           m_monitoringPfoEnergyToDisplay=1.0;     ///< Minimum pfo energy in order to display monitoring information
  int             m_correctHitTimesForTimeOfFlight=0;     ///< Correct hit times for straight line time of flight
  int             m_checkProtonCorrection=0;              ///< Check proton hypothesis
  int             m_checkKaonCorrection=0;                ///< Check charged hypothesis
  int             m_keepKShorts=1;                        ///< Keep kshorts
  int             m_useNeutronTiming=0;                   ///< Attempt to make a (dubious) neutron timing correction
  float           m_minimumEnergyForNeutronTiming=1.0;    ///< Minimum energy for attempted neutron timing correction

  float           m_farForwardCosTheta=0.975;             ///< Value of cos theta identifying the detector forward region
  float           m_ptCutForTightTiming=0.75;             ///< The pt value below which tight timing cuts are used
  
  float           m_photonPtCut=0.0;                      ///< The basic pt cut for a photon pfo
  float           m_photonPtCutForLooseTiming=4.0;        ///< The photon pt value below which tight timing cuts are used
  float           m_photonLooseTimingCut=2.0;             ///< The photon loose high timing cut
  float           m_photonTightTimingCut=1.0;             ///< The photon tight high timing cut
  
  float           m_chargedPfoPtCut=0.0;                  ///< The basic pt cut for a charged hadron pfo
  float           m_chargedPfoPtCutForLooseTiming=4.0;    ///< The charged hadron pt value below which tight timing cuts are used
  float           m_chargedPfoLooseTimingCut=3.0;         ///< The charged hadron loose high timing cut
  float           m_chargedPfoTightTimingCut=1.5;         ///< The charged hadron tight high timing cut
  float           m_chargedPfoNegativeLooseTimingCut=-1.0;///< The charged hadron loose low timing cut
  float           m_chargedPfoNegativeTightTimingCut=-0.5;///< The charged hadron tight low timing cut
  
  float           m_neutralHadronPtCut=0.0;               ///< The basic pt cut for a neutral hadron pfo
  float           m_neutralHadronPtCutForLooseTiming=8.0; ///< The neutral hadron pt value below which tight timing cuts are used
  float           m_neutralHadronLooseTimingCut=2.5;      ///< The neutral hadron loose high timing cut
  float           m_neutralHadronTightTimingCut=1.5;      ///< The neutral hadron tight high timing cut
  float           m_neutralFarForwardLooseTimingCut=2.0;  ///< The neutral hadron loose high timing cut for the forward region
  float           m_neutralFarForwardTightTimingCut=1.0;  ///< The neutral hadron tight high timing cut for the forward region
  float           m_photonFarForwardLooseTimingCut=2.0;   ///< The photon loose high timing cut for the forward region
  float           m_photonFarForwardTightTimingCut=1.0;   ///< The photon tight high timing cut for the forward region

  float           m_hCalBarrelLooseTimingCut=20.0;        ///< The loose timing cut for hits predominantly in hcal barrel
  float           m_hCalBarrelTightTimingCut=10.0;        ///< The tight timing cut for hits predominantly in hcal barrel
  float           m_hCalEndCapTimingFactor=1.0;           ///< Factor by which high timing cut is multiplied for hcal barrel hits
  float           m_neutralHadronBarrelPtCutForLooseTiming=3.5;///< pt above which loose timing cuts are applied to neutral hadrons in barrel
  
  int             m_minECalHitsForTiming=5;               ///< Minimum ecal hits in order to use ecal timing info
  int             m_minHCalEndCapHitsForTiming=5;         ///< Minimum hcal endcap hits in order to use hcal endcap timing info
  
  int             m_useClusterLessPfos=1;                 ///< Whether to accept any cluster-less pfos
  float           m_minMomentumForClusterLessPfos=0.5;    ///< Minimum momentum for a cluster-less pfo
  float           m_maxMomentumForClusterLessPfos=2.0;    ///< Minimum momentum for a cluster-less pfo
  float           m_minPtForClusterLessPfos=0.5;          ///< Minimum pT for a cluster-less pfo
  float           m_clusterLessPfoTrackTimeCut=10.0;      ///< Maximum arrival time at Ecal for cluster-less pfo
  float           m_forwardCosThetaForHighEnergyNeutralHadrons=0.95;   ///< Forward region of HCAL where timing cuts will be applied to all neutral hadrons
  float           m_forwardHighEnergyNeutralHadronsEnergy=10.0;        ///< Energy cut for specific HCAL timing requirements cuts for neutral hadrons

} ;

#endif



