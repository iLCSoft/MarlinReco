#ifndef CLICPFOSELECTOR_H
#define CLICPFOSELECTOR 1

#include "marlin/Processor.h"
#include <EVENT/ReconstructedParticle.h>
#include "lcio.h"
#include "TrackHitPair.h"
#include "HelixClass.h"
#include <string>
#include <map>
#include <set>
#include <algorithm>

#define FORMATTED_OUTPUT_TRACK_CLUSTER(N1, E1,E2,E3,N2,E4,N3,E5,E6,E7)	                        \
    std::cout  <<                                                                               \
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

#define FORMATTED_OUTPUT_CLUSTER(N1, E1,E2,E3,N3,E5,E6,E7)	                                \
    std::cout  <<                                                                               \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthInt)      <<    " - "     <<                                   \
    std::right << std::setw(widthFloat)    <<    " - "     <<                                   \
    std::right << std::setw(widthInt  )    <<    N3        <<                                   \
    std::right << std::setw(widthFloat)    <<    E5        <<                                   \
    std::right << std::setw(widthFloat)    <<    E6        <<                                   \
    std::right << std::setw(widthFloat)    <<    E7        << std::endl

#define FORMATTED_OUTPUT_TRACK(N1, E1,E2,E3,N2,E4)                             			\
    std::cout  <<                                                                               \
    std::right << std::setw(widthInt)      <<    N1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E1        <<                                   \
    std::right << std::setw(widthFloat)    <<    E2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E3        <<                                   \
    std::right << std::setw(widthInt)      <<    N2        <<                                   \
    std::right << std::setw(widthFloat)    <<    E4        << std::endl

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
  void   GetClusterTimes(const Cluster* cluster, float &meanTime, int &nCaloHitsUsed, float &meanTimeEcal, int &nEcal, float &meanTimeHcalEndcap, int &nHcalEnd);
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

 protected:

  void CleanUp(); 

  int _nRun ;
  int _nEvt ;

  float PI, PIOVER2, TWOPI;
  float _bField;
  int   m_debug;
  LCEvent * _evt;

  std::string     m_inputPfoCollection;                           ///< Input PFO collection name
  std::string     m_selectedPfoCollection;                        ///< Output PFO collection name
  int             m_monitoring;                                   ///< Whether to display monitoring information
  int             m_displaySelectedPfos;                          ///< Whether to display monitoring information concerning selected pfos
  int             m_displayRejectedPfos;                          ///< Whether to display monitoring information concerning rejected pfos
  float           m_monitoringPfoEnergyToDisplay;                 ///< Minimum pfo energy in order to display monitoring information 
  int             m_correctHitTimesForTimeOfFlight;               ///< Correct hit times for straight line time of flight  
  int             m_checkProtonCorrection;                        ///< Check proton hypothesis
  int             m_checkKaonCorrection;                          ///< Check charged hypothesis
  int             m_keepKShorts;                                  ///< Keep kshorts
  int             m_useNeutronTiming;                             ///< Attempt to make a (dubious) neutron timing correction
  float           m_minimumEnergyForNeutronTiming;                ///< Minimum energy for attempted neutron timing correction

  float           m_farForwardCosTheta;                           ///< Value of cos theta identifying the detector forward region
  float           m_ptCutForTightTiming;                          ///< The pt value below which tight timing cuts are used
  
  float           m_photonPtCut;                                  ///< The basic pt cut for a photon pfo
  float           m_photonPtCutForLooseTiming;                    ///< The photon pt value below which tight timing cuts are used
  float           m_photonLooseTimingCut;                         ///< The photon loose high timing cut
  float           m_photonTightTimingCut;                         ///< The photon tight high timing cut
  
  float           m_chargedPfoPtCut;                              ///< The basic pt cut for a charged hadron pfo
  float           m_chargedPfoPtCutForLooseTiming;                ///< The charged hadron pt value below which tight timing cuts are used
  float           m_chargedPfoLooseTimingCut;                     ///< The charged hadron loose high timing cut
  float           m_chargedPfoTightTimingCut;                     ///< The charged hadron tight high timing cut
  float           m_chargedPfoNegativeLooseTimingCut;             ///< The charged hadron loose low timing cut
  float           m_chargedPfoNegativeTightTimingCut;             ///< The charged hadron tight low timing cut
  
  float           m_neutralHadronPtCut;                           ///< The basic pt cut for a neutral hadron pfo
  float           m_neutralHadronPtCutForLooseTiming;             ///< The neutral hadron pt value below which tight timing cuts are used
  float           m_neutralHadronLooseTimingCut;                  ///< The neutral hadron loose high timing cut
  float           m_neutralHadronTightTimingCut;                  ///< The neutral hadron tight high timing cut
  float           m_neutralFarForwardLooseTimingCut;              ///< The neutral hadron loose high timing cut for the forward region
  float           m_neutralFarForwardTightTimingCut;              ///< The neutral hadron tight high timing cut for the forward region
  float           m_photonFarForwardLooseTimingCut;               ///< The photon loose high timing cut for the forward region
  float           m_photonFarForwardTightTimingCut;               ///< The photon tight high timing cut for the forward region

  float           m_hCalBarrelLooseTimingCut;                     ///< The loose timing cut for hits predominantly in hcal barrel
  float           m_hCalBarrelTightTimingCut;                     ///< The tight timing cut for hits predominantly in hcal barrel
  float           m_hCalEndCapTimingFactor;                       ///< Factor by which high timing cut is multiplied for hcal barrel hits
  float           m_neutralHadronBarrelPtCutForLooseTiming;       ///< pt above which loose timing cuts are applied to neutral hadrons in barrel
  
  int             m_minECalHitsForTiming;                         ///< Minimum ecal hits in order to use ecal timing info
  int             m_minHCalEndCapHitsForTiming;                   ///< Minimum hcal endcap hits in order to use hcal endcap timing info
  
  int             m_useClusterLessPfos;                           ///< Whether to accept any cluster-less pfos
  float           m_minMomentumForClusterLessPfos;                ///< Minimum momentum for a cluster-less pfo
  float           m_maxMomentumForClusterLessPfos;                ///< Minimum momentum for a cluster-less pfo
  float           m_minPtForClusterLessPfos;                      ///< Minimum pT for a cluster-less pfo  
  float           m_clusterLessPfoTrackTimeCut;                   ///< Maximum arrival time at Ecal for cluster-less pfo
  float           m_forwardCosThetaForHighEnergyNeutralHadrons;   ///< Forward region of HCAL where timing cuts will be applied to all neutral hadrons
  float           m_forwardHighEnergyNeutralHadronsEnergy;        ///< Energy cut for specific HCAL timing requirements cuts for neutral hadrons

 private:
  typedef std::vector<ReconstructedParticle*> PfoList;
  static bool PfoSortFunction(ReconstructedParticle* lhs,ReconstructedParticle* rhs); 

} ;

#endif



