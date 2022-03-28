#ifndef EMShowerFinder_h
#define EMShowerFinder_h 1

#include <iostream>
#include <vector>
#include <map>

#include <marlin/Processor.h>
#include <lcio.h>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ClusterImpl.h>

#include <EVENT/CalorimeterHit.h>
#include "UTIL/CellIDDecoder.h"
#include <LCRTRelations.h>


#include "ClusterShapes.h"
#include "KITutil.h"
#include "Phys_Geom_Database.h" 


// MarlinCED is only used for debugging
#include <MarlinCED.h>


using namespace lcio ;
using namespace marlin ;

/**
\addtogroup Clustering Clustering
@{
*/

typedef struct {

  CalorimeterHit* ECALHit{};
  std::vector<PROTSEED2*> relatedCores{};
  std::vector<double> probabilitiesForThisECALHit{};
  std::vector<double> distancesToCoresForThisECALHit{};
  std::vector<double> estimatedEnergyPerCore{};

} ECALHitWithAttributes;



// integer runtime extension which flags CalorimeterHits in EMShowers
struct isPartOfEMShowerCandidate : LCIntExtension<isPartOfEMShowerCandidate> {};



/**
\addtogroup EMShowerFinder EMShowerFinder
@{
Initial version of a processor to find electro-magnetic showers.
 *    It is based on the KIT package and takes only ECAL hits into account.
 *    The output is a collection of clusters in which the electro-magnetic showers are stored. More docu will come soon.
 *    <br>
 *
 *    @author O. Wendt (DESY)
 *    @version
 *
 */
class EMShowerFinder : public Processor {

 public:

  virtual Processor* newProcessor() { return new EMShowerFinder ; }

  EMShowerFinder() ;

  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ;   
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

  
 private:

  
 protected:

  std::string _colNameECAL{};
  std::string _collectionNameOfEMShowerCandidates{};
  std::string _ToClean{};
  int _CleanCut{};
  int _N{};
  vector<float> _miipstep{};
  int _MinHit0{};
  int _MinHitSplit{};
  double _Rcut{};
  double _Distcut{}; 
  double _Coscut{};
  double _energyDeviationCut{};
  double _probabilityDensityCut{};

  int _debugLevel{};  
  int _drawOnCED{};

  int _nRun{};
  int _nEvt{};

} ;

/** @} @} */

#endif
