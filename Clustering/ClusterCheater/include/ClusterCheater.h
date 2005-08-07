#ifndef FORMTRUECLUSTERSAR_H
#define FORMTRUECLUSTERSAR_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include "HelixClass.h"
#include <EVENT/CalorimeterHit.h>
#include <EVENT/MCParticle.h>

using namespace lcio ;
using namespace marlin ;


/** === Cluster Cheater === <br>
 *  This processor constructs true clusters, where a true <br>
 *  cluster is considered to comprise all hits attrubuted to <br>
 *  the same generator particle entering calorimeter. <br>
 *  Hits produced by back-scattered particles are assigned <br>
 *  to the primary which initiated this back-scattered particle. <br>
 *  Processor requires CalorimeterHit collections to be <br>
 *  specified with processor parameter CaloCollections. <br>
 *  Processor produces collection of clusters with the name <br>
 *  TrueClusters. Furthermore collection of MCParticle-to-Cluster <br>
 *  relations is stored in an LCIO event object. The name of <br>
 *  of the relation collection is TrueClusterToMCP. <br>
 *  In order to perform construction of true clusters <br>
 *  collection of relations between SimCalorimeterHits and CalorimeterHits <br>
 *  is required. The name of this collection is specified with <br>
 *  processor parameter RelCollection. <br> 
 *  To avoid collecting distant hits produced by neutrons or backscattered <br>
 *  particles, a proximity criteria is introduced in terms <br>
 *  of particle-to-hit distance cut. This is specified <br>
 *  via processor parameter ProximityCut. <br>
 *  Hit is assigned to the cluster produced by neutral particle if <br>
 *  the following condition is fulfilled <br>
 *  D{hit-to-Vertex}*Alpha < ProximityCut <br>.
 *  D{hit-to-IP} is the distance of the hit to the point <br>
 *  where neutral particle has been produced, and Alpha is the <br>
 *  angle between particle momentum vector and vector connecting <br>
 *  particle vertex and the hit. <br>
 *  For charged particles, the distance from the hit to <br>
 *  the helix associated with particle trajectory should 
 *  should not exceed ProximityCut. Clusters with the <br>
 *  number of hits less than predefined value (specified <br>
 *  with processor parameter MinimalHits) are not considered. <br>
 *  Magnetic field (in TESLA) is specified with processor <br>
 *  parameter MagneticField. <br>
 *    @author A. Raspereza (DESY)<br>
 *    @version $Id: ClusterCheater.h,v 1.4 2005-08-07 16:22:39 gaede Exp $<br>
 */
class ClusterCheater : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ClusterCheater ; }
  
  
  ClusterCheater() ;

  /** Initialization
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

  int _nRun ;
  int _nEvt ;
  int _ifBrahms;

  float _proximityCut;
  float _bField;
  int _minimal_hits;

  std::string _trueClustCollection;
  std::vector<std::string> _caloCollections;
  std::string _relCollection;

  float DistanceToChargeParticle(HelixClass * helix, CalorimeterHit * hit);
  float DistanceToNeutralParticle(MCParticle * par, CalorimeterHit * hit);
  HelixClass* AssignHelixToMCP( MCParticle * par);

} ;

#endif



