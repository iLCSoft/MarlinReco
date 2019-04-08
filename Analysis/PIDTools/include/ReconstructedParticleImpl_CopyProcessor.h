#ifndef ReconstructedParticleImpl_CopyProcessor_H
#define ReconstructedParticleImpl_CopyProcessor_H 1

#include "marlin/Processor.h"
#include "lcio.h"

using namespace lcio;
using namespace marlin;

/** ReconstructedParticleImpl_CopyProcessor <br>
 *
 * Copies an LCIOCollection of ReconstructedParticleImpl including copying each member attribute and each element of the member vectors.
 * In addition, an LCRelation is created, linking the original ReconstructedParticleImpl to their copies.
 * This indirect accessibility allows for principal compatibility with existing LCRelations, like recoMCTruthLink.
 * Necessary input are the names of the one input and two output collections.
 *
 * For every member attribute or member vector there is a boolean optional parameter with which the copying of the respective element can be switched off by setting the parameter to *false*.
 * This way, certain Marlin processors, that only modify the ReconstructedParticleImpl can be selectively run again instead of running the whole processor chain from the creation of the ReconstructedParticleImpl again.
 * It was initially written to selectively re-run the LikelihoodPIDProcessor (part of HighLevelReco of MarlinStdReco) on existing DST-files after certain changes in the processor parameters were made post production.
 * If the input collection does not exist in an event, or does not contain any elements, an empty instance of the output collection is added to the event.
 *
 * To use pre-existing LCRelations, use that relation to access the old ReconstructedParticleImpl, and then use the LCRelation created by this processor to access the new one (or the other way around).
 *
 * @author U. Einhaus, DESY
 * @version $1$
 */

class ReconstructedParticleImpl_CopyProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ReconstructedParticleImpl_CopyProcessor ; }
  
  ReconstructedParticleImpl_CopyProcessor();
  
  virtual ~ReconstructedParticleImpl_CopyProcessor() = default;
  
  ReconstructedParticleImpl_CopyProcessor(const ReconstructedParticleImpl_CopyProcessor&) = delete;
  
  ReconstructedParticleImpl_CopyProcessor& operator=(const ReconstructedParticleImpl_CopyProcessor&) = delete;
  
  
  virtual void init();
  
  //virtual void processRunHeader( LCRunHeader* run ) {};
  
  virtual void processEvent( LCEvent * evt );
  
  //virtual void check( LCEvent * evt ) {};
  
  virtual void end();
  

 protected:

  std::string _description = "";
  int _nEvt{};

  std::string _InputColName{};
  std::string _OutputColName{};
  std::string _RelationColName{};

  bool _copyType = true;
  bool _copyMomentum = true;
  bool _copyEnergy = true;
  bool _copyCovMatrix = true;
  bool _copyMass = true;
  bool _copyCharge = true;
  bool _copyReferencePoint = true;
  bool _copyParticleIDs = true;
  bool _copyParticleIDUsed = true;
  bool _copyGoodnessOfPID = true;
  bool _copyParticles = true;
  bool _copyClusters = true;
  bool _copyTracks = true;
  bool _copyStartVertex = true;


};

#endif
