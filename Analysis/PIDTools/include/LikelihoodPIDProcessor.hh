#ifndef LikelihoodPIDProcessor_hh
#define LikelihoodPIDProcessor_hh 1

#include <string>

#include <marlin/Processor.h>

#include <EVENT/LCCollection.h>

using namespace lcio ;
using namespace marlin ;

class LikelihoodPID;

class LikelihoodPIDProcessor : public Processor{
public:
  virtual Processor*  newProcessor() { return new LikelihoodPIDProcessor ; }
  LikelihoodPIDProcessor();
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run);
  virtual void processEvent( LCEvent * evt );
  virtual void check( LCEvent * evt );
  virtual void end();
 
private:
  void createParticleIDClass(int parttype, ReconstructedParticle *part, PIDHandler &pidh, int algoID, float MVAoutput);
  
  LikelihoodPID *_myPID;
  std::string _description;
  std::string _inputPFOsCollection;
  std::string _PDFName;
  std::string _weightFileName;
  EVENT::FloatVec _energyBoundary;
  LCCollection* _pfoCol;
  std::vector<int> _pdgTable;
  std::vector<std::string> _particleNames;
  std::vector<std::string> _dEdxNames;

  bool _basicFlg;
  bool _dEdxFlg;
  bool _showerShapesFlg;
};

#endif 
