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
  LikelihoodPID *_myPID;
  std::string _description;
  std::string _inputPFOsCollection;
  std::string _PDFName;
  EVENT::FloatVec _energyBoundary;
  LCCollection* _pfoCol;
  std::vector<int> _pdgTable;
  std::vector<std::string> _particleNames;
};

#endif 
