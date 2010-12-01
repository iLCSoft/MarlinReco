#ifndef CreatePDFs_h
#define CreatePDFs_h 1

#include "Configure.h"

#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/ReconstructedParticle.h>

#include "PDF.h"

#include <vector>
#include <string>
#include <fstream>

using namespace lcio ;
using namespace marlin ;


/** === CreatePDFs Processor === <br>
 * written by Martin Ohlerich
 */


class PDF;

class CreatePDFs : public Processor {

 public:

  virtual Processor* newProcessor() { return new CreatePDFs ; }

  CreatePDFs() ;

  virtual void init() ;

  virtual void processRunHeader( LCRunHeader * run ) ;

  virtual void processEvent( LCEvent * evt ) ;

  virtual void check( LCEvent * evt ) ;

  virtual void end() ;


 protected:

  int _nRun ;
  int _nEvt ;
  int _nParticle ;

  std::string _recoCol ;     // reconstructed particle Collection
  std::string _mcCol ;       // MC particle Collection
  std::string _filename_c ;  // Name of file containing the pdfs (charged)
  std::string _filename_n ;  // Name of file containing the pdfs (neutral)
  std::string _relCol ;      // calohit to MC relation
  std::string _MCTracksRelCol; // Track to MC relation
  std::vector<int> _chStart, _nStart; // arrays of category event start
  std::vector<int> _pidCol, _NoOfPID;  // pid collection

  std::string curr_cat;
  std::vector<int> MyPidCol;
  int myindex, noClusterParticle;
  float _bField;

  PDF *pdf, *npdf;
  struct info_t{
    double px, py, pz;   // momentum
    double Eecal, Ehcal; // subdetector energies
    double ex;           // excentricity
    double dmean;        // mean distance of hits to helix
    double Edmean;        // energy weighted mean distance of hits to helix
    double Necal, Nhcal;    // numbers of hits in subdetectors
    double L1,L2,L3;        // deepest Layers of HCal with hits
    double EtoN_ecal, EtoN_hcal; // ratios
    double EecalToEtot; // ratio
    double EL1,EL2,EL3;        // first Layers of ECal with hits
  };
  info_t info;

  void init_info();


};

#endif


