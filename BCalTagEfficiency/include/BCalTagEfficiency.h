#ifndef BCalTagEfficiency_h
#define BCalTagEfficiency_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


// includes for the BcEnergyDensity

#include "BcEnergyDensity.h"

// root includes 

#include <TTree.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TClonesArray.h>


using namespace lcio ;
using namespace marlin ;


/**  Marlin processor to calculate BCAL tagging efficiency
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: BCalTagEfficiency.h,v 1.2 2009-01-16 11:09:12 engels Exp $ 
 */

class BCalTagEfficiency : public Processor {
 
 public:
  
  virtual Processor*  newProcessor() { return new BCalTagEfficiency ; }
  
  
  BCalTagEfficiency() ;
  
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
  std::string _colName ;
  std::string _BCALcolName ;
  std::string _particle;

  int _nRun ;
  int _nEvt ;


  // root variables to create tree 

  TFile* rootfile;
  TTree* tree;
  std::string rootFileName;

  // per particle
  int mcp;
  enum {MCP = 100};
  int pdg[MCP];
  float energy[MCP], ePrime[MCP], pxPrime[MCP], pxIP[MCP], pyIP[MCP], pzIP[MCP];
  float phiIP[MCP], theIP[MCP];
  float posx[MCP], posy[MCP], posz[MCP]; 
  float radius[MCP], phi[MCP], ebkg[MCP]; 
  float efficiency[MCP], rand[MCP];
  int tag[MCP];
  float scaleP[MCP];

  // parameters of detector
  std::string backgroundfilename;
  float densityScaling;
  float bField;
  float eBeam;
  float zbcal;
  float thresholdMin, thresholdMax;

  //
  bool detectAll;
  bool smearEnergy;
  
  // background handler
  BcEnergyDensity *bc_en;

} ;

#endif



