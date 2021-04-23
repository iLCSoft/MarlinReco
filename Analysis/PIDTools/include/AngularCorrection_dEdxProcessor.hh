#ifndef AngularCorrection_dEdxProcessor_hh
#define AngularCorrection_dEdxProcessor_hh 1

#include <string>
#include <vector>
#include <random>
#include <marlin/Processor.h>
#include <EVENT/LCCollection.h>
#include <TH2.h>

using namespace lcio ;
using namespace marlin ;

/** AngularCorrection_dEdxProcessor <br>
 *  This processor calculates an extra correction to be applied to the computed dE/dx for every track.<br>
 * <h4>Input collections and prerequisites</h4>
 *  The processor requires a MarlinTrk Collection.<br>
 *  <h4>Output</h4>
 *  The calculated dE/dx is rewritten (brute force) in the track collection.<br>
 *  This is only possible if Marlin setting AllowToModifyEvent is set to true.<br>
 *  @param _LDCTrackCollection - name of the input track collection <br>
 *  default: MarlinTrkTracks
 *  @param _writedEdx - flag indicating if calculated dE/dx should be attached to track<br>
 *  If fully reconstructed tracks are used as input collection this can be switched off to only generate histograms.<br>
 *  default: false<br>
 *  @param _par - parameter for angular correction<br>
 *  @author A. Irles, IFIC
 *  @version $Id$
 */

class AngularCorrection_dEdxProcessor : public Processor{
public:
  virtual Processor*  newProcessor() { return new AngularCorrection_dEdxProcessor ; }
  AngularCorrection_dEdxProcessor();
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run);
  virtual void processEvent( LCEvent * evt );
  virtual void check( LCEvent * evt );
  virtual void end();
 
private:
  AngularCorrection_dEdxProcessor(const AngularCorrection_dEdxProcessor&) = delete;
  AngularCorrection_dEdxProcessor& operator=(const AngularCorrection_dEdxProcessor&) = delete;

  //std::pair<float,float> getExtraCorrection(float dedx, float dedx_error, float trkcos);

  std::string _description = "";
  std::string _LDCTrackCollection = "";
  LCCollection* _LDCCol = NULL;
  bool _writedEdx = true;

  std::vector<float> _par = {};

};

#endif
