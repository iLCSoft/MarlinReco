#ifndef AngularCorrection_dEdxProcessor_hh
#define AngularCorrection_dEdxProcessor_hh 1

#include "marlin/Processor.h"
#include "marlin/EventModifier.h"
#include "lcio.h"
#include <string>
#include <IMPL/LCEventImpl.h>
#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <IMPL/TrackImpl.h>

#include <marlin/Global.h>

using namespace lcio ;
using namespace marlin ;

/** AngularCorrection_dEdxProcessor <br>
 *  This processor calculates an extra correction to be applied to the computed dE/dx for every track.<br>
 * <h4>Input collections and prerequisites</h4>
 *  The processor requires a MarlinTrk Collection.<br>
 *  <h4>Output</h4>
 *  The calculated dE/dx is rewritten (ATTENTION) in the track collection.<br>
 *  @param _LDCTrackCollection - name of the input track collection <br>
 *  default: MarlinTrkTracks
 *  @param _par - parameter for angular correction<br>
 *  @author A. Irles, IFIC
 *  @version $Id$
 */

class AngularCorrection_dEdxProcessor : public Processor, public EventModifier  {

public:

  virtual Processor*  newProcessor() { return new AngularCorrection_dEdxProcessor ; }
  AngularCorrection_dEdxProcessor();
  virtual const std::string & name() const { return Processor::name() ; }

  virtual void modifyEvent( LCEvent * evt ) ;

  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run);
  //  virtual void processEvent( LCEvent * evt );
  virtual void check( LCEvent * evt );
  virtual void end();

private:
  AngularCorrection_dEdxProcessor(const AngularCorrection_dEdxProcessor&) = delete;
  AngularCorrection_dEdxProcessor& operator=(const AngularCorrection_dEdxProcessor&) = delete;

  //std::pair<float,float> getExtraCorrection(float dedx, float dedx_error, float trkcos);

  std::string _description = "";
  std::string _LDCTrackCollection = "";
  LCCollection* _LDCCol = NULL;

  std::vector<float> _par = {};

};


#endif
