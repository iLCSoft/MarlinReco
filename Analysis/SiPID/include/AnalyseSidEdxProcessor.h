/*
 * AnalyseSidEdxProcessor.h
 *
 *  Created on: Dec 15, 2016
 *      Author: strahinja
 */

#ifndef INCLUDE_ANALYSESIDEDXPROCESSOR_H_
#define INCLUDE_ANALYSESIDEDXPROCESSOR_H_


#include <string>

#include "marlin/Processor.h"

#include "lcio.h"
#include <TTree.h>
#include <TFile.h>

using namespace lcio ;
using namespace marlin ;


/**  AnalyseSidEdxProcessor for marlin.
 *
 *  Reads Si tracker dEdx data from the .slcio file
 *  and stores them in a root file for analysis.
 *
 * @author S. Lukic, Vinca, Belgrade
 * December 2016
 */

class AnalyseSidEdxProcessor : public Processor {

 public:
  virtual Processor*  newProcessor() { return new AnalyseSidEdxProcessor ; }

  AnalyseSidEdxProcessor() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   * Really?
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

  /*** Steerable parameters ***/
  std::string m_rootFileName;
  std::string m_trackColName;
  std::string m_linkColName;
  StringVec m_trkHitCollNames;

  TFile* rootfile;
  TTree* tree;

  /** ROOT output **/
  FloatVec pMC, thetaMC;
  FloatVec eTrack, dEdxTrack, dEdxError, eEvt;
  FloatVec d0, m; // Impact factor and mass of the associated particle
  FloatVec nTrkHits, nTrkRelatedParticles;
  FloatVec zTrackHit, xTrackHit, yTrackHit, eTrackHit, typeTrackHit;
  FloatVec zHit, xHit, yHit, eHit, typeHit;
  int nTracks;

  int lastRunHeaderProcessed;
} ;




#endif /* INCLUDE_ANALYSESIDEDXPROCESSOR_H_ */
