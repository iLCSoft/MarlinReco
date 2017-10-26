#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cmath>
#include <math.h>

#include <marlin/Global.h>
#include <marlin/Processor.h>

#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/TrackImpl.h>

#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>

#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"

#include "Compute_dEdxProcessor.hh"

Compute_dEdxProcessor aCompute_dEdxProcessor ;

Compute_dEdxProcessor::Compute_dEdxProcessor()
  : Processor("Compute_dEdxProcessor") {
  
  // Processor description
  _description = "dE/dx calculation using truncation method" ;
  
  registerProcessorParameter( "EnergyLossErrorTPC",
			      "Fractional error of dEdx in the TPC",
			      _energyLossErrorTPC,
			      float(0.054));
  
  registerProcessorParameter( "isSmearing",
			      "Flag for dEdx Smearing",
			      _isSmearing,
			      bool(0));

  registerProcessorParameter( "smearingFactor",
			      "Smearing factor for dEdx smearing",
			      _smearingFactor,
			      float(0.035));

  std::vector<float> apar;
  apar.push_back(0.635762);
  apar.push_back(-0.0573237);
  registerProcessorParameter( "AngularCorrectionParameters",
			      "parameters for angular correction",
			      _acorrpar,
			      apar);

  registerProcessorParameter( "NumberofHitsCorrectionParameters",
			      "parameters for number of hits correction",
			      _ncorrpar,
			      float(1.468));
  
  registerInputCollection(LCIO::TRACK,
			  "LDCTrackCollection",
			  "LDC track collection name",
			  _LDCTrackCollection,
			  std::string("MarlinTrkTracks"));
} 

void Compute_dEdxProcessor::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
  //_mydEdx = new ComputeddEdx(); 
  
  //get TPC inner radius
  dd4hep::Detector& lcdd = dd4hep::Detector::getInstance();
  dd4hep::DetElement tpcDE = lcdd.detector("TPC") ;
  const dd4hep::rec::FixedPadSizeTPCData* tpc = tpcDE.extension<dd4hep::rec::FixedPadSizeTPCData>() ;
  _TPC_inner = tpc->rMinReadout ;

  
  printParameters();
  
}

void Compute_dEdxProcessor::processRunHeader( LCRunHeader* ) { 
} 

void Compute_dEdxProcessor::processEvent( LCEvent * evt ) { 
  _LDCCol = evt->getCollection( _LDCTrackCollection ) ;
  int nTrkCand = _LDCCol->getNumberOfElements();
  
  for (int iTRK=0;iTRK<nTrkCand;++iTRK) {
    
    TrackImpl * trkCand = (TrackImpl*) _LDCCol->getElementAt( iTRK );
    TrackerHitVec trkHits = trkCand->getTrackerHits();
    
    //calculate dEdx here when the subtrack is the TPC track!
    //silicon dE/dx can also be calculated. is it necessary??
    float dedx=0.0;
    float unum=0.0;
    
    //get tpc hits only
    TrackerHitVec tpcHitVec;
    for(unsigned int sdet=0;sdet<trkHits.size(); sdet++){
      if(trkHits[sdet]->getType() == 0){   //it is very suspicious condition...
	//check whether this hit is in TPC
	float x=trkHits[sdet]->getPosition()[0];
	float y=trkHits[sdet]->getPosition()[1];
	if(std::sqrt(x*x+y*y)>=_TPC_inner){
	  //add this hit
	  tpcHitVec.push_back(trkHits[sdet]);
	}
      }
    }

    float* dummy=CalculateEnergyLoss(tpcHitVec, trkCand);
    dedx=dummy[0];
    unum=dummy[1]*_energyLossErrorTPC;

    //fill values
    trkCand->setdEdx(dedx);
    trkCand->setdEdxError(unum);
  }
}

void Compute_dEdxProcessor::check( LCEvent * ) { 
}

void Compute_dEdxProcessor::end() { 
  //delete _mydEdx;
}
 
float* Compute_dEdxProcessor::CalculateEnergyLoss(TrackerHitVec& hitVec, Track* trk){
  int tH=hitVec.size();

  //estimate truncated point(s) - (largest energy point(s)?)
  int nTruncate=0;
  //first, eatimate max. and min. enegy which shold be truncated
  //making deposit energy vectors
  float *depe= new float[tH];
  //float *depe2= new float[tH];
  int itpc=0; 
  double dx=0.0,dy=0.0,dz=0.0,length=0.0,totlength=0.0;
  for(int iH=1;iH<tH;iH++){
    TrackerHit * hit = hitVec[iH];
    //std::cout << "hit type: " << hit->getType() << std::endl;
    //if( BitSet32( hit->getType() )[  lcio::ILDDetID::TPC   ] )  {
    if( true )  {
      dx=hit->getPosition()[0]-hitVec[iH-1]->getPosition()[0];
      dy=hit->getPosition()[1]-hitVec[iH-1]->getPosition()[1];
      dz=hit->getPosition()[2]-hitVec[iH-1]->getPosition()[2];

      //cal. length of curvature
      double dxy=sqrt(dx*dx+dy*dy);
      double r=1.0/trk->getOmega();
      double theta=2.0*asin(0.5*dxy/r);
      double l=r*theta;
      //total helix length
      length=sqrt(l*l+dz*dz);

      //sum of flight length
      totlength += length;

      depe[itpc]=hit->getEDep()/length;  //20150421   changed to the distance between 2 hits
      //depe2[iH]=depe[iH];
      itpc++;
    }
  }

  //second, get truncated points vector
  //changing vectors in descending order
  float tmpe=0.0;
  for(int iH=0;iH<itpc;iH++){
    for(int jH=iH+1;jH<itpc;jH++){
      if(depe[iH]<depe[jH]){
        tmpe=depe[iH];
        depe[iH]=depe[jH];
        depe[jH]=tmpe;
      }
    }
  }

  //check max. and min energy to be truncated(min 8% and max 30%) (aleph: min 8% max 40%)
  int maxi=std::max((int)(0.30*itpc)-1,0);
  int mini=std::min((int)(0.92*itpc)-1,itpc-1);
  float maxe=depe[maxi];
  float mine=depe[mini];
  //float unum=(float)(mini-maxi+1);

  //std::cout << "check max and min: " << maxe << " " << mine << std::endl;
  //
  //    Reset
  //
  double dedx = 0.0;
  //int biH=-1;
  for(int iH=0;iH<itpc;iH++){
    //TrackerHitExtended * hitExt = hitVec[iH];
    //TrackerHit * hit = hitExt->getTrackerHit();
    if(float(depe[iH])>mine && float(depe[iH])<maxe){   //need to check here
      //if(float(depe[iH])<maxe){   //need to check here
      dedx+=depe[iH];
      nTruncate++;
      //biH=iH;
    }
  }

  //cal. truncated mean
  dedx = dedx/(float)(nTruncate);
  //std::cout << "check dedx: " << dedx << " " << unum << std::endl;

  double trkcos=sqrt(trk->getTanLambda()*trk->getTanLambda()/(1.0+trk->getTanLambda()*trk->getTanLambda()));
  float normdedx=getNormalization(dedx, (float) nTruncate, trkcos);
  if(_isSmearing) normdedx += getSmearing(normdedx);   //additional smearing 20170901

  float *ret=new float[2];
  ret[0]=normdedx;
  ret[1]=normdedx*std::pow(0.001*totlength, -0.34)*std::pow(nTruncate, -0.45);   //todo: what is gas pressure in TPC?

  if(ret[0]!=ret[0]){
    ret[0]=0.0;
    ret[1]=0.0;
  }

  return ret;
}

//correct polar angle dependence and number of hits dependence
float Compute_dEdxProcessor::getNormalization(double dedx, float hit, double trkcos){
  //cal. hit dep.
  double f1=1.0+std::exp(-hit/_ncorrpar);
  //cal. polar angle dep.
  // double c=1.0/sqrt(1.0-trkcos*trkcos);
  // double f2=1.0/(1.0-0.08887*std::log(c)); 

  //cal. polar angle dep.   20160702
  //double c = std::acos(trkcos);
  //if(c>3.141592/2.0) c= 3.141592-c;
  //double f2 = 1.0/std::pow(c, 0.0703);

  //change polar angle dep. 20170901
  double f2 = _acorrpar[0] / (_acorrpar[0] + _acorrpar[1] * trkcos * trkcos);
   
  return dedx/(f1*f2);
}

//create smearing factor
float Compute_dEdxProcessor::getSmearing(float dEdx){
  double z=sqrt( -2.0*std::log(dist(*engine)) ) * std::sin( 2.0*3.141592*dist(*engine) );

  //std::cout << dist(*engine) << " " << z << std::endl;

  return 0.0 + dEdx * _smearingFactor * z;
}
