/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "VTXDigiProcessor.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>
#include <EVENT/MCParticle.h>
// #include "random.h"
// #include <CLHEP/Random/RandGauss.h>
#include <gsl/gsl_randist.h>
#include "marlin/ProcessorEventSeeder.h"


#include "CLHEP/Vector/TwoVector.h"

#include <cmath>
#include <algorithm>
#include <sstream>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

VTXDigiProcessor aVTXDigiProcessor ;


VTXDigiProcessor::VTXDigiProcessor() : Processor("VTXDigiProcessor") {
  
  // modify processor description
  _description = "VTXDigiProcessor should create VTX TrackerHits from SimTrackerHits" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "SmearAlongLadders" ,
                              "Points smeared along the ladders"  ,
                              _smearAlongLadders ,
                              int(1)) ;

  registerProcessorParameter( "PointResolutionRPhi_VTX" ,
                              "R-Phi Resolution in VTX"  ,
                              _pointResoRPhi_VTX ,
                              float(0.0040)) ;
	
  registerProcessorParameter( "PointResolutionZ_VTX" , 
                              "Z Resolution in VTX" ,
                              _pointResoZ_VTX ,
                              float(0.0040));

  registerProcessorParameter( "PointResolutionRPhi_SIT" ,
                              "R-Phi Resolution in SIT"  ,
                              _pointResoRPhi_SIT ,
                              float(0.010)) ;
	
  registerProcessorParameter( "PointResolutionZ_SIT" , 
                              "Z Resolution in SIT" ,
                              _pointResoZ_SIT ,
                              float(0.010));

  registerProcessorParameter( "PointResolutionRPhi_SET" ,
                              "R-Phi Resolution in SET"  ,
                              _pointResoRPhi_SET ,
                              float(0.010)) ;
	
  registerProcessorParameter( "PointResolutionZ_SET" , 
                              "Z Resolution in SET" ,
                              _pointResoZ_SET ,
                              float(0.010));

  std::vector<int> activeSETLayers ;
  activeSETLayers.push_back( 1 ) ;
  
  registerProcessorParameter( "ActiveSETLayers",
                              "only SET hits from active layers are digitized (mimic stereo layers)",
                              _activeSETLayers,
                              activeSETLayers );
  
  registerProcessorParameter( "RemoveDrays" ,
                              "Remove D-rays ?",
                              _removeDRays,
                              int(0));

  registerProcessorParameter( "MomentumCutForDRays" , 
                              "Momentum Cut For D Rays (MeV)",
                              _momCut ,
                              float(10.0));

  registerProcessorParameter( "Debug",
                              "Debugging option",
                              _debug,
                              int(0)); 
 

  FloatVec effDefault(6) ;
  for(int i=0 ;i<6; i++ ) 
    effDefault[i] = 1.0 ;
  
  
  registerProcessorParameter( "HitEfficiencyPerLayer_VTX" ,
                              "hit efficiencies per VXD layer (default: 1.0)"  ,
                              _vxdEff ,
                              effDefault ) ;


  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "VTXCollectionName" , 
                           "Name of the VTX SimTrackerHit collection"  ,
                           _colNameVTX ,
                           std::string("VXDCollection") ) ;
  
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "SITCollectionName" , 
                           "Name of the SIT SimTrackerHit collection"  ,
                           _colNameSIT ,
                           std::string("SITCollection") ) ;
  
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "SETCollectionName" , 
                           "Name of the SET SimTrackerHit collection"  ,
                           _colNameSET ,
                           std::string("SETCollection") ) ;
  
  // Output collections
  registerOutputCollection( LCIO::TRACKERHIT,
                            "VTXHitCollection" , 
                            "Name of the vxd TrackerHit output collection"  ,
                            _outColNameVTX ,
                            std::string("VTXTrackerHits") ) ;

  registerOutputCollection( LCIO::TRACKERHIT,
                            "SITHitCollection" , 
                            "Name of the sit TrackerHit output collection"  ,
                            _outColNameSIT ,
                            std::string("SITTrackerHits") ) ;
  
  registerOutputCollection( LCIO::TRACKERHIT,
                            "SETHitCollection" , 
                            "Name of the set TrackerHit output collection"  ,
                            _outColNameSET ,
                            std::string("SETTrackerHits") ) ;
  
}


void VTXDigiProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);
  Global::EVENTSEEDER->registerProcessor(this);

  // check that we have the efficiencies for all layers
  const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters() ;
  const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout(); 
  const unsigned int nLayer = layerVXD.getNLayers() ;
  
  if( _vxdEff.size() < nLayer ){

    std::stringstream s ;
    s << " VTXDigiProcessor::init - wrong number of VXD efficiencies given : " << _vxdEff.size()
      << " expected one for every " <<  nLayer  << " layers " << std::endl ;
    throw Exception( s.str() ) ;
    
  }

  _vxdCount.resize( nLayer ) ;
}

void VTXDigiProcessor::processRunHeader( LCRunHeader*  /*run*/) { 
  _nRun++ ;
} 

void VTXDigiProcessor::processEvent( LCEvent * evt ) { 

  gsl_rng_set( _rng, Global::EVENTSEEDER->getSeed(this) ) ;   
  streamlog_out( DEBUG ) << "seed set to " << Global::EVENTSEEDER->getSeed(this) << std::endl;

  for (int iColl=0;iColl<3;++iColl) {

    LCCollection* STHcol = 0 ;
    try{
      if (iColl==0)
        STHcol = evt->getCollection( _colNameVTX ) ;
      else if (iColl==1)
        STHcol = evt->getCollection( _colNameSIT ) ;
      else 
        STHcol = evt->getCollection( _colNameSET ) ;
    }
    catch(DataNotAvailableException &e){

      if (iColl==0){
        streamlog_out(DEBUG) << "Collection " << _colNameVTX.c_str() << " is unavailable in event " << _nEvt << std::endl;
      }
      else if (iColl==1){
        streamlog_out(DEBUG) << "Collection " << _colNameSIT.c_str() << " is unavailable in event " << _nEvt << std::endl;
      }
      else{
        streamlog_out(DEBUG) << "Collection " << _colNameSET.c_str() << " is unavailable in event " << _nEvt << std::endl; 
      }
    }

    if( STHcol != 0 ){    
    
      LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;
    
      int nSimHits = STHcol->getNumberOfElements()  ;
    
      for(int i=0; i< nSimHits; i++){
        
        SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( STHcol->getElementAt( i ) ) ;
        
        //fg: --- accept only activeSETLayers :
        if( iColl==2 ) { // SET 
          
          // fg: the layer number is the cellID (for Mokka at least) 
          int cellID = SimTHit->getCellID0() ;
          
          if( find(_activeSETLayers.begin(),_activeSETLayers.end(),cellID)==_activeSETLayers.end() ) {
            continue ;   // ----------------- ignore hit 
          } 
        }
        
        if (_removeDRays) { // check if hit originates from delta-electron 
          float totMomentum = 0;
          for (int j=0;j<3;++j) 
            {
              totMomentum+=SimTHit->getMomentum()[j]*SimTHit->getMomentum()[j];
            }
          totMomentum = sqrt(totMomentum);
          
          if (totMomentum < _momCut)
            {
              streamlog_out( DEBUG ) << " removeDRays enabled: hit originates from delta-electron, hit dropped" << std::endl ;
              continue ;  // ----------------- ignore hit 
            }
        }
        

        const int celId = SimTHit->getCellID0() ;
        
        const double *pos ;
        pos =  SimTHit->getPosition() ;  
        
        double smearedPos[3];
        
        //VXD smearing
        if (iColl==0) {        

          streamlog_out( DEBUG ) << " processing collection " << _colNameVTX 
                                 << " with " <<  nSimHits  << " hits ... " << std::endl ;
            
          //find which layer hit is in - encoded in cell ID
          int layer = SimTHit->getCellID0() - 1;


          _vxdCount[layer].second++ ;
          // drop hit due to inefficiency ?
          double urand = gsl_rng_uniform( _rng ) ;

          if( urand > _vxdEff[  layer ] ){
            streamlog_out( DEBUG ) << " dropping hit in layer " << layer << std::endl ;
            continue ;  // ----------------- ignore hit 
          } 
          _vxdCount[layer].first++ ;


          //get VXD geometry info
          const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters() ;
          const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout(); 

          gear::Vector3D hitvec(pos[0],pos[1],pos[2]);
          gear::Vector3D smearedhitvec(pos[0],pos[1],pos[2]);

          streamlog_out(DEBUG) <<"Position of hit before smearing = "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
            
          //check detector geometry to decide how to smear hits    
          //ladders in layer -> need to smear hits along ladder plane
          
          if( _smearAlongLadders !=0 && layerVXD.getNLadders(0) !=0 ) {

            streamlog_out(DEBUG) << "start smearing along ladders for: " << layer << std::endl;
            
            int hitLayer = SimTHit->getCellID0() - 1;
              
            //phi between each ladder
            double deltaPhi = ( 2 * M_PI ) / layerVXD.getNLadders(hitLayer) ;
              
            double PhiInLocal(0.0);
            //find the ladder that the hit is in
            int ladderIndex = -1;
            double ladderPhi=999;
              
            for (int ic=0; ic < layerVXD.getNLadders(hitLayer); ++ic) {
                
              ladderPhi = correctPhiRange( layerVXD.getPhi0( hitLayer ) + ic*deltaPhi ) ;
                
              PhiInLocal = hitvec.phi() - ladderPhi;
              double RXY = hitvec.rho();
                
              // check if point is in range of ladder
              if (RXY*cos(PhiInLocal) - layerVXD.getSensitiveDistance(hitLayer) > -layerVXD.getSensitiveThickness(hitLayer) && 
                  RXY*cos(PhiInLocal) - layerVXD.getSensitiveDistance(hitLayer) <  layerVXD.getSensitiveThickness(hitLayer) )
                {
                  ladderIndex = ic;
                  break;
                }
            }

            double sensitive_width  = layerVXD.getSensitiveWidth(hitLayer);
            double sensitive_offset = layerVXD.getSensitiveOffset(hitLayer);

            double ladder_incline = correctPhiRange( (M_PI/2.0 ) + ladderPhi );
              
            double u = (hitvec.rho() * sin(PhiInLocal) - sensitive_offset );
            
            streamlog_out(DEBUG) << ":" 
                                 << " Event: " << _nEvt 
                                 << " hit: " << i 
                                 << " of "   << nSimHits
                                 << "  layer: " << hitLayer 
                                 << "  ladderIndex: " << ladderIndex 
                                 << "  half ladder width " << sensitive_width * 0.5 
                                 << "  u: " <<  u
                                 << "  layer sensitive_offset " << sensitive_offset
                                 << "  layer phi0 " << layerVXD.getPhi0( hitLayer )
                                 << "  phi: " <<  hitvec.phi()
                                 << "  PhiInLocal: " << PhiInLocal
                                 << "  ladderPhi: " << ladderPhi
                                 << "  ladder_incline: " << ladder_incline
                                 << std::endl;
              
            if( u > sensitive_width * 0.5 || u < -sensitive_width * 0.5)
              {
                streamlog_out(DEBUG) << "hit not in sensitive: u: " << u << " half ladder width = " << sensitive_width * 0.5 << std::endl;
                continue; // hit is not in sensitive so drop this hit and go on to the next
              }

            int  tries = 0;              
            // try to smear the hit within the ladder
            bool accept_hit = false;
            while( tries < 100 )
              {
                  
                if(tries > 0) streamlog_out(DEBUG) << "retry smearing for " << hitLayer << " " << ladderIndex << " : retries " << tries << std::endl;
                  
                _pointResoRPhi = _pointResoRPhi_VTX;
                _pointResoZ    = _pointResoZ_VTX;

                double rPhiSmear  = gsl_ran_gaussian(_rng, _pointResoRPhi);
                  
                if( (u+rPhiSmear) < sensitive_width * 0.5 && (u+rPhiSmear) > -sensitive_width * 0.5)
                  {
                    accept_hit =true;
                    double zSmear  = gsl_ran_gaussian(_rng, _pointResoZ);
                    
                    //find smearing for x and y, so that hit is smeared along ladder plane
                    smearedPos[0] = hitvec.x() + rPhiSmear * cos(ladder_incline);
                    smearedPos[1] = hitvec.y() + rPhiSmear * sin(ladder_incline); 
                    smearedPos[2] = hitvec.z() + zSmear;
                    break;
                    
                  }
                ++tries;
              } 
            if( accept_hit == false ) 
              {
                streamlog_out(DEBUG) << "hit could not be smeared within ladder after 100 tries: hit dropped"  << std::endl;
                continue; 
              } // 
          }
            
          else  // no ladders in layers -> just smear around cylinders
            {
                
              streamlog_out(DEBUG) << "start simple smearing for: " << layer << std::endl;
              CLHEP::Hep3Vector point(pos[0], pos[1], pos[2]);
                
              _pointResoRPhi = _pointResoRPhi_VTX;
              _pointResoZ    = _pointResoZ_VTX;

              double rphiSmear = gsl_ran_gaussian(_rng, _pointResoRPhi);
              double zSmear = gsl_ran_gaussian(_rng, _pointResoZ);
              
              point.setPhi( point.phi() + rphiSmear / point.perp() );
              point.setZ( point.z() + zSmear );
              
              smearedPos[0] = point.x();
              smearedPos[1] = point.y();
              smearedPos[2] = point.z();
            }
          
        }
        // SIT/SET Smearing
        else {
            
          if (iColl==1) {
            _pointResoRPhi = _pointResoRPhi_SIT;
            _pointResoZ = _pointResoZ_SIT;               
          }
          else {
            _pointResoRPhi = _pointResoRPhi_SET;
            _pointResoZ = _pointResoZ_SET; 
          }
            
          double xSmear = gsl_ran_gaussian(_rng, _pointResoRPhi);
          double zSmear = gsl_ran_gaussian(_rng, _pointResoZ);
            
          double phi = atan2(pos[1],pos[0]);
          double rad = sqrt(pos[1]*pos[1]+pos[0]*pos[0]);
          double phi_new = phi + xSmear/rad;
          smearedPos[0] = rad*cos(phi_new);
          smearedPos[1] = rad*sin(phi_new);
          smearedPos[2] = pos[2] + zSmear;      
            
        }
          
        float edep ;
        float dedxSmear = 0.0 ;
        edep = SimTHit->getEDep() ;
          
        edep = edep + dedxSmear ; 
          
        MCParticle *mcp ;
        mcp = SimTHit->getMCParticle() ;
          
        //store hit variables
        TrackerHitImpl* trkHit = new TrackerHitImpl ;
          
        //FIXME: SJA: this is a temporary work around the set'er should take a const double * 
        trkHit->setPosition( smearedPos ) ;
          
        trkHit->setEDep( edep ) ;
        if (iColl==0) 
          trkHit->setType(100+celId ); 
        else {
          trkHit->setType(400+celId);
        }
          
        float covMat[TRKHITNCOVMATRIX]={0.,0.,_pointResoRPhi*_pointResoRPhi,0.,0.,_pointResoZ*_pointResoZ};
        trkHit->setCovMatrix(covMat);      
        // 	  push back the SimTHit for this TrackerHit
        // fg: only if we have a sim hit with proper link to MC truth
        if( mcp != 0 )  {
          trkHit->rawHits().push_back( SimTHit ) ;
        }
        //else{
        //  streamlog_out( DEBUG ) << " ignore simhit pointer as MCParticle pointer is NULL ! " << std::endl ;
        //}


        trkhitVec->addElement( trkHit ) ; 
      }      
      
    
      if (iColl==0) 
        evt->addCollection( trkhitVec , _outColNameVTX ) ;
      else if (iColl == 1)
        evt->addCollection( trkhitVec , _outColNameSIT ) ;
      else 
        evt->addCollection( trkhitVec , _outColNameSET );
    }
  }
  _nEvt ++ ;
}



void VTXDigiProcessor::check( LCEvent *  /*evt*/ ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void VTXDigiProcessor::end(){ 
  
  streamlog_out(MESSAGE) << " end()  " << name() 
                         << " processed " << _nEvt << " events in " << _nRun << " runs "
                         << std::endl ;

  streamlog_out(DEBUG) << " VXD hits - efficiency : "  << std::endl ;
  for(unsigned i=0 ; i<_vxdCount.size(); ++i) {
    
    streamlog_out(DEBUG) << "     layer " << i << " kept " 
                           << _vxdCount[i].first << " hits of " 
                           << _vxdCount[i].second << " -> eff = " 
                           <<  double(_vxdCount[i].first)/ double(_vxdCount[i].second) 
                           << std::endl ;
      
  }
}


double VTXDigiProcessor::correctPhiRange( double Phi ) const {

  while( (Phi < -1.0*M_PI) || (Phi > 1.0*M_PI) )
    {
      if( Phi > 1.0*M_PI )
        {
          Phi -= 2.0 * M_PI;
        }
      else
        {
          Phi += 2.0 * M_PI;
        }
    }
  
  return Phi ;
  
} // function correctPhiRange
