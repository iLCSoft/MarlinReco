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

#include <cmath>
#include <algorithm>
//#include <math.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

VTXDigiProcessor aVTXDigiProcessor ;


VTXDigiProcessor::VTXDigiProcessor() : Processor("VTXDigiProcessor") {
  
  // modify processor description
  _description = "VTXDigiProcessor should create VTX TrackerHits from SimTrackerHits" ;
  

  // register steering parameters: name, description, class-variable, default value

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
  r = gsl_rng_alloc(gsl_rng_ranlxs2);
}

void VTXDigiProcessor::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void VTXDigiProcessor::processEvent( LCEvent * evt ) { 

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
      if (_debug == 1) {
        if (iColl==0)
          std::cout << "Collection " << _colNameVTX.c_str() << " is unavailable in event " << _nEvt << std::endl;
        else if (iColl==1)
          std::cout << "Collection " << _colNameSIT.c_str() << " is unavailable in event " << _nEvt << std::endl;
        else 
          std::cout << "Collection " << _colNameSET.c_str() << " is unavailable in event " << _nEvt << std::endl;
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
          int cellID = SimTHit->getCellID() ;
          
          if( find(_activeSETLayers.begin(),_activeSETLayers.end(),cellID)==_activeSETLayers.end() ) {
            
            continue ;   // ----------------- ignore hit 
          } 
        }
        
        bool accept = 1;
        if (_removeDRays) { // check if hit originates from delta-electron 
          float totMomentum = 0;
          for (int i=0;i<3;++i) 
            totMomentum+=SimTHit->getMomentum()[i]*SimTHit->getMomentum()[i];
          totMomentum = sqrt(totMomentum);
          if (totMomentum < _momCut)
            accept = 0;
        }
        
        if (accept == 1) {
          
          const int celId = SimTHit->getCellID() ;
          
          const double *pos ;
          pos =  SimTHit->getPosition() ;  

          double newPos[3];
          
          //VXD smearing
          if (iColl==0) {        
          
            //get VXD geometry info
            const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters() ;
            const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout(); 

            gear::Vector3D hitvec(pos[0],pos[1],pos[2]);

            if(_debug ==1)
              {
                cout<<"Position of hit before smearing = "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
                cout<<"Hit in sensitive volume? "<<gearVXD.isPointInSensitive(hitvec)<<endl;
              }

            //check detector geometry to decide how to smear hits    
            //ladders in layer -> need to smear hits along ladder plane
            bool UseLadders=true;
            if(layerVXD.getNLadders(0) !=0){
              
              //find which layer hit is in - encoded in cell ID
              int layer = SimTHit->getCellID() - 1;
              
              //check that this is a valid layer
              if (layer < 0 || layer > layerVXD.getNLayers() ) 
                {
                  if(_debug == 1)
                    cout<<"layer of hit not found, smearing in cylinder "<<endl;
                  UseLadders=false;
                }
              //if the layer is valid, smear along the ladder
              if(UseLadders){
                float PhiInLocal = 0;
 
                //phi between each ladder
                double deltaPhi = ( 2 * M_PI ) / layerVXD.getNLadders(layer) ;
                

                // get phi of hit in projection 2D
                double pPhi = getPhiPoint( hitvec ) ;

                //find the ladder that the hit is in
                int ladderIndex = -1;
                double lPhi=999;
                for (int ic=0; ic < layerVXD.getNLadders(layer); ++ic) 
                  {
                    lPhi = layerVXD.getPhi0( layer ) + ic*deltaPhi ;
                    lPhi = correctPhiRange( lPhi ) ;

                    PhiInLocal = pPhi - lPhi;
                    float RXY = sqrt((pos[0]*pos[0]+pos[1]*pos[1]));

                    // check if point is in range of ladder
                    if (RXY*cos(PhiInLocal)- layerVXD.getSensitiveDistance(layer) > -layerVXD.getSensitiveThickness(layer) && 
                        RXY*cos(PhiInLocal)-layerVXD.getSensitiveDistance(layer) < layerVXD.getSensitiveThickness(layer) )
                      {
                        ladderIndex = ic;
                        break;
                      }
                  }
                _pointResoRPhi = _pointResoRPhi_VTX/cos(PhiInLocal);
                _pointResoZ    = _pointResoZ_VTX;

                //finding the smearing constant
                double rSmear  = gsl_ran_gaussian(r,_pointResoRPhi_VTX);
                
                //find smearing for x and y, so that hit is smeared along ladder plane
                double xSmear = rSmear * cos(lPhi);
                double ySmear = rSmear * sin(lPhi);

                newPos[0] = pos[0] + xSmear;
                newPos[1] = pos[1] - ySmear;

                //smear in z
                double zSmear = gsl_ran_gaussian(r,_pointResoZ_VTX);

                newPos[2] = pos[2] + zSmear;

                //check that hits are still on ladder, if they aren't put them on end of ladder. 
                gear::Vector3D smearedhitvec(newPos[0],newPos[1],newPos[2]);

                //info for debug
                if(_debug==1)
                  {
                    cout<<"Position of hit after smearing = "<<newPos[0]<<" "<<newPos[1]<<" "<<newPos[2]<<endl;
                    cout<<"Smeared hits on sensitive volume? "<<gearVXD.isPointInSensitive(smearedhitvec)<<endl;
                  }

                if(gearVXD.isPointInSensitive(smearedhitvec)==0)
                  {
                    
                    //start phi for first ladder in layer
                    double startPhi = layerVXD.getPhi0(layer) + 
                      atan( (-layerVXD.getSensitiveWidth(layer) /2 + layerVXD.getSensitiveOffset(layer)) / (layerVXD.getSensitiveDistance(layer)));
                    //end phi for first ladder in layer
                    double endPhi = layerVXD.getPhi0(layer) +  
                      atan( (layerVXD.getSensitiveWidth(layer) /2 + layerVXD.getSensitiveOffset(layer)) / (layerVXD.getSensitiveDistance(layer)));

                    // get start and end phi for the ladder that this hit is on 
                    float sPhi = correctPhiRange( startPhi + ladderIndex*deltaPhi ) ;
                    float ePhi = correctPhiRange( endPhi + ladderIndex*deltaPhi) ;

                    pPhi = getPhiPoint( smearedhitvec ) ;

                    //point in ladder where line perpendicular to ladder, from interaction point crosses ladder in x and y
                    float xladder = sin(lPhi)*(layerVXD.getSensitiveDistance(layer)+(0.5*layerVXD.getSensitiveThickness(layer)));
                    float yladder = cos(lPhi)*(layerVXD.getSensitiveDistance(layer)+(0.5*layerVXD.getSensitiveThickness(layer)));

                    //end of ladder in z
                    float endz = layerVXD.getSensitiveLength(layer);

                    //if point has been smeared further than start of ladder move it to start of ladder in x and y
                    if(pPhi < sPhi)
                      {
                        newPos[0] = xladder - (cos(lPhi)*((layerVXD.getSensitiveWidth(layer)*0.5)-layerVXD.getSensitiveOffset( layer )));
                        newPos[1] = yladder + (sin(lPhi)*((layerVXD.getSensitiveWidth(layer)*0.5)-layerVXD.getSensitiveOffset( layer )));
                      }
                    //if point has been smeared further than end of ladder move it to end of ladder in x and y
                    if(pPhi > ePhi)
                      {
                         newPos[0] = xladder + (cos(lPhi)*((layerVXD.getSensitiveWidth(layer)*0.5)+layerVXD.getSensitiveOffset( layer )));
                         newPos[1] = yladder - (sin(lPhi)*((layerVXD.getSensitiveWidth(layer)*0.5)+layerVXD.getSensitiveOffset( layer )));
                       }
                    
                  
                    //if point has been smeared further than ends of ladder in z, move it to the nearest end.
                    if(newPos[2] > endz)
                      newPos[2] = endz;
                    if(newPos[2] < -endz)
                      newPos[2] = -endz;
                
                    //info for debug
                    gear::Vector3D repohitvec(newPos[0],newPos[1],newPos[2]);
                    if(_debug==1)
                      {
                        cout<<"Position of hit after repositioning = "<<newPos[0]<<" "<<newPos[1]<<" "<<newPos[2]<<endl;
                        cout<<"repositioned hits on sensitive volume? "<<gearVXD.isPointInSensitive(repohitvec)<<endl;
                      }

                  }
              }
              
            }
     

            // no ladders in layers -> just smear around cylinders
            if(layerVXD.getNLadders(0) ==0 || UseLadders==false){
              
              _pointResoRPhi = _pointResoRPhi_VTX;
              _pointResoZ = _pointResoZ_VTX;

              double xSmear = gsl_ran_gaussian(r,_pointResoRPhi_VTX);
              double zSmear = gsl_ran_gaussian(r,_pointResoZ_VTX);
              
              double phi = atan2(pos[1],pos[0]);
              double rad = sqrt(pos[1]*pos[1]+pos[0]*pos[0]);
              double phi_new = phi + xSmear/rad;
              newPos[0] = rad*cos(phi_new);
              newPos[1] = rad*sin(phi_new);
              newPos[2] = pos[2] + zSmear;
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

            double xSmear = gsl_ran_gaussian(r,_pointResoRPhi);
            double zSmear = gsl_ran_gaussian(r,_pointResoZ);

            double phi = atan2(pos[1],pos[0]);
            double rad = sqrt(pos[1]*pos[1]+pos[0]*pos[0]);
            double phi_new = phi + xSmear/rad;
            newPos[0] = rad*cos(phi_new);
            newPos[1] = rad*sin(phi_new);
            newPos[2] = pos[2] + zSmear;      

          }
        
          float de_dx ;
          float dedxSmear = 0.0 ;
          de_dx = SimTHit->getdEdx() ;
        
          de_dx = de_dx + dedxSmear ; 
        
          MCParticle *mcp ;
          mcp = SimTHit->getMCParticle() ;
          
          //store hit variables
          TrackerHitImpl* trkHit = new TrackerHitImpl ;
          
          //FIXME: SJA: this is a temporary work around the set'er should take a const double * 
          trkHit->setPosition( newPos ) ;
          
          trkHit->setdEdx( de_dx ) ;
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



  void VTXDigiProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void VTXDigiProcessor::end(){ 
  
  std::cout << "VTXDigiProcessor::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    << std::endl ;
}


double VTXDigiProcessor::correctPhiRange( double Phi ) const {
    
    if( Phi > M_PI ) {
      return ( Phi - 2 * M_PI ) ;
    }
    if( Phi <= -M_PI ) {
      return ( Phi + 2 * M_PI ) ;
    } 
    
    return Phi ;

  } // function correctPhiRange


double VTXDigiProcessor::getPhiPoint( gear::Vector3D p ) const {

    //fg: definition of phi - seems like this is the the angle with the negative y-axis ????
    //    return correctPhiRange( p.phi() ) ;

    // get phi of point in projection 2D
    double pPhi = 0. ;
    if( ( p[0] >= 0 ) && ( p[1] == 0 ) )
      pPhi = -M_PI/2 ;
    
    if( ( p[0] < 0 ) && ( p[1] == 0 ) )
      pPhi = M_PI/2 ;
    
    if( ( p[0] == 0 ) && ( p[1] < 0 ) ) 
      pPhi = 0 ;
    
    if( ( p[0] != 0 ) && ( p[1] < 0 ) )
      pPhi = atan( p[0] / p[1] ) + M_PI ;
    
    else
      pPhi = atan( p[0] / p[1] ) ;
    
    pPhi = correctPhiRange( pPhi ) ;  
    
    return pPhi ;
    
} // function getPhiPoint
