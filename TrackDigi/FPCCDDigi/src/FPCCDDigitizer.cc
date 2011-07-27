/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "FPCCDDigitizer.h"
#include "FPCCDPixelHit.h"
#include "FPCCDData.h"

#include <iostream>
#include <cstring>
#include <cstdio>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>
#include "UTIL/LCTOOLS.h"

#include <cmath>
#include <algorithm>
#include <sstream>
#include <map>
#include <TRandom.h>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

namespace
{
  template<class T>
  class SortByPosX
  {
  public:
    bool operator()( const T& object1, const T& object2 )
    {
      return object1->x() < object2->x();
    }
  };

}


FPCCDDigitizer aFPCCDDigitizer ;



// =====================================================================
FPCCDDigitizer::FPCCDDigitizer() : Processor("FPCCDDigitizer") {
  
  // modify processor description
  _description = "FPCCDDigitizer should create FPCCD's VTXPixelHits from SimTrackerHits" ;
  
  // register steering parameters: name, description, class-variable, default value

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

  registerProcessorParameter( "Debug",
                              "Debugging option",
                              _debug,
                              int(0)); 

  registerProcessorParameter( "modifySimTrackerHit",
                              "Modify SimTrackerHit into layer",
                              _modifySimTHit,
                              bool(true)); 

  registerProcessorParameter( "FPCCD_pixelSize(mm)" ,
                              "Pixel size of FPCCD (unit:mm) (default: 0.005)"  ,
                              _pixelSize ,
                              float(0.005) ) ;
  
  registerProcessorParameter( "PixelHeight(mm)" , 
                              "Pixel Height(mm)",
                              _pixelheight ,
                              float(0.015));

  registerProcessorParameter("sigmaConstant",
                             "ration of sigma of Landau distribution and path length",
                             _sigmaConst,
                             double(1600/0.05));

  registerProcessorParameter("HitQuality",
                             "Hit Quality of the Event",
                             _isSignal,
                             bool(true));

  registerProcessorParameter( "HelixTMomentumCriteria(MeV)" ,
                              "Transverse momentum criteria for helix approximation (MeV)",
                              _momCut ,
                              float(100.0));

  registerProcessorParameter( "Ladder_Number_encoded_in_cellID" ,
                              "Mokka has encoded the ladder number in the cellID" ,
                              _ladder_Number_encoded_in_cellID ,
                              bool(true));

  registerProcessorParameter( "VTX_Only" ,
                              "Use digitizer for SIT and SET another",
                              _vtx_only,
                              bool(false));
  
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
  registerOutputCollection( LCIO::LCGENERICOBJECT,
                            "VTXPixelHitCollection" , 
                            "Name of the VXD PixelHit output collection"  ,
                            _outColNameVTX ,
                            std::string("VTXPixelHits") ) ;
  
  registerOutputCollection( LCIO::TRACKERHIT,
                            "SITHitCollection" , 
                            "Name of the SIT TrackerHit output collection"  ,
                            _outColNameSIT ,
                            std::string("SITTrackerHits") ) ;
  
  registerOutputCollection( LCIO::TRACKERHIT,
                            "SETHitCollection" , 
                            "Name of the SET TrackerHit output collection"  ,
                            _outColNameSET ,
                            std::string("SETTrackerHits") ) ;
  
  
}


// =====================================================================
void FPCCDDigitizer::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  InitGeometry();
}

// =====================================================================
void FPCCDDigitizer::InitGeometry() 
{ 
  // Save frequently used parameters.

//   const gear::VXDParameters &gearVXD = Global::GEAR->getVXDParameters();
//   const gear::VXDLayerLayout &layerVXD = gearVXD.getVXDLayerLayout();
  const gear::VXDParameters &gearVXD = Global::GEAR->getVXDParameters();
  const gear::VXDLayerLayout &layerVXD = gearVXD.getVXDLayerLayout();
  _nLayer = layerVXD.getNLayers() ;
  _geodata.resize(_nLayer);
  _maxLadder = 0;
  
  for(int ly=0;ly<_nLayer;ly++){
    _geodata[ly].nladder = layerVXD.getNLadders(ly);           // Number of ladders in this layer
    if( _maxLadder < _geodata[ly].nladder ) { _maxLadder = _geodata[ly].nladder; }
    _geodata[ly].rmin = layerVXD.getSensitiveDistance(ly);     // Distance of sensitive area from
    _geodata[ly].dphi = (2*M_PI)/(double)_geodata[ly].nladder;
    _geodata[ly].phi0 = layerVXD.getPhi0(ly);                  // phi offset.
    _geodata[ly].sthick = layerVXD.getSensitiveThickness(ly);
    _geodata[ly].sximin = -layerVXD.getSensitiveOffset(ly)
                          -layerVXD.getSensitiveWidth(ly)/2.0;
    _geodata[ly].sximax = -layerVXD.getSensitiveOffset(ly)
                          +layerVXD.getSensitiveWidth(ly)/2.0;
    _geodata[ly].hlength = layerVXD.getSensitiveLength(ly);
    _geodata[ly].cosphi.resize( _geodata[ly].nladder );
    _geodata[ly].sinphi.resize( _geodata[ly].nladder );
    for(int ld = 0;ld<_geodata[ly].nladder;ld++) {
      double phi = _geodata[ly].phi0 + _geodata[ly].dphi*ld;
      _geodata[ly].cosphi[ld] = cos(phi);
      _geodata[ly].sinphi[ld] = sin(phi);
    }
  }
} 


// =====================================================================
void FPCCDDigitizer::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

// =====================================================================
void FPCCDDigitizer::modifyEvent( LCEvent * evt )
{

  int nCol = 3;
  if(_vtx_only) nCol = 1;
  
  for (int iColl=0; iColl<nCol; ++iColl) {

    LCCollection* col = 0 ;

    try{
      if (iColl==0)
        col = evt->getCollection( _colNameVTX ) ;
      else if (iColl==1)
        col = evt->getCollection( _colNameSIT ) ;
      else 
        col = evt->getCollection( _colNameSET ) ;
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

    if( col != 0 ){    
      
      LCCollectionVec* trkhitVec = new LCCollectionVec( LCIO::TRACKERHIT )  ;

      if(iColl==0){
      
        LCCollectionVec* fpccdDataVec = new LCCollectionVec( LCIO::LCGENERICOBJECT )  ;
        
        FPCCDData theData(_nLayer, _maxLadder);  // prepare object to make pixelhits
        
        int nSimHits = col->getNumberOfElements()  ;
        
        // Loop over all VXD hits
        
        for(int i=0; i< nSimHits; i++){
          SimTrackerHitImpl* SimTHitImpl = dynamic_cast<SimTrackerHitImpl*>( col->getElementAt( i ) ) ;
          
          makePixelHits(SimTHitImpl,  theData);        

        } // End of loop over all SimTrackerHits

        theData.packPixelHits( *fpccdDataVec );

        evt->addCollection( fpccdDataVec , _outColNameVTX ) ;
        if(_debug == 1){
        std::cout << "dumpEvent Digitizer" << std::endl;
        LCTOOLS::dumpEvent( evt ) ;
        }
      }
      else{
        
        if( iColl == 1 ){
          _pointResoRPhi = _pointResoRPhi_SIT;
          _pointResoZ = _pointResoZ_SIT;               
        }
        else{
          _pointResoRPhi = _pointResoRPhi_SET;
          _pointResoZ = _pointResoZ_SET;               
        }
        for(int j = 0; j<col->getNumberOfElements(); j++){
          SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( col->getElementAt( j ) );
          
          TrackerHitImpl* trkHit = new TrackerHitImpl ;
          
          trkHit->setPosition( (double*)SimTHit->getPosition() );
          trkHit->setEDep( SimTHit->getEDep() );
          trkHit->setType( 400 + SimTHit->getCellID0() );
          float covMat[TRKHITNCOVMATRIX] = {0.,0.,_pointResoRPhi*_pointResoRPhi,0.,0.,_pointResoZ*_pointResoZ};
          trkHit->setCovMatrix( covMat );
          
          trkhitVec->addElement( trkHit );
        }
      
        if( iColl == 1)
          evt->addCollection( trkhitVec , _outColNameSIT ) ;
        else  if( iColl == 2)
          evt->addCollection( trkhitVec , _outColNameSET );
        
    }
      } // End of process when VXD has hits
  }  
  _nEvt ++ ;
}

// =====================================================================
// void FPCCDDigitizer::makePixelHits(const SimTrackerHit *simHit,
//                                    FPCCDData &hitVec)
// {
//   // Obtain ladderID, xiID and zetaID and fill them into FPCCDPixelHits,
//   // which is an array of FPCCDLadder_t

//   int layer = simHit->getCellID0() - 1;
//   const double *pos =  simHit->getPosition() ;  
//   gear::Vector3D posvec(pos[0],pos[1],pos[2]);

//   // Convert SimTrackerHit position into pixelID;

//   FPCCDID_t fpccdID=encodeFPCCDID(layer, posvec);
//   if( fpccdID.layer >= 0 ) {  // When valid fpccdid is obtained, ..
    
//     float de_dx ;
//     float dedxSmear = 0.0 ;
//     de_dx = simHit->getEDep() ;

//     MCParticle *mcp ;
//     mcp = simHit->getMCParticle() ;
//     FPCCDPixelHit::HitQuality_t quality=FPCCDPixelHit::kSingle;

//     FPCCDPixelHit aHit(fpccdID.layer, fpccdID.ladder, fpccdID.xi, fpccdID.zeta,
//                        de_dx, quality, mcp);
    
//     hitVec.addPixelHit(aHit, true);

//   }

// }

void FPCCDDigitizer::makePixelHits(SimTrackerHitImpl *SimTHit,  FPCCDData &hitVec)
{
  gear::Vector3D* HitPosInMokka = new gear::Vector3D(SimTHit->getPosition()[0],SimTHit->getPosition()[1],SimTHit->getPosition()[2]);

  //----
  // get basic info.
  //----
  const gear::BField& gearBField = Global::GEAR->getBField();
  double posphi = HitPosInMokka->phi();
  if(HitPosInMokka->y()<0) posphi += 2*M_PI;

  int layer = 0 ;
  int ladderID = 0 ;
  const int celId = SimTHit->getCellID0();
  
  if(_ladder_Number_encoded_in_cellID) {
    layer  = ( celId / 10000 ) - 1 ;
  }
  else{
    layer = celId  - 1 ;
  }

  //----
  // check which ladder hit on
  //----

  if(_ladder_Number_encoded_in_cellID) {
    ladderID = ( celId % ( 10000 * (layer + 1) ) ) -1 ;
  }
  else{
    
    ladderID = getLadderID( HitPosInMokka, layer);    
  }
    
  double  sximin = _geodata[layer].sximin;
  double hlength = _geodata[layer].hlength;
  
  //----
  // get hit dir(mom) at hit points and other info.
  //----
  
  gear::Vector3D* MomAtHitPos = new gear::Vector3D(SimTHit->getMomentum()[0],SimTHit->getMomentum()[1],SimTHit->getMomentum()[2]);
  
  double momphi = MomAtHitPos->phi();
  if(MomAtHitPos->y()<0) momphi += 2*M_PI;
  gear::Vector3D origin;
  gear::Vector3D bfield = gearBField.at(origin);
  gear::Vector3D* BField = new gear::Vector3D(bfield.x(),bfield.y(),bfield.z());
  MCParticle* mcp = SimTHit -> getMCParticle();
  if( mcp == 0 ) return ;
  float charge = mcp->getCharge();
  
  //----
  // get local pos and dir on each ladder
  //----
  
  gear::Vector3D*      LocalHitPos = getLocalPos(HitPosInMokka, layer, ladderID);
  gear::Vector3D* MomAtLocalHitPos = new gear::Vector3D(MomAtHitPos->x()*_geodata[layer].sinphi[ladderID]-MomAtHitPos->y()*_geodata[layer].cosphi[ladderID],
                                                        MomAtHitPos->z(),
                                                        MomAtHitPos->x()*_geodata[layer].cosphi[ladderID]+MomAtHitPos->y()*_geodata[layer].sinphi[ladderID]);
  gear::Vector3D*      LocalBField = new gear::Vector3D( _geodata[layer].sinphi[ladderID]*BField->x()-_geodata[layer].cosphi[ladderID]*BField->y(),
                                                         BField->z(),
                                                         _geodata[layer].cosphi[ladderID]*BField->x()+_geodata[layer].sinphi[ladderID]*BField->y());
  gear::Vector3D* PosOutFromLadder = new gear::Vector3D(0,0,0);
  gear::Vector3D*    PosInToLadder = new gear::Vector3D(0,0,0);
  
  if( sqrt(MomAtLocalHitPos->x()*MomAtLocalHitPos->x() + MomAtLocalHitPos->z()*MomAtLocalHitPos->z()) < _momCut*1e-03){ // transvers momentum criteria for helix approximation
    //----
    // Helix Fitting
    //----             
    if( _debug == 1 ) cout << "========== Track of this hit is trated as a helix ==========" << endl;
    getInOutPosOfHelixOnLadder(SimTHit, PosOutFromLadder, PosInToLadder, LocalHitPos, MomAtLocalHitPos, LocalBField,charge);
  }
  else{
    //----
    //Line Fittng
    //----  
    getInOutPosOnLadder(SimTHit, PosOutFromLadder,PosInToLadder,LocalHitPos,MomAtLocalHitPos); 
    *LocalHitPos = gear::Vector3D((PosOutFromLadder->x()+PosInToLadder->x())/2 ,
                                  (PosOutFromLadder->y()+PosInToLadder->y())/2 ,
                                  (PosOutFromLadder->z()+PosInToLadder->z())/2);
  }

  if( inSensitiveRegion( PosOutFromLadder, layer) && inSensitiveRegion( PosInToLadder, layer) ){ // check if the particle through the sensitive region.
  
    //----
    // make new SimTrackerHit for seinsitive thickness of FPCCD 
    //----
    if(_modifySimTHit){
      
      gear::Vector3D Path(PosOutFromLadder->x()-PosInToLadder->x(),
                          PosOutFromLadder->y()-PosInToLadder->y(),
                          PosOutFromLadder->z()-PosInToLadder->z()); 
      double PathLength = Path.r();
      
      makeNewSimTHit(SimTHit, LocalHitPos, MomAtLocalHitPos, layer, ladderID, PathLength);
    }
    
    //----
    // get hit points incoming and outgoiong each pixel (@ edge of pixel)
    //----
    
    vector<gear::Vector3D*> EdgeOfPixel = getIntersectionOfTrkAndPix(layer,PosOutFromLadder,PosInToLadder);
    
    //----
    // calculate length that particle pass through pixel & central point of pixel
    //----
    
    sort( EdgeOfPixel.begin(), EdgeOfPixel.end(), SortByPosX<gear::Vector3D*>() );
    map< pair< double, double>*, double> pixel_local = getLocalPixel(SimTHit,EdgeOfPixel);
    map< pair< double, double>*, double>::iterator i_pixel_local=pixel_local.begin();

    while( i_pixel_local != pixel_local.end() ){

      double   xi = ((*i_pixel_local).first->first-sximin-_pixelSize/2);
      double zeta = ((*i_pixel_local).first->second+hlength-_pixelSize/2);
      double mergin = 1e-2;
      unsigned short int   xiid = (unsigned short int)(xi/_pixelSize + mergin);
      unsigned short int zetaid = (unsigned short int)(zeta/_pixelSize + mergin);
      double de_dx = (float)(*i_pixel_local).second; //  [mm]   
      FPCCDPixelHit::HitQuality_t quality;

      if( _isSignal && mcp->getParents().size() == 0 ){ quality =FPCCDPixelHit::kSingle; }
      else{ quality = FPCCDPixelHit::kBKG;}
        FPCCDPixelHit aHit(layer,ladderID,xiid,zetaid,de_dx,quality,mcp);
        hitVec.addPixelHit(aHit,_isSignal);
       if(_debug ==1 ) aHit.print();

        i_pixel_local++;
    }
  }
  else{
    if(_debug == 1) cout << "This particle didn't through the sensitive region." << endl;
  }

  if(_debug == 1){
    std::cout << std::endl;
    std::cout << "============================================================================" << std::endl;
    std::cout << std::endl;
  }
    
  delete HitPosInMokka;
  delete LocalHitPos;
  delete MomAtLocalHitPos;
  delete LocalBField;
  delete PosOutFromLadder;
  delete PosInToLadder;
  
  return ;
}

// =====================================================================
void FPCCDDigitizer::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}
  
// =====================================================================
void FPCCDDigitizer::end(){ 
  
  streamlog_out(MESSAGE) << " end()  " << name() 
                         << " processed " << _nEvt << " events in " << _nRun << " runs "
                         << endl ;
}




// =====================================================================
// =====================================================================
int FPCCDDigitizer::getLadderID(const gear::Vector3D* pos, const int layer){
  double layerthickness = _geodata[layer].sthick;
  double         radius = _geodata[layer].rmin+0.5*layerthickness;
  int           nladder = _geodata[layer].nladder;
  double           dphi = _geodata[layer].dphi;
  double     offset_phi = _geodata[layer].phi0;
  double         posphi = pos->phi();
  double           posR = pos->rho();

  int ladderID = 50;
  double local_phi = 50.;
  for(int j=0; j<nladder; j++){
    local_phi = posphi - dphi*j - offset_phi;
    if((posR*cos(local_phi)-radius >= -layerthickness) && (posR*cos(local_phi)-radius <= layerthickness)){
      ladderID = j;
      break;
    }
  }
  if(ladderID>100) cout << "LADDERID>100!: " << ladderID << endl;
  local_phi = posphi - dphi*(ladderID+1) - offset_phi;
  if((posR*cos(local_phi)-radius >= -layerthickness) && (posR*cos(local_phi)-radius <= layerthickness)){
    ladderID++;
  }
  if(ladderID==50){
    for(int j=0; j<nladder; j++){
      local_phi = posphi - dphi*j - offset_phi;
      if((posR*cos(local_phi)-radius >= -1.*layerthickness) && (posR*cos(local_phi)-radius <= 1.*layerthickness)){
        ladderID = j;
        break;
      }
    }
  }
  if(ladderID==50) cout << "LADDERID==50!" << endl;
  if(ladderID>100) cout << "LADDERID>100!..." << ladderID << endl;

  return ladderID;
}

// =====================================================================
gear::Vector3D* FPCCDDigitizer::getLocalPos(const gear::Vector3D* pos, const int layer, const int ladder){
  double   posR = pos->rho();
  double posphi = pos->phi();
  double radius = _geodata[layer].rmin+0.5*_pixelheight + (layer%2)*(_geodata[layer].sthick-_pixelheight);
  gear::Vector3D* LocalPos = new gear::Vector3D(posR*(cos(posphi)*_geodata[layer].sinphi[ladder]-sin(posphi)*_geodata[layer].cosphi[ladder]),
                                           pos->z(),
                                           posR*(cos(posphi)*_geodata[layer].cosphi[ladder]+sin(posphi)*_geodata[layer].sinphi[ladder])-radius);
  
  return LocalPos;
}

// =====================================================================
void FPCCDDigitizer::getInOutPosOnLadder(IMPL::SimTrackerHitImpl* simthit,gear::Vector3D* outpos, gear::Vector3D* inpos, gear::Vector3D* pos,gear::Vector3D* mom){
  int layer = simthit->getCellID0() - 1;
  double f_z = 0.5 * _pixelheight;

  double top_y = (mom->y()/mom->z())*(f_z-pos->z())+pos->y();
  double top_x = (mom->x()/mom->z())*(f_z-pos->z())+pos->x();

  double bottom_y = (mom->y()/mom->z())*(-f_z-pos->z())+pos->y();
  double bottom_x = (mom->x()/mom->z())*(-f_z-pos->z())+pos->x();

  if( (mom->x() > 0 && top_x > bottom_x) || (mom->x() < 0 && top_x < bottom_x ) ) {
    *outpos = gear::Vector3D( top_x, top_y, f_z);
    *inpos = gear::Vector3D( bottom_x, bottom_y, -f_z);
  }
  else{
    *inpos = gear::Vector3D( top_x, top_y, f_z);
    *outpos = gear::Vector3D( bottom_x, bottom_y, -f_z);
  }

  ModifyIntoLadder( inpos, layer, inpos, mom);
  ModifyIntoLadder( outpos, layer, outpos, mom);
  
  return ;
}

// =====================================================================
void FPCCDDigitizer::ModifyIntoLadder(gear::Vector3D* bemodifiedpos,const int f_layer,gear::Vector3D* pos,gear::Vector3D* mom){
  double   sximin =  _geodata[f_layer].sximin;
  double   sximax =  _geodata[f_layer].sximax;
  double szetamin = -_geodata[f_layer].hlength;
  double szetamax =  _geodata[f_layer].hlength;

  double preposx = bemodifiedpos->x();
  double preposy = bemodifiedpos->y();
  double preposz = bemodifiedpos->z();

  if(preposx < sximin){
    if(_debug ==1 ) cout << "this hit was modified Sx" << endl; 
    preposx = sximin;
    preposy = (mom->y()/mom->x())*(preposx-pos->x())+pos->y();
    preposz = (mom->z()/mom->x())*(preposx-pos->x())+pos->z();
  }else if(preposx > sximax){
    if(_debug ==1 ) cout << "this hit was modified Ex" << endl; 
    preposx = sximax;
    preposy = (mom->y()/mom->x())*(preposx-pos->x())+pos->y();
    preposz = (mom->z()/mom->x())*(preposx-pos->x())+pos->z();
  }
  if(preposy < szetamin){
    if(_debug ==1 ) cout << "this hit was modified Sy" << endl; 
    preposy = szetamin;
    preposx = (mom->x()/mom->y())*(preposy-pos->y())+pos->x();
    preposz = (mom->z()/mom->y())*(preposy-pos->y())+pos->z();
  }else if(preposy > szetamax){
    if(_debug ==1 ) cout << "this hit was modified Ey" << endl; 
    preposy = szetamax;
    preposx = (mom->x()/mom->y())*(preposy-pos->y())+pos->x();
    preposz = (mom->z()/mom->y())*(preposy-pos->y())+pos->z();
  }
  gear::Vector3D prepos(preposx,preposy,preposz);
  *bemodifiedpos = prepos;

  return ;
}

// =====================================================================
void FPCCDDigitizer::getInOutPosOfHelixOnLadder(SimTrackerHitImpl* simthit,gear::Vector3D* outpos, gear::Vector3D* inpos, gear::Vector3D* pos,gear::Vector3D* mom,gear::Vector3D* BField,float Charge){
  double out[3];
  double in[3];
  
  int layer = simthit->getCellID0() - 1;
  double   outerZ =  0.5*_pixelheight;
  double   innerZ = -0.5*_pixelheight;
  double   sximin =  _geodata[layer].sximin;
  double   sximax =  _geodata[layer].sximax;
  double szetamin = -_geodata[layer].hlength;
  double szetamax =  _geodata[layer].hlength;
  gear::Vector3D tmom(mom->x(),0,mom->z());
  //     gear::Vector3D* tmom = new gear::Vector3D((1-BField->x()/BField->r())*mom->x(),
  //                                               (1-BField->y()/BField->r())*mom->y(),
  //                                               (1-BField->z()/BField->r())*mom->z());    
  double PI = M_PI;
  
  double      radius = (tmom.r()/299.79*BField->r())*1e+01;
  double         v_y = radius*mom->y()/sqrt(mom->x()*mom->x()+mom->z()*mom->z());
  double init_pos[3] = {pos->x(),pos->y(),pos->z()};
  double init_dir[3] = {mom->x()/sqrt(mom->x()*mom->x()+mom->z()*mom->z()), 0, mom->z()/sqrt(mom->x()*mom->x()+mom->z()*mom->z())};

  double offset_pos[3];
  double offset_phi;
  double outphi;
  double inphi;
  double tmpasin = -(1/Charge)*asin(init_dir[0]/(Charge));    
  double tmpacos = (1/Charge)*acos(init_dir[2]/(Charge));
  
  if( tmpasin >= 0){
    tmpacos >= PI/2 ? offset_phi = tmpacos : offset_phi = tmpasin;
  }
  else{
    tmpacos >= PI/2 ? offset_phi =-tmpacos : offset_phi = tmpasin;
  }
  
  offset_pos[0] = init_pos[0]-radius*cos(Charge*offset_phi);
  offset_pos[1] = init_pos[1]-v_y*offset_phi;
  offset_pos[2] = init_pos[2]-radius*sin(Charge*offset_phi);
  
  double outerZ_phi = 0;
  double innerZ_phi = 0;
  //Check if the track intesects with x-y plane. 
  if( fabs((outerZ-offset_pos[2])/radius) <= 1. ) outerZ_phi = 1 / Charge*asin((outerZ-offset_pos[2])/radius) - offset_phi;
  if( fabs((innerZ-offset_pos[2])/radius) <= 1. ) innerZ_phi = 1 / Charge*asin((innerZ-offset_pos[2])/radius) - offset_phi;
  if(outerZ_phi >= innerZ_phi ){ outphi = innerZ_phi;  inphi = outerZ_phi;}
  else{ outphi = outerZ_phi; inphi = innerZ_phi;}

  int flag = 0;
  if( radius*cos(Charge*(outphi + offset_phi)) + offset_pos[0] <= sximin
      || radius*cos(Charge*(outphi + offset_phi)) + offset_pos[0] >= sximax
      || v_y*(outphi + offset_phi) + offset_pos[1] <= szetamin
      || v_y*(outphi + offset_phi) + offset_pos[1] >= szetamax ) flag = flag + 1;
  
  if( radius*cos(Charge*(inphi + offset_phi)) + offset_pos[0] <= sximin
      || radius*cos(Charge*(inphi + offset_phi)) + offset_pos[0] >= sximax
      || v_y*(inphi + offset_phi) + offset_pos[1] <= szetamin
      || v_y*(inphi + offset_phi) + offset_pos[1] >= szetamax ) flag = flag + 2;
  
  if(flag != 0){
    
    double tmpphi;
    double   sximin_phi = 1 / Charge*acos((sximin - offset_pos[0])/radius) - offset_phi;
    double   sximax_phi = 1 / Charge*acos((sximax - offset_pos[0])/radius) - offset_phi;  
    double szetamin_phi = (szetamin - offset_pos[1]) / v_y - offset_phi;
    double szetamax_phi = (szetamax - offset_pos[1]) / v_y - offset_phi;
    
    if(fabs(sximin_phi) <= fabs(sximax_phi) && fabs(sximin_phi) <= fabs(szetamin_phi) && fabs(sximin_phi) <= fabs(szetamax_phi)){
      tmpphi = sximin_phi;
      sximin_phi = 1e+8;
    }
    else if(fabs(sximax_phi) <= fabs(szetamin_phi) && fabs(sximax_phi) <= fabs(szetamax_phi)){
      tmpphi = sximax_phi;
      sximax_phi = 1e+8;
    }
    else if(fabs(szetamin_phi) <= fabs(szetamax_phi)){
      tmpphi = szetamin_phi;
      szetamin_phi = 1e+8;
    }
    else{
      tmpphi = szetamax_phi;
      szetamax_phi = 1e+8;
    }      
    if(flag == 1) outphi = tmpphi;
    if(flag == 2) inphi = tmpphi;
    if(flag == 3){
      inphi = tmpphi;
      if(fabs(sximin_phi) <= fabs(sximax_phi) && fabs(sximin_phi) <= fabs(szetamin_phi) && fabs(sximin_phi) <= fabs(szetamax_phi)) tmpphi = sximin_phi;
      else if(fabs(sximax_phi) <= fabs(szetamin_phi) && fabs(sximax_phi) <= fabs(szetamax_phi)) tmpphi = sximax_phi;
      else if(fabs(szetamin_phi) <= fabs(szetamax_phi))tmpphi = szetamin_phi;
      else  tmpphi = szetamax_phi;
      outphi = tmpphi;
    }
  }
  out[0] = radius*cos(Charge*(outphi+offset_phi)) + offset_pos[0];
  out[1] = v_y * (outphi+offset_phi)+ offset_pos[1];
  out[2] = radius*sin(Charge*(outphi+offset_phi)) + offset_pos[2];
  
  in[0] = radius*cos(Charge*(inphi+offset_phi)) + offset_pos[0];
  in[1] = v_y * (inphi+offset_phi)+ offset_pos[1];
  in[2] = radius*sin(Charge*(inphi+offset_phi)) + offset_pos[2];
  
  double center_phi = (outphi+inphi)/2;
  *pos = gear::Vector3D(radius*cos(Charge*(offset_phi+center_phi)) + offset_pos[0],
                        v_y*(offset_phi+center_phi) + offset_pos[1],
                        radius*sin(Charge*(offset_phi+center_phi)) + offset_pos[2]);
  *mom = gear::Vector3D(-radius*sin(Charge*(offset_phi+center_phi)), v_y , radius*(Charge*(offset_phi+center_phi))); 
  *outpos = gear::Vector3D(out[0], out[1], out[2]);
  *inpos = gear::Vector3D(in[0], in[1], in[2]);
  
  return ;
}

// =====================================================================
void FPCCDDigitizer::setLoopRange(int* looprange){
  _sloopx = looprange[0];
  _eloopx = looprange[1];
  _sloopy = looprange[2];
  _eloopy = looprange[3];
  return;
}

// =====================================================================
vector<gear::Vector3D*> FPCCDDigitizer::getIntersectionOfTrkAndPix(const int f_layer, gear::Vector3D* top,gear::Vector3D* bottom){
  vector<gear::Vector3D*> EdgeOfPixel;
  double sximin = _geodata[f_layer].sximin;
  double szetamin = -_geodata[f_layer].hlength;
  double small_localx;
  double large_localx;
  double small_localy;
  double large_localy;
  if(top->x()>=bottom->x()){ large_localx = top->x(); small_localx = bottom->x(); }
  else { large_localx = bottom->x(); small_localx = top->x(); }
  if(top->y()>=bottom->y()) {large_localy = top->y(); small_localy = bottom->y(); }
  else { large_localy = bottom->y(); small_localy = top->y(); }
  int nloopx = 0;
  int nloopy = 0;
  int sloopx = 0;
  int sloopy = 0;
  int eloopx = 0;
  int eloopy = 0;
  sloopx = (int)ceil((small_localx-sximin)/_pixelSize) ;
  eloopx = (int)ceil((large_localx-sximin)/_pixelSize) - 1;
  sloopy = (int)ceil((small_localy-szetamin)/_pixelSize) ;
  eloopy = (int)ceil((large_localy-szetamin)/_pixelSize) - 1;

  int looprange[4] = {sloopx,eloopx,sloopy,eloopy};
  setLoopRange(looprange);
  nloopx = eloopx - sloopx + 1;
  nloopy = eloopy - sloopy + 1;
   
  for(int j=0; j<nloopx; j++){
    double XEdgeOfPixelZ = ((top->z()-bottom->z())/(top->x()-bottom->x()))*(sximin+(sloopx+j)*_pixelSize-(top->x()+bottom->x())/2.)+(top->z()+bottom->z())/2.;
    double XEdgeOfPixelX = sximin+(sloopx+j)*_pixelSize;
    double XEdgeOfPixelY = ((top->y()-bottom->y())/(top->x()-bottom->x()))*(sximin+(sloopx+j)*_pixelSize-(top->x()+bottom->x())/2.)+(top->y()+bottom->y())/2.;
    
    gear::Vector3D* XEdgeOfPixel = new gear::Vector3D(XEdgeOfPixelX,XEdgeOfPixelY,XEdgeOfPixelZ);
    EdgeOfPixel.push_back(XEdgeOfPixel);
  }
  for(int j=0; j<nloopy; j++){
    double YEdgeOfPixelZ = ((top->z()-bottom->z())/(top->y()-bottom->y()))*(szetamin+(sloopy+j)*_pixelSize-(top->y()+bottom->y())/2.)+(top->z()+bottom->z())/2.;
    double YEdgeOfPixelX = ((top->x()-bottom->x())/(top->y()-bottom->y()))*(szetamin+(sloopy+j)*_pixelSize-(top->y()+bottom->y())/2.)+(top->x()+bottom->x())/2.;
    double YEdgeOfPixelY = szetamin+(sloopy+j)*_pixelSize;

    gear::Vector3D* YEdgeOfPixel = new gear::Vector3D(YEdgeOfPixelX,YEdgeOfPixelY,YEdgeOfPixelZ);
    EdgeOfPixel.push_back(YEdgeOfPixel);
  }
  EdgeOfPixel.push_back(top);
  EdgeOfPixel.push_back(bottom);
  return EdgeOfPixel;
}

// =====================================================================
map< pair<double, double>*, double> FPCCDDigitizer::getLocalPixel(SimTrackerHitImpl* simthit, vector<gear::Vector3D*> edgeofpixel){
  map<pair< double, double>*, double> pixel_local;
  int f_layer = simthit->getCellID0() - 1;
  double L_through_pixel = 0.;
  double dEdx = simthit->getEDep()*1e+9; // Energy deposit of 50um thickness ladder. [eV]
  double path_length = simthit->getPathLength(); // Path length of 50um thickness ladder. [mm]

  double mpv = 0.;
  double dE = 0.;
  double sigma = 0.;
  unsigned int count = 0;
  vector<gear::Vector3D*>::iterator nxt=edgeofpixel.begin();
  while( nxt != edgeofpixel.end() ){
    if( count == (unsigned int)edgeofpixel.size() - 1 ) break;
    vector<gear::Vector3D*>::iterator fst=nxt;
    nxt++;

    double diffx = (*nxt)->x() - (*fst)->x();
    double diffy = (*nxt)->y() - (*fst)->y();
    double diffz = (*nxt)->z() - (*fst)->z();
    double mag = sqrt(diffx*diffx+diffy*diffy+diffz*diffz);

    L_through_pixel = mag;
    mpv = (dEdx/path_length)*L_through_pixel;
    sigma = (_sigmaConst)*L_through_pixel;

    dE = gRandom->Landau( mpv, sigma)*1e-9; // Energy deposit is smeared by Landau distribution. [GeV]

    pair< double, double>* center_pixel = FindPixel((*fst),(*nxt),f_layer); // estimate central point of pixels 
    pixel_local.insert(map< pair<double, double>*, double>::value_type( center_pixel, dE));
    count++;
  }
  return pixel_local;
}

// =====================================================================
pair< double, double>* FPCCDDigitizer::FindPixel(gear::Vector3D* fst, gear::Vector3D* nxt, int layer){
  double hlength = _geodata[layer].hlength;
  double sximin = _geodata[layer].sximin;
  
  double fst_x = fst->x();
  double fst_y = fst->y();
  double nxt_x = nxt->x();
  double nxt_y = nxt->y();
  int i_x = -1000;
  int i_fst_x = 1000;
  int i_nxt_x = -1000;
  int i_y = -1000;
  int i_fst_y = 1000;
  int i_nxt_y = -1000;
  
  double limit = -1e-5;
  for(int i=-1; i<=_eloopx-_sloopx+1; i++){
    if( ( fst_x - sximin-(_sloopx + i)*_pixelSize) >= limit
     && (-fst_x + sximin+(_sloopx+i+1)*_pixelSize) >= limit )
      i_fst_x = i;
    if( ( nxt_x - sximin-(_sloopx + i)*_pixelSize) >= limit
     && (-nxt_x + sximin+(_sloopx+i+1)*_pixelSize) >= limit )
      i_nxt_x = i;
    if( i_fst_x == i_nxt_x ){ i_x = i; break; }
  }
  if( i_x == -1000 ){
    cout << "i_x == -1000" << endl;
  }
  for(int i=-1; i<=_eloopy-_sloopy+1; i++){
    if( ( fst_y + hlength-(_sloopy + i)*_pixelSize) >= limit
     && (-fst_y - hlength+(_sloopy+i+1)*_pixelSize) >= limit )
      i_fst_y = i;
    if( ( nxt_y + hlength-(_sloopy + i)*_pixelSize) >= limit
     && (-nxt_y - hlength+(_sloopy+i+1)*_pixelSize) >= limit )
      i_nxt_y = i;
    if( i_fst_y == i_nxt_y ){ i_y = i; break; }
  }
  if( i_y == -1000 ){
    cout << "i_y == -1000" << endl;
  }
  double pxl_cx = sximin+(_sloopx+i_x + 1)*_pixelSize+_pixelSize/2. ;
  double pxl_cy =-hlength+(_sloopy+i_y + 1)*_pixelSize+_pixelSize/2.;

  pair< double, double>* pixel_local = new pair< double, double>(pxl_cx, pxl_cy);

  return pixel_local;
}

// =====================================================================
void FPCCDDigitizer::makeNewSimTHit(SimTrackerHitImpl* simthit, gear::Vector3D* newpos, gear::Vector3D* newmom, int layer, int ladder, double newPathLength){
  float newdEdx = simthit->getEDep()*newPathLength/simthit->getPathLength();
  double   eta0 = _geodata[layer].rmin+0.5*_pixelheight + ((layer)%2)*(_geodata[layer].sthick-_pixelheight);
  double pos[3] = {newpos->x()*_geodata[layer].sinphi[ladder]+eta0*_geodata[layer].cosphi[ladder],
                  -newpos->x()*_geodata[layer].cosphi[ladder]+eta0*_geodata[layer].sinphi[ladder],
                   newpos->y()};
  float  mom[3] = {(float)(newmom->x()*_geodata[layer].cosphi[ladder] - newmom->z()*_geodata[layer].sinphi[ladder]),
                   (float)(newmom->x()*_geodata[layer].sinphi[ladder] + newmom->z()*_geodata[layer].cosphi[ladder]),
                   (float)newmom->y()};
  
  simthit->setEDep( newdEdx );
  simthit->setPathLength( (float)newPathLength );
  simthit->setPosition( pos );
  simthit->setMomentum( mom );
  return ;
}

// =====================================================================
int FPCCDDigitizer::inSensitiveRegion( gear::Vector3D* pos, int layer){
  if(fabs(pos->z()) > 0.5*_pixelheight + 1e-5
     || fabs(pos->y()) > _geodata[layer].hlength + 1e-5
     || pos->x() > _geodata[layer].sximax + 1e-5
     || pos->x() < _geodata[layer].sximin - 1e-5){
 return 0;
  }
  else{
 return 1;
  } 
}
