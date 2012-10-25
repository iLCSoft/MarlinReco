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
#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
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
  template<class T1,class T2,class T3>
  class SortByPosX
  {
  public:
    bool operator()( const T1& object1, const T1& object2 )
    {
      return object1.first->x() < object2.first->x();
    }
  };
}


// Start of New edit


 // This is a initialization of Each ladder's PixelSize.
                                                      // After this, PixelSizeVec may be changed by r.P.P..

// End of New edit




FPCCDDigitizer aFPCCDDigitizer ;



// =====================================================================
FPCCDDigitizer::FPCCDDigitizer() : Processor("FPCCDDigitizer") {
  
  // modify processor description
  _description = "FPCCDDigitizer should create FPCCD's VTXPixelHits from SimTrackerHits" ;
  
  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "Debug",
                              "Debugging option",
                              _debug,
                              int(0)); 

  registerProcessorParameter( "modifySimTrackerHit",
                              "Modify SimTrackerHit into layer",
                              _modifySimTHit,
                              bool(true)); 

  FloatVec PixelSizeVec;
  for(int i=0;i<6;i++){PixelSizeVec.push_back(0.005);}
  
  registerProcessorParameter( "Each_FPCCD_pixelSize(mm)",
                              "Each ladder's Pixel size of FPCCD (unit:mm) (default:0.005)",
                              _PixelSizeVec,
                              PixelSizeVec );

  
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
                              float(1.0));

  registerProcessorParameter( "Ladder_Number_encoded_in_cellID" ,
                              "Mokka has encoded the ladder number in the cellID" ,
                              _ladder_Number_encoded_in_cellID ,
                              bool(false));

  
  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "VTXCollectionName" , 
                           "Name of the VTX SimTrackerHit collection"  ,
                           _colNameVTX ,
                           std::string("VXDCollection") ) ;
  
  
  // Output collections
  registerOutputCollection( LCIO::LCGENERICOBJECT,
                            "VTXPixelHitCollection" , 
                            "Name of the VXD PixelHit output collection"  ,
                            _outColNameVTX ,
                            std::string("VTXPixelHits") ) ;
    
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
  const gear::ZPlanarParameters &gearVXD = Global::GEAR->getVXDParameters();
  const gear::ZPlanarLayerLayout &layerVXD = gearVXD.getVXDLayerLayout();
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
  
  col = 0 ;
  
  try{
    col = evt->getCollection( _colNameVTX ) ;
  }
  catch(DataNotAvailableException &e){
    if (_debug == 1) {   
      std::cout << "Collection " << _colNameVTX.c_str() << " is unavailable in event " << _nEvt << std::endl;
    }
  }
  
  if( col != 0 ){    

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
  } // End of process when VXD has hits
  _nEvt ++ ;
}

// =====================================================================
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
  const int cellId = SimTHit->getCellID0();

  if(_ladder_Number_encoded_in_cellID) {
    streamlog_out( DEBUG3 ) << "Get Layer Number using StandardILD Encoding from ILDConf.h : cellId = " << cellId << std::endl;
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;
    encoder.setValue(cellId) ;
    layer    = encoder[lcio::ILDCellID0::layer] ;
    ladderID = encoder[lcio::ILDCellID0::module] ;
    streamlog_out( DEBUG3 ) << "layer Number = " << layer << std::endl;
    streamlog_out( DEBUG3 ) << "ladder Number = " << ladderID << std::endl;
  }
  else{
    layer = cellId  - 1 ;
  }

  //----
  // check which ladder hit on
  //----

  if( !_ladder_Number_encoded_in_cellID ) {    
    ladderID = getLadderID( HitPosInMokka, layer);
  }
  
  //----
  // change pixel size due to the dependency of which ladder is used
  //----

  _pixelSize = _PixelSizeVec[layer];

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

 //----
  // get the intersections with the surface of sensitive region
  //----
  if( sqrt(MomAtLocalHitPos->x()*MomAtLocalHitPos->x() + MomAtLocalHitPos->z()*MomAtLocalHitPos->z()) < _momCut*1e-03){ // transvers momentum criteria for helix approximation
    if( _debug == 1 ) cout << "========== Track of this hit is trated as a helix ==========" << endl;
    getInOutPosOfHelixOnLadder(layer, PosOutFromLadder, PosInToLadder, LocalHitPos, MomAtLocalHitPos, LocalBField,charge);
  }
  else{
    getInOutPosOnLadder(layer, PosOutFromLadder,PosInToLadder,LocalHitPos,MomAtLocalHitPos); 
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
    
    std::vector<std::pair<const gear::Vector3D*,int> >EdgeOfPixel = getIntersectionOfTrkAndPix(layer,PosOutFromLadder,PosInToLadder);
    
    //----
    // calculate length that particle pass through pixel & central point of pixel
    //----
    
    sort( EdgeOfPixel.begin(), EdgeOfPixel.end(), SortByPosX<std::pair<const gear::Vector3D*, int>, const gear::Vector3D*, int >() );

    std::map< std::pair< int, int>*, double> pixel_local = getLocalPixel(SimTHit,EdgeOfPixel);

    std::map< std::pair< int, int>*, double>::iterator i_pixel_local = pixel_local.begin();
    while( i_pixel_local != pixel_local.end() ){
      int    xiID = (*i_pixel_local).first->first;
      int  zetaID = (*i_pixel_local).first->second;
      float    dE = (float)(*i_pixel_local).second; //  [mm]   
      FPCCDPixelHit::HitQuality_t quality;
      if( _isSignal && mcp->getParents().size() == 0 ){ quality = FPCCDPixelHit::kSingle; }
      else{ quality = FPCCDPixelHit::kBKG;}
      FPCCDPixelHit aHit(layer, ladderID, xiID, zetaID, dE, quality, mcp);
        hitVec.addPixelHit(aHit,quality);
       if(_debug ==1 ) aHit.print();

        i_pixel_local++;
    }
  }
  else{
    if(_debug == 1) cout << "This particle didn't through the sensitive region." << endl;
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
void FPCCDDigitizer::getInOutPosOnLadder(int layer, gear::Vector3D* outpos, gear::Vector3D* inpos, gear::Vector3D* pos,gear::Vector3D* mom){
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

  *pos = gear::Vector3D((outpos->x()+inpos->x())/2 ,
                        (outpos->y()+inpos->y())/2 ,
                        (outpos->z()+inpos->z())/2);
  
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
void FPCCDDigitizer::getInOutPosOfHelixOnLadder(int layer, gear::Vector3D* outpos, gear::Vector3D* inpos, gear::Vector3D* pos,gear::Vector3D* mom, gear::Vector3D* BField,float Charge){
  double out[3];
  double in[3];
  
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
  
  double      radius = (tmom.r()/(0.29979*BField->r()*fabs(Charge)))*1e+03;
  if( radius/2 < _pixelheight ){
    *outpos = gear::Vector3D( pos->x(), pos->y(), outerZ);
    *inpos  = gear::Vector3D( pos->x(), pos->y(), innerZ);
    *pos    = gear::Vector3D( pos->x(), pos->y(), 0);
    return ;
  }
  double init_pos[3] = {pos->x(),pos->y(),pos->z()};
  double init_dir[3] = {mom->x()/tmom.r(),
                        mom->y()/tmom.r(),
                        mom->z()/tmom.r()};
  double         v_y = radius*init_dir[1];
  
  double offset_pos[3];
  double offset_phi;
  double asin_offset = -(1/Charge)*asin(init_dir[0]/(Charge));    
  double acos_offset =  (1/Charge)*acos(init_dir[2]/(Charge));
  
  asin_offset*Charge >= 0 ? offset_phi = acos_offset : offset_phi = -acos_offset;
  int sign;
  fabs(acos_offset*Charge) > PI/2 ? sign = -1 : sign = 1;
  
  offset_pos[0] = init_pos[0]-radius*cos(Charge*offset_phi);
  offset_pos[1] = init_pos[1]-v_y*offset_phi;
  offset_pos[2] = init_pos[2]-radius*sin(Charge*offset_phi);
  
  double outerZ_phi = 0;
  double innerZ_phi = 0;
  //Check if the track intesects with x-y plane. 
  if( fabs((outerZ-offset_pos[2])/radius) <= 1. ) outerZ_phi = 1 / Charge*asin((outerZ-offset_pos[2])/radius) - asin_offset;
  if( fabs((innerZ-offset_pos[2])/radius) <= 1. ) innerZ_phi = 1 / Charge*asin((innerZ-offset_pos[2])/radius) - asin_offset;

  double outphi = sign*outerZ_phi ;
  double  inphi = sign*innerZ_phi ; 

  int flag = 0;
  if( radius*cos(Charge*(outphi + offset_phi)) + offset_pos[0] <= sximin
   || radius*cos(Charge*(outphi + offset_phi)) + offset_pos[0] >= sximax
   || v_y*(outphi + offset_phi) + offset_pos[1] <= szetamin
   || v_y*(outphi + offset_phi) + offset_pos[1] >= szetamax ) flag = flag + 1; // outerZ_phi will be modified.
  
  if( radius*cos(Charge*(inphi + offset_phi)) + offset_pos[0] <= sximin
   || radius*cos(Charge*(inphi + offset_phi)) + offset_pos[0] >= sximax
   || v_y*(inphi + offset_phi) + offset_pos[1] <= szetamin
   || v_y*(inphi + offset_phi) + offset_pos[1] >= szetamax ) flag = flag + 2; // innerZ_phi will be modified.
  
  if(flag != 0){
    
    double tmpphi;
    double   sximin_phi = 1e+8;
    double   sximax_phi = 1e+8;
    if(fabs((sximin-offset_pos[0])/radius) <= 1.) sximin_phi = sign*(1 / Charge*acos((sximin - offset_pos[0])/radius) - acos_offset);
    if(fabs((sximax-offset_pos[0])/radius) <= 1.) sximax_phi = sign*(1 / Charge*acos((sximax - offset_pos[0])/radius) - acos_offset);  
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
    if(flag == 2)  inphi = tmpphi;
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
  *pos = gear::Vector3D(radius*cos(Charge*(center_phi+offset_phi)) + offset_pos[0],
                        v_y*(center_phi+offset_phi) + offset_pos[1],
                        radius*sin(Charge*(center_phi+offset_phi)) + offset_pos[2]);
  *mom = gear::Vector3D(-radius*sin(Charge*(center_phi+offset_phi)),
                        v_y ,
                        radius*cos(Charge*(center_phi+offset_phi))); 
  *outpos = gear::Vector3D(out[0], out[1], out[2]);
  *inpos  = gear::Vector3D( in[0],  in[1],  in[2]);
  
  return ;
}

// =====================================================================
std::vector<std::pair<const gear::Vector3D*, int> > FPCCDDigitizer::getIntersectionOfTrkAndPix(const int f_layer, gear::Vector3D* top,gear::Vector3D* bottom){
  std::vector<std::pair<const gear::Vector3D*, int> > EdgeOfPixel;
  double   sximin =  _geodata[f_layer].sximin;
  double szetamin = -_geodata[f_layer].hlength;
  double small_localx;
  double large_localx;
  double small_localy;
  double large_localy;
  if(top->x()>=bottom->x()){ large_localx = top->x(); small_localx = bottom->x(); }
  else { large_localx = bottom->x(); small_localx = top->x(); }
  if(top->y()>=bottom->y()) {large_localy = top->y(); small_localy = bottom->y(); }
  else { large_localy = bottom->y(); small_localy = top->y(); }

  int sloopx = (int)ceil((small_localx-sximin)/_pixelSize) ;
  int eloopx = (int)ceil((large_localx-sximin)/_pixelSize) - 1;
  int sloopy = (int)ceil((small_localy-szetamin)/_pixelSize) ;
  int eloopy = (int)ceil((large_localy-szetamin)/_pixelSize) - 1;

  int nloopx = eloopx - sloopx + 1;
  int nloopy = eloopy - sloopy + 1;
   
  for(int j=0; j<nloopx; j++){
    double XEdgeOfPixelZ = ((top->z()-bottom->z())/(top->x()-bottom->x()))*(sximin+(sloopx+j)*_pixelSize-(top->x()+bottom->x())/2.)+(top->z()+bottom->z())/2.;
    double XEdgeOfPixelX = sximin+(sloopx+j)*_pixelSize;
    double XEdgeOfPixelY = ((top->y()-bottom->y())/(top->x()-bottom->x()))*(sximin+(sloopx+j)*_pixelSize-(top->x()+bottom->x())/2.)+(top->y()+bottom->y())/2.;
    
    gear::Vector3D* XEdgeOfPixel = new gear::Vector3D(XEdgeOfPixelX,XEdgeOfPixelY,XEdgeOfPixelZ);
    EdgeOfPixel.push_back(std::pair<gear::Vector3D*,int>(XEdgeOfPixel,1));
  }
  for(int j=0; j<nloopy; j++){
    double YEdgeOfPixelZ = ((top->z()-bottom->z())/(top->y()-bottom->y()))*(szetamin+(sloopy+j)*_pixelSize-(top->y()+bottom->y())/2.)+(top->z()+bottom->z())/2.;
    double YEdgeOfPixelX = ((top->x()-bottom->x())/(top->y()-bottom->y()))*(szetamin+(sloopy+j)*_pixelSize-(top->y()+bottom->y())/2.)+(top->x()+bottom->x())/2.;
    double YEdgeOfPixelY = szetamin+(sloopy+j)*_pixelSize;

    gear::Vector3D* YEdgeOfPixel = new gear::Vector3D(YEdgeOfPixelX,YEdgeOfPixelY,YEdgeOfPixelZ);
    EdgeOfPixel.push_back(std::pair<gear::Vector3D*,int>(YEdgeOfPixel,2));
  }
  EdgeOfPixel.push_back(std::pair<gear::Vector3D*,int>(top,0));
  EdgeOfPixel.push_back(std::pair<gear::Vector3D*,int>(bottom,0));
  return EdgeOfPixel;
}

// =====================================================================
std::map< pair<int, int>*, double> FPCCDDigitizer::getLocalPixel(SimTrackerHitImpl* simthit, vector<std::pair<const gear::Vector3D*,int> > edgeofpixel){
  std::map<pair< int, int>*, double> pixel_local;
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;
  encoder.setValue(simthit->getCellID0()) ;
  int f_layer    = encoder[lcio::ILDCellID0::layer] ;
  double dEdx = simthit->getEDep()*1e+9; // Energy deposit of 50um thickness ladder. [eV]
  double path_length = simthit->getPathLength(); // Path length of 50um thickness ladder. [mm]

  double   mpv = 0.;
  double    dE = 0.;
  double sigma = 0.;
  std::vector<std::pair<const gear::Vector3D*,int> >::iterator nxt = edgeofpixel.begin();
  while( nxt != edgeofpixel.end()-1 ){
    std::vector<std::pair<const gear::Vector3D*,int> >::iterator fst = nxt;
    nxt++;

    double diffx = (*nxt).first->x() - (*fst).first->x();
    double diffy = (*nxt).first->y() - (*fst).first->y();
    double diffz = (*nxt).first->z() - (*fst).first->z();
    double L_through_pixel = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

    mpv = (dEdx/path_length)*L_through_pixel;
    sigma = (_sigmaConst)*L_through_pixel;

    dE = gRandom->Landau( mpv, sigma)*1e-9; // Energy deposit is smeared by Landau distribution. [GeV]

    std::pair< int, int>* pixID = FindPixel((*fst),(*nxt),f_layer); // estimate central point of pixels 
    pixel_local.insert(std::map< std::pair<int, int>*, double>::value_type( pixID, dE));
  }
  return pixel_local;
}

// =====================================================================
std::pair< int, int>* FPCCDDigitizer::FindPixel(std::pair<const gear::Vector3D*,int> fst, std::pair<const gear::Vector3D*,int> nxt, int layer){

  std::pair<int,int> fst_array[4] ;
  std::pair<int,int> nxt_array[4] ;

  makeCandidates(fst, fst_array, layer);    
  makeCandidates(nxt, nxt_array, layer);
  
  std::pair<int,int>* pixID = new std::pair< int, int>(0,0);
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      if(fst_array[i].first >= 0 && fst_array[i] == nxt_array[j] ){
        pixID->first = fst_array[i].first;
        pixID->second = fst_array[i].second;
        goto GotPixID;
      }
    }
  }
 GotPixID:

  return pixID;
}
// =====================================================================
void FPCCDDigitizer::makeCandidates( std::pair<const gear::Vector3D*,int> edge, std::pair<int,int>* cand_array, int layer){

  std::pair<int,int> candID1(-1,-1);
  std::pair<int,int> candID2(-1,-1);
  std::pair<int,int> candID3(-1,-1);
  std::pair<int,int> candID4(-1,-1);

  double hlength = _geodata[layer].hlength;
  double sximin  = _geodata[layer].sximin;

  double quotient_x = (edge.first->x()-sximin)/_pixelSize ;
  double quotient_y = (edge.first->y()+hlength)/_pixelSize ;

  double mergin = 1e-5;

  switch(edge.second)
    {
    case 0 : // intersection with ladder surface.
      candID1.first = (int)floor(quotient_x);
      if     ( quotient_x - candID1.first < mergin/_pixelSize )     candID2.first = candID1.first - 1;
      else if( quotient_x - candID1.first > 1 - mergin/_pixelSize ) candID2.first = candID1.first + 1;
      else candID2.first = candID1.first;
      
      candID1.second = (int)floor(quotient_y);
      if     ( quotient_y - candID1.second < mergin/_pixelSize )    candID2.second = candID1.second - 1;
      else if( quotient_y - candID1.second > 1 - mergin/_pixelSize )candID2.second = candID1.second + 1;      
      else candID2.second = candID1.second;
      break;

    case 1 : // intersection with pixel border in xi
      candID1.first  = (int)floor(quotient_x + 0.5);
      candID2.first  = candID1.first - 1;
      candID1.second = (int)floor(quotient_y);
      candID2.second = candID1.second;
      
      if     ( quotient_y - candID1.second < mergin/_pixelSize){
        candID3.first  = candID1.first;
        candID3.second = candID1.second - 1;
        candID4.first  = candID2.first;
        candID4.second = candID2.second - 1;
      }
      else if( quotient_y - candID1.second > 1 - mergin/_pixelSize){
        candID3.first  = candID1.first;
        candID3.second = candID1.second + 1;
        candID4.first  = candID2.first;
        candID4.second = candID2.second + 1;
      }
      break;

    case 2 : // intersection with pixel border in zeta
      candID1.first  = (int)floor(quotient_x);
      candID2.first  = candID1.first;
      candID1.second = (int)floor(quotient_y + 0.5);
      candID2.second = candID1.second - 1;

      if     ( quotient_x - candID1.first < mergin/_pixelSize){
        candID3.first  = candID1.first - 1;
        candID3.second = candID1.second;
        candID4.first  = candID2.first - 1;
        candID4.second = candID2.second;
      }
      else if( quotient_x - candID1.first > 1 - mergin/_pixelSize){
        candID3.first  = candID1.first + 1;
        candID3.second = candID1.second;
        candID4.first  = candID2.first + 1;
        candID4.second = candID2.second;
      }break;

    default :
      if(_debug == 1) std::cout << "pixel candidates cannot found. : " << edge.second << std::endl;
      break;
    }
  *cand_array = candID1;
  *(cand_array+1) = candID2;
  *(cand_array+2) = candID3;
  *(cand_array+3) = candID4;
    
  return ;
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

  CellIDEncoder<SimTrackerHitImpl> cellid_encoder( lcio::ILDCellID0::encoder_string , col) ;
  cellid_encoder[ lcio::ILDCellID0::subdet ] = lcio::ILDDetID::VXD ;
  cellid_encoder[ lcio::ILDCellID0::side   ] = 0 ;
  cellid_encoder[ lcio::ILDCellID0::layer  ] = layer ;
  cellid_encoder[ lcio::ILDCellID0::module ] = ladder ;
  cellid_encoder[ lcio::ILDCellID0::sensor ] = 0 ;
  cellid_encoder.setCellID( simthit );
  simthit->setEDep( newdEdx );
  simthit->setPathLength( (float)newPathLength );
  simthit->setPosition( pos );
  simthit->setMomentum( mom );
  return ;
}

// =====================================================================
bool FPCCDDigitizer::inSensitiveRegion( gear::Vector3D* pos, int layer){
  if(fabs(pos->z()) > 0.5*_pixelheight + 1e-5
     || fabs(pos->y()) > _geodata[layer].hlength + 1e-5
     || pos->x() > _geodata[layer].sximax + 1e-5
     || pos->x() < _geodata[layer].sximin - 1e-5){
 return false;
  }
  else return true;
}
