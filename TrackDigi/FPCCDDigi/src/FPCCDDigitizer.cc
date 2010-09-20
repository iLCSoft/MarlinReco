/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "FPCCDDigitizer.h"
#include "FPCCDPixelHit.h"
#include "FPCCDData.h"

#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gearimpl/Vector3D.h>
#include "UTIL/LCTOOLS.h"

#include <cmath>
#include <algorithm>
#include <sstream>
#include <set>
#include <map>
#include <vector>

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
  
  registerProcessorParameter( "Debug",
                              "Debugging option",
                              _debug,
                              int(0)); 

  registerProcessorParameter( "FPCCD_PixelSize" ,
                              "Pixel size of FPCCD (unit:mm) (default: 0.005)"  ,
                              _pixelSize ,
                              float(0.005) ) ;

  registerProcessorParameter( "PixelSizeX" , 
                              "Pixel Size of local X direction(mm)",
                              _pixelsizex ,
                              float(0.005));

  registerProcessorParameter( "PixelSizeY" , 
                              "Pixel Size of local Y direction(mm)",
                              _pixelsizey ,
                              float(0.005));
  
  registerProcessorParameter( "PixelHeight" , 
                              "Pixel Height(mm)",
                              _pixelheight ,
                              float(0.015));

  registerProcessorParameter("ElectronsPerKeV",
                             "Electrons per keV",
                             _electronsPerKeV,
                             double(3620));
  
  registerProcessorParameter("Threshold",
                             "Cell Threshold in electrons",
                             _threshold,
                             30.);
  
  registerProcessorParameter("HitQuality",
                             "Hit Quality of the Event",
                             _isSignal,
                             bool(true));
  
  // Input collections
  registerInputCollection( LCIO::SIMTRACKERHIT,
                           "VTXCollectionName" , 
                           "Name of the VTX SimTrackerHit collection"  ,
                           _colNameVTX ,
                           string("VXDCollection") ) ;
  
  // Output collections
  registerOutputCollection( LCIO::LCGENERICOBJECT,
                           "VTXPixelHitCollection" , 
                           "Name of the vxd PixelHit output collection"  ,
                           _outColNameVTX ,
                           string("VTXPixelHits") ) ;

  
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

  const gear::VXDParameters &gearVXD = Global::GEAR->getVXDParameters();
  const gear::VXDLayerLayout &layerVXD = gearVXD.getVXDLayerLayout();
  
  _nLayer = layerVXD.getNLayers() ;
  _geodata.resize(_nLayer);
  _maxLadder = 0;
  
  for(int ly = 0;ly<_nLayer;ly++){
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
    _geodata[ly].cosphi.resize(_geodata[ly].nladder);
    _geodata[ly].sinphi.resize(_geodata[ly].nladder);
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
void FPCCDDigitizer::processEvent( LCEvent* evt )
{

  LCCollection* col=0;
  try {
    col=evt->getCollection( _colNameVTX ) ;
  }
  catch(DataNotAvailableException &e){
    if (_debug == 1) {
      cout << "Collection " << _colNameVTX.c_str() << " is unavailable in event " 
                << _nEvt << endl;
    }
  }

  
  int iColl=0;
  if( col != 0 ){    
    LCCollectionVec* fpccdDataVec = new LCCollectionVec( LCIO::LCGENERICOBJECT )  ;

    FPCCDData theData(_nLayer, _maxLadder);  // prepare object to make pixelhits

    int nSimHits = col->getNumberOfElements()  ;
    
    // Loop over all VXD hits

    
    for(int i=0; i< nSimHits; i++){
      SimTrackerHit* SimTHit = dynamic_cast<SimTrackerHit*>( col->getElementAt( i ) ) ;
      
      makePixelHits(SimTHit, theData);
      
    } // End of loop over all SimTrackerHits
    
    theData.packPixelHits( *fpccdDataVec );
    
    evt->addCollection( fpccdDataVec , _outColNameVTX ) ;
  } // End of process when VXD has hits
  
  std::cout << "dumpEvent Digitizer" << std::endl;
  LCTOOLS::dumpEvent( evt ) ;
  _nEvt ++ ;
  
}

// =====================================================================
// void FPCCDDigitizer::makePixelHits(const SimTrackerHit *simHit,
//                                    FPCCDData &hitVec)
// {
//   // Obtain ladderID, xiID and zetaID and fill them into FPCCDPixelHits,
//   // which is an array of FPCCDLadder_t

//   int layer = simHit->getCellID() - 1;
//   const double *pos =  simHit->getPosition() ;  
//   gear::Vector3D posvec(pos[0],pos[1],pos[2]);

//   // Convert SimTrackerHit position into pixelID;

//   FPCCDID_t fpccdID=encodeFPCCDID(layer, posvec);
//   if( fpccdID.layer >= 0 ) {  // When valid fpccdid is obtained, ..
    
//     float edep ;
//     float dedxSmear = 0.0 ;
//     edep = simHit->getEDep() ;

//     MCParticle *mcp ;
//     mcp = simHit->getMCParticle() ;
//     FPCCDPixelHit::HitQuality_t quality=FPCCDPixelHit::kSingle;

//     FPCCDPixelHit aHit(fpccdID.layer, fpccdID.ladder, fpccdID.xi, fpccdID.zeta,
//                        edep, quality, mcp);
    
//     hitVec.addPixelHit(aHit, true);

//   }

// }

void FPCCDDigitizer::makePixelHits(const SimTrackerHit *SimTHit,
                                   FPCCDData &hitVec)
{
  int iColl = 0;
  
  if( iColl==0 ){
    const int celId = SimTHit->getCellID() ;

    gear::Vector3D* HitPosInMokka = new gear::Vector3D(SimTHit->getPosition()[0],SimTHit->getPosition()[1],SimTHit->getPosition()[2]);
    cout << "HitPosInMokka: " << HitPosInMokka->x() << ", " << HitPosInMokka->y() << ", " << HitPosInMokka->z() << endl;
    
    //----
    // get basic info.
    //----
    const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters();
    const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout();
    const gear::BField& gearBField = Global::GEAR->getBField();
    
    double PI = M_PI;
    double TWOPI = 2*PI;
    double posphi = HitPosInMokka->phi();
    
    if(HitPosInMokka->y()<0) posphi += TWOPI;
    
    int cellID_c;
    cellID_c = SimTHit->getCellID();
    
    int layer;
    layer = SimTHit->getCellID() - 1;
    
    int nladder = _geodata[layer].nladder;     // Number of ladders in this layer
    double rmin = _geodata[layer].rmin;        // Distance of sensitive area from IP
    double dphi = _geodata[layer].dphi; 
    double sthick = _geodata[layer].sthick;
    double sximin = _geodata[layer].sximin;
    double sximax = _geodata[layer].sximax;
    double sxiwidth = sximax-sximin;
    double hlength = _geodata[layer].hlength;

    //----
    // get hit dir(mom) at hit points
    //----

    gear::Vector3D* MomAtHitPos = new gear::Vector3D(SimTHit->getMomentum()[0],SimTHit->getMomentum()[1],SimTHit->getMomentum()[2]);
    double momphi = MomAtHitPos->phi();
    if(MomAtHitPos->y()<0) momphi += TWOPI;
    
    gear::Vector3D origin;
    gear::Vector3D bfield = gearBField.at(origin);
    gear::Vector3D* BField = new gear::Vector3D(bfield.x(),bfield.y(),bfield.z());
    MCParticle* mcp = SimTHit -> getMCParticle();
    float charge = mcp->getCharge();
    //----
    //  print Basic Info. to Global variables
    //----
    //             printBasicInfoToGlobal(nhits,SimTHit,&layerVXD);
    //             cout << "printBasicInfoToGlobal" << endl;
    
    //----
    // check which ladder hit on
    //----
    
    int ladderID = getLadderID(SimTHit);
     
    //----
    // get local pos and dir on each ladder
    //----
    double tolocalrad = getToLocalRad(SimTHit, ladderID);
    gear::Vector3D* LocalHitPos = getLocalPos(SimTHit, tolocalrad);
    gear::Vector3D* MomAtLocalHitPos = getMomAtLocalPos(SimTHit, tolocalrad);
    gear::Vector3D* LocalBField = new gear::Vector3D(sin(tolocalrad)*BField->x()-cos(tolocalrad)*BField->y(),BField->z(),cos(tolocalrad)*BField->x()+sin(tolocalrad)*BField->y());
    
    gear::Vector3D* PosTopOnLadder;
    gear::Vector3D* PosBottomOnLadder;
    
    bool treatashelix = 0;
    
    if(MomAtLocalHitPos->z()/MomAtLocalHitPos->r() < sthick/sxiwidth || MomAtLocalHitPos->r() < 0.1){
      //----
      // Helix Fitting
      //----             
      
      cout << "========== Track of this hit is trated as a helix ==========" << endl;
      
      PosTopOnLadder = getInOutPosOfHelixOnLadder(SimTHit,LocalHitPos,MomAtLocalHitPos,LocalBField,charge,"in");
      PosBottomOnLadder = getInOutPosOfHelixOnLadder(SimTHit,LocalHitPos,MomAtLocalHitPos,LocalBField,charge,"out");
      if(mcp->isStopped()){
        gear::Vector3D* topOfHit = new gear::Vector3D(LocalHitPos->x(),LocalHitPos->y(),0);
        PosTopOnLadder = topOfHit;
        gear::Vector3D* bottomOfHit = new gear::Vector3D(LocalHitPos->x(),LocalHitPos->y(),0);
        PosBottomOnLadder = bottomOfHit;
        cout << "NOTICE:This particle stopped in ladder." << endl;
      }
      treatashelix = 1;
      
      if(fabs(PosTopOnLadder->z()) > 0.5*_pixelheight + 1e-05
         || fabs(PosBottomOnLadder->z()) > 0.5*_pixelheight + 1e-05
         || fabs(PosTopOnLadder->y()) > hlength + 1.E-02
         || fabs(PosBottomOnLadder->y()) > hlength + 1.E-02
         || PosTopOnLadder->x() > sximax
         || PosTopOnLadder->x() < sximin
         || PosBottomOnLadder->x() > sximax
         || PosBottomOnLadder->x() < sximin) treatashelix = 0;
    }
    if(treatashelix == 0 ){

      //----
      //Line Fittng
      //----
      // get hit 2 point on upper or lower plane of pixel
      //----
      
      PosTopOnLadder = getTopBottomPosOnLadder(LocalHitPos,MomAtLocalHitPos,"top");
      PosBottomOnLadder = getTopBottomPosOnLadder(LocalHitPos,MomAtLocalHitPos,"bottom");

      //----
      // modify the 2 point on upper or lower plane of pixel if the 2 points is not on the plane of ladder
      //----
      ModifyIntoLadder(PosTopOnLadder,SimTHit,LocalHitPos,MomAtLocalHitPos);
      ModifyIntoLadder(PosBottomOnLadder,SimTHit,LocalHitPos,MomAtLocalHitPos);
      
    }                   
    
    //----
    // get hit points incoming and outgoiong each pixel (@ edge of pixel)
    //----
    vector<gear::Vector3D*> EdgeOfPixel = getIntersectionOfTrkAndPix(SimTHit,LocalHitPos,MomAtLocalHitPos,PosTopOnLadder,PosBottomOnLadder);
    vector<gear::Vector3D*>::iterator i_EdgeOfPixel=EdgeOfPixel.begin();
    
    //----
    // calculate length that particle pass through pixel & central point of pixel
    //----
    sort( EdgeOfPixel.begin(), EdgeOfPixel.end(), SortByPosX<gear::Vector3D*>() );
    map<gear::Vector3D*, double> pixel_local = getLocalPixel(SimTHit,EdgeOfPixel);
    map<gear::Vector3D*, double>::iterator i_pixel_local=pixel_local.begin();
    while( i_pixel_local != pixel_local.end() ){
      double xi = ((*i_pixel_local).first->x()-sximin-_pixelsizex/2);
      double zeta = ((*i_pixel_local).first->y()+hlength-_pixelsizey/2);
      unsigned short int xiid = (unsigned short int)(xi/_pixelsizex);
      unsigned short int zetaid = (unsigned short int)(zeta/_pixelsizey);
      
      
//       cout << "Local_X : " << (*i_pixel_local).first->x() << endl;
//       cout << "Local_Y : " << (*i_pixel_local).first->y() << endl;
      
//       cout << "xi : " << xi << " xiid : " << xiid << endl;
//       cout << "zeta : " << zeta << " zetaid : " << zetaid << endl;
      
      float de_dx = (*i_pixel_local).second;
      FPCCDPixelHit::HitQuality_t quality;
      if(de_dx*_electronsPerKeV*1e+06 >= _threshold ){            // the calculation of energy deposit should be improved. 
        if( _isSignal ){ quality = FPCCDPixelHit::kSingle; }
        else{ quality = FPCCDPixelHit::kBKG; }
        FPCCDPixelHit aHit(layer,ladderID,xiid,zetaid,de_dx,quality,mcp);
        hitVec.addPixelHit(aHit,_isSignal);
        aHit.print();
      }
      i_pixel_local++;
    }
    cout << endl;
    cout << "============================================================================" << endl;
    cout << endl;
  }
}

// =====================================================================
FPCCDID_t FPCCDDigitizer::encodeFPCCDID(const int layer, const gear::Vector3D hitvec)
{
  /*  
   * Get hit ID (ladder, xiID, zetaID) from layer number and hit position (hitvec)
   * in a global coordinate system.  (xi, eta, zeta) is a coordinate in a local 
   * coordinate system defined in each ladder.
   *
   *  (X,Y,Z) are converted to (xi, eta, zeta) as follows.
   *  vector(center) is a position vector of the center of ladder.
   *  vector(center)=(r_min*cos(phi), r_min*sin(phi), 0);
   *    phi is an angle defined by ladderID and r_min is a distance between 
   * IP and ladder
   *  vector(eta_axis)= ( cos(phi), sin(phi), 0 ) : 
   *    same direction as vector(ladder_center)
   *  vector(xi_axis)=( +sin(phi), -cos(phi), 0)  : 
   *    perpendicular to vector(eta_axis) 
   *  vector(zeta_axis)=(0, 0, 1);
   *
   *  then 
   *    position in lab coordinate, pos(X,Y,Z) is given by
   *  pos(X,Y,Z) = vector(center) + xi*vector(xi_axis) 
   *             + eta*vector(eta_axis) + zeta+vector(zeta_axis)
   *  Solving this vector equation, we obtain
   *   xi=sin(phi)*(X - r_min*cos(phi)) - cos(phi)*( Y - r_min*sin(phi) )
   *   eta=cos(phi)*(X - r_min*cos(phi)) + sin(phi)*( Y - r_min*sin(phi) )
   *   zeta=Z
   *  
   *  LadderID is decided by checking that (xiID, etaID, zetaID) is 
   * in a sensitive region of a ladder.
   *
   *  xiID and zetaID are pixelID.  pixelID is 0 or positive. 
   * not negative. a pixel with minimum (xi, zeta) has 
   * pixelID (0,0).  xiID and zetaID is calculated as follows.
   *    xiID=(xi-xi_minimum)/pixel_size
   *    zetaID=(zeta-zeta_minimum)/pixel_size
   *  where
   *   xi_minimum = - half_width - xi_offset.
   *   zeta_minimum = - half_length
   *  Note that the center of ladder is shifted by offset.
   */

  int nladder=_geodata[layer].nladder;     // Number of ladders in this layer
  double rmin=_geodata[layer].rmin;        // Distance of sensitive area from IP
  double dphi=_geodata[layer].dphi; 
  double sthick=_geodata[layer].sthick;
  double sximin=_geodata[layer].sximin;
  double sximax=_geodata[layer].sximax;
  double hlength=_geodata[layer].hlength;
  
  int ladder=-999;
  int xiID=-999;
  double detamin=99999.0;
  double dximin=99999.0;
  double mergin=1.E-8;
  for(int iphi=0;iphi<nladder;iphi++) {
    double cosphi=_geodata[layer].cosphi[iphi];
    double sinphi=_geodata[layer].sinphi[iphi];
    double xi=sinphi*(hitvec.x()-rmin*cosphi)-cosphi*(hitvec.y()-rmin*sinphi);
    double eta=cosphi*(hitvec.x()-rmin*cosphi)+sinphi*(hitvec.y()-rmin*sinphi);
    if( ( eta >= -mergin && eta <= sthick+mergin  ) 
        && ( xi >= sximin && xi <= sximax ) ) {
      ladder=iphi;
      xiID=(int)((xi-sximin)/_pixelSize);
      break;
    }
    if( eta < 0 ) { detamin = -eta < detamin ? -eta : detamin ; }
    if( eta > sthick ) { detamin = eta -sthick < detamin ? eta-sthick : detamin ; }
    if( xi<sximin ) { dximin= sximin-xi < dximin  ?  sximin-xi : dximin ; }
    if( xi>sximax ) { dximin= xi-sximax < dximin  ?  xi-sximax : dximin ; }
  }
  FPCCDID_t fpccdID;
  if( ladder < 0 ) {
    fpccdID.layer=-1;
    cout << "=== Warning FPCCDDigitizer::encodeFPCCDID  given point is not within FPCCD sensitive region ... " << endl;
    cout << " layer =" << layer 
         << " hitpos=" << hitvec.x() << "," << hitvec.y() 
         << ", " << hitvec.z() << endl;
    cout << " rmin=" << rmin << " dphi=" << dphi << endl;
    cout << " sthick=" << sthick 
         << " hlength=" << hlength << endl; 
    cout << " sximin=" << sximin << " sximax=" << sximax << endl;
    cout << " detamin=" << detamin << " dximin=" << dximin << endl;

    return fpccdID;
  }
  fpccdID.layer=layer;
  fpccdID.ladder=ladder;
  fpccdID.xi=xiID;
  fpccdID.zeta=(int)((hitvec.z()+hlength)/_pixelSize);

  return fpccdID;

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

void FPCCDDigitizer::printBasicInfoToGlobal(const int nhits, const SimTrackerHit* simthit){
  gear::Vector3D* f_pos = new gear::Vector3D(simthit->getPosition()[0],simthit->getPosition()[1],simthit->getPosition()[2]);
  gear::Vector3D* f_mom = new gear::Vector3D(simthit->getMomentum()[0],simthit->getMomentum()[1],simthit->getMomentum()[2]);
  int f_layer = simthit->getCellID()-1;
  
  delete f_pos;
  delete f_mom;
  return;
}

// =====================================================================
int FPCCDDigitizer::getLadderID(const SimTrackerHit* simthit){
  gear::Vector3D* f_pos = new gear::Vector3D(simthit->getPosition()[0],simthit->getPosition()[1],simthit->getPosition()[2]);
  int f_layer = simthit->getCellID() - 1;
  double f_layerthickness = _geodata[f_layer].sthick;
  double f_radius = _geodata[f_layer].rmin+0.5*f_layerthickness;
  int f_nladder = _geodata[f_layer].nladder;
  double f_aphi = _geodata[f_layer].dphi;
  double f_offsetphi = _geodata[f_layer].phi0;
  double f_posphi = f_pos->phi();
  double f_posR = f_pos->rho();
  
  int ladderID = 50;
  double local_phi = 50.;
  for(int j=0; j<f_nladder; j++){
    local_phi = f_posphi - f_aphi*j - f_offsetphi;
    if((f_posR*cos(local_phi)-f_radius >= -f_layerthickness)&&(f_posR*cos(local_phi)-f_radius <= f_layerthickness)){
      ladderID = j;
      break;
    }
  }
  if(ladderID>100) cout << "LADDERID>100!: " << ladderID << endl;
  local_phi = f_posphi - f_aphi*(ladderID+1) - f_offsetphi;
  if((f_posR*cos(local_phi)-f_radius >= -f_layerthickness)&&(f_posR*cos(local_phi)-f_radius <= f_layerthickness)){
    ladderID++;
  }
  if(ladderID==50){
    for(int j=0; j<f_nladder; j++){
      local_phi = f_posphi - f_aphi*j - f_offsetphi;
      if((f_posR*cos(local_phi)-f_radius >= -1.*f_layerthickness)&&(f_posR*cos(local_phi)-f_radius <= 1.*f_layerthickness)){
        ladderID = j;
        break;
      }
    }
  }
  if(ladderID==50) cout << "LADDERID==50!" << endl;
  if(ladderID>100) cout << "LADDERID>100!..." << ladderID << endl;
  delete f_pos;
  return ladderID;
}

// =====================================================================
double FPCCDDigitizer::getToLocalRad(const SimTrackerHit* simthit,const int ladderID){
  int f_layer = simthit->getCellID() - 1;
  int f_nladder = _geodata[f_layer].nladder;
  double f_aphi = _geodata[f_layer].dphi;
  double f_offsetphi = _geodata[f_layer].phi0;
  double f_tolocalrad = ladderID*f_aphi + f_offsetphi;
  
  return f_tolocalrad;
}

// =====================================================================
gear::Vector3D* FPCCDDigitizer::getLocalPos(const SimTrackerHit* simthit,const double tolocalrad){
  gear::Vector3D* f_pos = new gear::Vector3D(simthit->getPosition()[0],simthit->getPosition()[1],simthit->getPosition()[2]);
  double f_posR = f_pos->rho();
  double f_posphi = f_pos->phi();
  int f_layer = simthit->getCellID() - 1;
  double f_radius = _geodata[f_layer].rmin+0.5*_geodata[f_layer].sthick;
  
  gear::Vector3D* Pos = new gear::Vector3D(-1.*f_posR*sin(f_posphi-tolocalrad), f_pos->z(), f_posR*cos(f_posphi-tolocalrad)-f_radius);
  
  delete f_pos;
  return Pos;
}

// =====================================================================
gear::Vector3D* FPCCDDigitizer::getMomAtLocalPos(const SimTrackerHit* simthit,const double tolocalrad){
  gear::Vector3D* f_mom = new gear::Vector3D(simthit->getMomentum()[0],simthit->getMomentum()[1],simthit->getMomentum()[2]);

  gear::Vector3D* Mom = new gear::Vector3D(sin(tolocalrad)*f_mom->x()-cos(tolocalrad)*f_mom->y(), f_mom->z(), cos(tolocalrad)*f_mom->x()+sin(tolocalrad)*f_mom->y());
  
  delete f_mom;
  return Mom;
}

// =====================================================================
gear::Vector3D* FPCCDDigitizer::getTopBottomPosOnLadder(gear::Vector3D* pos,gear::Vector3D* mom,const char* topbottom){
  double f_z = 0.5*_pixelheight;
  double f_sign = 0.;
  if(strcmp(topbottom,"top")==0) f_sign = 1.;
  else if(strcmp(topbottom,"bottom")==0) f_sign = -1.;
  f_z = f_sign*f_z;
  
  //  y-z plane  //
  double f_y = (mom->y()/mom->z())*(f_z-pos->z())+pos->y();
  //  x-z plane  //
  double f_x = (mom->x()/mom->z())*(f_z-pos->z())+pos->x();
  
  gear::Vector3D* PosTopBottomOnLadder = new gear::Vector3D(f_x,f_y,f_z);
  return PosTopBottomOnLadder;
}

// =====================================================================
void FPCCDDigitizer::ModifyIntoLadder(gear::Vector3D* bemodifiedpos,const SimTrackerHit* simthit,gear::Vector3D* pos,gear::Vector3D* mom){
  int f_layer = simthit->getCellID() - 1;
  double f_LadderRphi = _geodata[f_layer].sximax-_geodata[f_layer].sximin;
  double f_LadderZOver2 = _geodata[f_layer].hlength;
  double f_stepx = f_LadderRphi/_pixelsizex;
  double f_stepy = 2*f_LadderZOver2/_pixelsizey;
  double f_LadderSensitiveOffset = (_geodata[f_layer].sximax+_geodata[f_layer].sximin)/2.;
  double StartLadderRphi = _geodata[f_layer].sximin;
  double StartLadderZ = - f_LadderZOver2;
  double LadderStartx = StartLadderRphi;
  double LadderEndx = StartLadderRphi+(f_stepx-1)*_pixelsizex;
  double LadderStarty = StartLadderZ;
  double LadderEndy = StartLadderZ+(f_stepy-1)*_pixelsizey;

  double preposx = bemodifiedpos->x();
  double preposy = bemodifiedpos->y();
  double preposz = bemodifiedpos->z();

  if(preposx<LadderStartx){
    cout << "this hit was modified Sx" << endl; 
    preposx = LadderStartx;
    preposy = (mom->y()/mom->x())*(preposx-pos->x())+pos->y();
    preposz = (mom->z()/mom->x())*(preposx-pos->x())+pos->z();
  }else if(preposx>LadderEndx){
    cout << "this hit was modified Ex" << endl; 
    preposx = LadderEndx;
    preposy = (mom->y()/mom->x())*(preposx-pos->x())+pos->y();
    preposz = (mom->z()/mom->x())*(preposx-pos->x())+pos->z();
  }
  if(preposy<LadderStarty){
    cout << "this hit was modified Sy" << endl; 
    preposy = LadderStarty;
    preposx = (mom->x()/mom->y())*(preposy-pos->y())+pos->x();
    preposz = (mom->z()/mom->y())*(preposy-pos->y())+pos->z();
  }else if(preposy>LadderEndy){
    cout << "this hit was modified Ey" << endl; 
    preposy = LadderEndy;
    preposx = (mom->x()/mom->y())*(preposy-pos->y())+pos->x();
    preposz = (mom->z()/mom->y())*(preposy-pos->y())+pos->z();
  }
  gear::Vector3D prepos(preposx,preposy,preposz);
  *bemodifiedpos = prepos;
  return ;
}

// =====================================================================
gear::Vector3D* FPCCDDigitizer::getInOutPosOfHelixOnLadder(const SimTrackerHit* simthit,gear::Vector3D* pos,gear::Vector3D* mom,gear::Vector3D* BField,float Charge,const char* inout){
  double f_x;
  double f_y;
  double f_z;
  if(mom->r() > 0.01){
    int f_layer = simthit->getCellID() - 1;
    double f_LadderRphi = _geodata[f_layer].sximax-_geodata[f_layer].sximin;
    double f_LadderSensitiveOffset = (_geodata[f_layer].sximax+_geodata[f_layer].sximin)/2.;
    double f_LadderZOver2 = _geodata[f_layer].hlength;
    double f_LayerThickness = _geodata[f_layer].sthick;
    double OuterLadderZ = 0.5*_pixelheight;
    double InnerLadderZ = -0.5*_pixelheight;
    double StartLadderX = _geodata[f_layer].sximin;
    double EndLadderX = _geodata[f_layer].sximax;
    double StartLadderY = f_LadderZOver2;
    double EndLadderY = -f_LadderZOver2;
        
    int f_sign = 0;
    if(strcmp(inout,"out")==0) f_sign = 1;
    else if(strcmp(inout,"in")==0) f_sign = -1;
    
    Charge = f_sign*Charge;
       
    gear::Vector3D* tmom = new gear::Vector3D(mom->x(),0,mom->z());
    //     gear::Vector3D* tmom = new gear::Vector3D((1-BField->x()/BField->r())*mom->x(),
    //                                               (1-BField->y()/BField->r())*mom->y(),
    //                                               (1-BField->z()/BField->r())*mom->z());
    
    double PI = M_PI;

    double radius = (tmom->r()/PI*BField->r())*1e+01;
    double v_y = f_sign*radius*mom->y()/sqrt(mom->x()*mom->x()+mom->z()*mom->z());
    double init_pos[3] = {pos->x(),pos->y(),pos->z()};
    double init_dir[3] = {f_sign*mom->x()/sqrt(mom->x()*mom->x()+mom->z()*mom->z()), 0, f_sign*mom->z()/sqrt(mom->x()*mom->x()+mom->z()*mom->z())};
    double offset_pos[3];
    double offset_phi;
    double outphi;
    double tmpasin = -(1/Charge)*asin(init_dir[0]/(Charge));    
    double tmpacos = (1/Charge)*acos(init_dir[2]/(Charge));
    
    if( tmpasin >= 0){ tmpacos >= PI/2 ? offset_phi = tmpacos: offset_phi = tmpasin; }
    else{ tmpacos >= PI/2 ? offset_phi = -tmpacos : offset_phi = tmpasin; }
    
    offset_pos[0] = init_pos[0]-radius*cos(Charge*offset_phi);
    offset_pos[1] = init_pos[1]-v_y*offset_phi;
    offset_pos[2] = init_pos[2]-radius*sin(Charge*offset_phi);

    double center_phi = 0;
    double OuterLadderZ_phi = 0;
    double invOuterLadderZ_phi = 0;
    double InnerLadderZ_phi = 0;
    double invInnerLadderZ_phi = 0;
    
    if(fabs((OuterLadderZ-offset_pos[2])/radius) <= 1 ){
      OuterLadderZ_phi = Charge*asin((OuterLadderZ-offset_pos[2])/radius)-offset_phi;
      invOuterLadderZ_phi = 1/OuterLadderZ_phi;
    } 
    if( fabs((InnerLadderZ-offset_pos[2])/radius) <=1 ){
      InnerLadderZ_phi = Charge*asin((InnerLadderZ-offset_pos[2])/radius)-offset_phi;
      invInnerLadderZ_phi = 1/InnerLadderZ_phi;
    }
    
    if(fabs(OuterLadderZ_phi + InnerLadderZ_phi) > fabs(OuterLadderZ_phi - InnerLadderZ_phi)){
      center_phi = - (OuterLadderZ_phi + InnerLadderZ_phi)/2.;
      OuterLadderZ_phi += center_phi;
      InnerLadderZ_phi += center_phi;
      invOuterLadderZ_phi = 1/OuterLadderZ_phi;
      invInnerLadderZ_phi = 1/InnerLadderZ_phi;
    }
    
    if(invOuterLadderZ_phi <= invInnerLadderZ_phi ) outphi = InnerLadderZ_phi;
    else outphi = OuterLadderZ_phi;
    
    if(outphi <= 0
       || radius*cos(Charge*(outphi + offset_phi)) + offset_pos[0] >= StartLadderX
       || radius*cos(Charge*(outphi + offset_phi)) + offset_pos[0] <= EndLadderX
       || v_y*(outphi+offset_phi) + offset_pos[1] >= StartLadderY
       || v_y*outphi + offset_pos[1] <= EndLadderY ){
      
      double StartLadderX_phi = Charge*acos((StartLadderX - offset_pos[0])/radius)-offset_phi+center_phi;
      double EndLadderX_phi = Charge*acos((EndLadderX - offset_pos[0])/radius)-offset_phi+center_phi;  
      double StartLadderY_phi = (StartLadderY - offset_pos[1])/v_y+center_phi;
      double EndLadderY_phi = (EndLadderY - offset_pos[1])/v_y+center_phi;
      
      double invStartLadderX_phi = 1/StartLadderX_phi;
      double invEndLadderX_phi = 1/EndLadderX_phi;
      double invStartLadderY_phi = 1/StartLadderY_phi;
      double invEndLadderY_phi = 1/EndLadderY_phi;    
      
      if(invStartLadderX_phi >= invEndLadderX_phi && invStartLadderX_phi >= invStartLadderY_phi && invStartLadderX_phi >= invEndLadderY_phi && invStartLadderX_phi >= 1./outphi) outphi = StartLadderX_phi;
      else if(invEndLadderX_phi >= invStartLadderY_phi && invEndLadderX_phi >= invEndLadderY_phi && invEndLadderX_phi >= 1./outphi) outphi = EndLadderX_phi;
      else if(invStartLadderY_phi >= invEndLadderY_phi && invStartLadderY_phi >= 1./outphi) outphi = StartLadderY_phi;
      else if(invEndLadderY_phi >= 1./outphi) outphi = EndLadderY_phi; 
    }
    f_x = radius*cos(Charge*(outphi+offset_phi-center_phi)) + offset_pos[0];
    f_y = v_y * (outphi+offset_phi-center_phi)+ offset_pos[1];
    f_z = radius*sin(Charge*(outphi+offset_phi-center_phi)) + offset_pos[2];
    delete tmom;
  }
  else{
    f_x = pos->x();
    f_y = pos->y();
    f_z = pos->z();
  }
  
  gear::Vector3D* InOutPosOfHelixOnLadder = new gear::Vector3D(f_x,f_y,f_z);
  
  return InOutPosOfHelixOnLadder;
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
vector<gear::Vector3D*> FPCCDDigitizer::getIntersectionOfTrkAndPix(const SimTrackerHit* simthit,gear::Vector3D* pos,gear::Vector3D* mom,gear::Vector3D* top,gear::Vector3D* bottom){
  vector<gear::Vector3D*> EdgeOfPixel;
  int f_layer = simthit->getCellID() - 1;
  double f_LadderRphi = _geodata[f_layer].sximax-_geodata[f_layer].sximin;
  double f_LadderZOver2 = _geodata[f_layer].hlength;
  double f_stepx = f_LadderRphi/_pixelsizex;
  double f_stepy = 2*f_LadderZOver2/_pixelsizey;
  double f_LadderSensitiveOffset = (_geodata[f_layer].sximax+_geodata[f_layer].sximin)/2.;
  double StartLadderRphi = _geodata[f_layer].sximin;
  double StartLadderZ = - f_LadderZOver2;
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
  sloopx = (int)ceil((small_localx-StartLadderRphi)/_pixelsizex) ;
  eloopx = (int)ceil((large_localx-StartLadderRphi)/_pixelsizex) - 1;
  sloopy = (int)ceil((small_localy-StartLadderZ)/_pixelsizey) ;
  eloopy = (int)ceil((large_localy-StartLadderZ)/_pixelsizey) - 1;

  int looprange[4] = {sloopx,eloopx,sloopy,eloopy};
  setLoopRange(looprange);
  nloopx = eloopx - sloopx + 1;
  nloopy = eloopy - sloopy + 1;
  const int nhit_edgepixel = nloopx + nloopy;
   
  for(int j=0; j<nloopx; j++){
    double XEdgeOfPixelZ = ((top->z()-bottom->z())/(top->x()-bottom->x()))*(StartLadderRphi+(sloopx+j)*_pixelsizex-(top->x()+bottom->x())/2.)+(top->z()+bottom->z())/2.;
    double XEdgeOfPixelX = StartLadderRphi+(sloopx+j)*_pixelsizex;
    double XEdgeOfPixelY = ((top->y()-bottom->y())/(top->x()-bottom->x()))*(StartLadderRphi+(sloopx+j)*_pixelsizex-(top->x()+bottom->x())/2.)+(top->y()+bottom->y())/2.;
    
    gear::Vector3D* XEdgeOfPixel = new gear::Vector3D(XEdgeOfPixelX,XEdgeOfPixelY,XEdgeOfPixelZ);
    EdgeOfPixel.push_back(XEdgeOfPixel);
  }
  for(int j=0; j<nloopy; j++){
    double YEdgeOfPixelZ = ((top->z()-bottom->z())/(top->y()-bottom->y()))*(StartLadderZ+(sloopy+j)*_pixelsizey-(top->y()+bottom->y())/2.)+(top->z()+bottom->z())/2.;
    double YEdgeOfPixelX = ((top->x()-bottom->x())/(top->y()-bottom->y()))*(StartLadderZ+(sloopy+j)*_pixelsizey-(top->y()+bottom->y())/2.)+(top->x()+bottom->x())/2.;
    double YEdgeOfPixelY = StartLadderZ+(sloopy+j)*_pixelsizey;
    gear::Vector3D* YEdgeOfPixel = new gear::Vector3D(YEdgeOfPixelX,YEdgeOfPixelY,YEdgeOfPixelZ);
    EdgeOfPixel.push_back(YEdgeOfPixel);
  }
  EdgeOfPixel.push_back(top);
  EdgeOfPixel.push_back(bottom);
  return EdgeOfPixel;
}

// =====================================================================
map<gear::Vector3D*, double> FPCCDDigitizer::getLocalPixel(const SimTrackerHit* simthit, vector<gear::Vector3D*> edgeofpixel){
  map<gear::Vector3D*, double> pixel_local;
  int f_layer = simthit->getCellID() - 1;
  double L_through_pixel = 0.;
  double edep = simthit ->getEDep();
  double dE = 0.;
  int count = 0;
  vector<gear::Vector3D*>::iterator nxt=edgeofpixel.begin();
  while( nxt != edgeofpixel.end() ){
    if( count == (unsigned int)edgeofpixel.size()-1 ) break;
    vector<gear::Vector3D*>::iterator fst=nxt;
    nxt++;
    double diffx = (*nxt)->x() - (*fst)->x();
    double diffy = (*nxt)->y() - (*fst)->y();
    double diffz = (*nxt)->z() - (*fst)->z();
    double mag = sqrt(diffx*diffx+diffy*diffy+diffz*diffz);
    L_through_pixel = mag;
    dE = edep*L_through_pixel;
    gear::Vector3D* center_pixel = FindPixel((*fst),(*nxt),f_layer,"central");      // estimate central point of pixels 
    pixel_local.insert(map< gear::Vector3D*, double>::value_type( center_pixel, dE));
    count++;
  }
  
  return pixel_local;
}

// =====================================================================
gear::Vector3D* FPCCDDigitizer::FindPixel(gear::Vector3D* f_fst, gear::Vector3D* f_nxt, int f_layer, char* updown){
  double wefer_rphi  = _geodata[f_layer].sximax-_geodata[f_layer].sximin;
  double wefer_z2 = _geodata[f_layer].hlength;
  double sioffset_rphi = (_geodata[f_layer].sximax+_geodata[f_layer].sximin)/2.;
  double swefer_rphi = _geodata[f_layer].sximin;
  double swefer_z2 = -wefer_z2;
  double layerthickness = _geodata[f_layer].sthick;
    
  double fst_x = f_fst->x();
  double fst_y = f_fst->y();
  double nxt_x = f_nxt->x();
  double nxt_y = f_nxt->y();
  int i_x = -1000;
  int i_fst_x = 1000;
  int i_nxt_x = -1000;
  char cnxt_x[500] = "";
  char cfst_x[500] = "";
  sprintf(cnxt_x, "%g", nxt_x);
  sprintf(cfst_x, "%g", fst_x);
  int i_y = -1000;
  int i_fst_y = 1000;
  int i_nxt_y = -1000;
  char cnxt_y[500] = "";
  char cfst_y[500] = "";
  sprintf(cnxt_y, "%g", nxt_y);
  sprintf(cfst_y, "%g", fst_y);
  
  double limit = -0.000001;
  for(int i=-1; i<=_eloopx-_sloopx+1; i++){
    if( (fst_x - swefer_rphi+(_sloopx+i)*_pixelsizex) >= limit
        && (swefer_rphi+(_sloopx+i+1)*_pixelsizex - fst_x) >= limit )
      i_fst_x = i;
    if( (nxt_x - swefer_rphi+(_sloopx+i)*_pixelsizex) >= limit
        && (swefer_rphi+(_sloopx+i+1)*_pixelsizex - nxt_x) >= limit )
      i_nxt_x = i;
    if( i_fst_x == i_nxt_x ){ i_x = i; break; }
  }
  if( i_x == -1000 ){
    cout << "i_x == -1000" << endl;
  }
  
  for(int i=-1; i<=_eloopy-_sloopy+1; i++){
    if( (fst_y - swefer_z2+(_sloopy+i)*_pixelsizey) >= limit
        && (swefer_z2+(_sloopy+i+1)*_pixelsizey - fst_y) >= limit )
      i_fst_y = i;
    if( (nxt_y - swefer_z2+(_sloopy+i)*_pixelsizey ) >= limit
        && (swefer_z2+(_sloopy+i+1)*_pixelsizey - nxt_y) >= limit )
      i_nxt_y = i;
    if( i_fst_y == i_nxt_y ){ i_y = i; break; }
  }
  if( i_y == -1000 ){
    cout << "i_y == -1000" << endl;
  }
    
  double pxl_cx = swefer_rphi+(_sloopx+i_x)*_pixelsizex+_pixelsizex/2.;
  double pxl_cy = swefer_z2+(_sloopy+i_y)*_pixelsizey+_pixelsizey/2.;
  double pxl_cz = 100.;
  if(strcmp(updown,"central")==0)
    pxl_cz = 0.;
  
  gear::Vector3D* pixel_local = new gear::Vector3D(pxl_cx,pxl_cy,pxl_cz);
  return pixel_local;
}
