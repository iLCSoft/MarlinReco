/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "FPCCDClustering.h"
#include "FPCCDPixelHit.h"
#include "FPCCDData.h"

#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>
#include <UTIL/CellIDEncoder.h>
#include "UTIL/LCTOOLS.h"

#include <ILDCellIDEncoding.h>

#include <gsl/gsl_randist.h>

#include <cmath>
#include <algorithm>
#include <sstream>
#include <map>
#include <vector>
#include <utility>
#include <stack>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

FPCCDClustering aFPCCDClustering ;

// =====================================================================
FPCCDClustering::FPCCDClustering() : Processor("FPCCDClustering") {
  
  // modify processor description
  _description = "FPCCDClustering icteats TrackerHits from FPCCDPixelHits" ;
  
  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "Debug",
                              "Debugging option",
                              _debug,
                              int(0)); 

  registerProcessorParameter( "EnergyDigitizatoin",
                              "switch for digitization of energy deposit",
                              _energyDigitization,
                              int(1)); 

  registerProcessorParameter( "RandomNoise",
                              "Random noise option",
                              _randomNoise,
                              int(0)); 

  registerProcessorParameter( "FPCCD_PixelSize" ,
                              "Pixel size of FPCCD (unit:mm) (default: 0.005)",
                              _pixelSize,
                               float(0.0050) ) ;
  
  registerProcessorParameter( "PixelHeight" , 
                              "Pixel Height(mm)",
                              _pixelheight ,
                              float(0.015));
  
  registerProcessorParameter( "PointResolutionRPhi" ,
                              "R-Phi Resolution in VTX"  ,
                              _pointResoRPhi ,
                              float(0.001440)) ; // 5/sqrt(12)
  
  registerProcessorParameter( "PointResolutionZ" ,
                              "Z Resolution in VTX" ,
                              _pointResoZ ,
                              float(0.001440));   // 5/sqrt(12)
  
  registerProcessorParameter("noiseElectronSigma",
                             "gaussian sigma of number of electron of noise",
                             _electronNoiseRate,
                             double(50));

  registerProcessorParameter("NumberOfBitsForEnergyDeposit",
                             "Number of Bits used for energy depsit",
                             _nbitsForEdep,
                             int(7));

  registerProcessorParameter("NumberOfElectronForStep",
                             "Number of electron for 1 step",
                             _electronsPerStep,
                             int(25));
  
  registerProcessorParameter("Threshold",
                             "Cell Threshold in electrons",
                             _threshold,
                             double(200));
  
  registerProcessorParameter("ElectronsPerKeV",
                             "Electrons per keV",
                             _electronsPerKeV,
                             double(276)); // 3.62eV/1electron
  

  registerProcessorParameter( "NewTrackingSystem",
                              "Use new tracking system",
                              _new_tracking_system,
                              bool(false)); 
  
  // Input collections
  registerInputCollection( LCIO::LCGENERICOBJECT,
                           "VTXPixelHitCollection" , 
                           "Name of the VTX PixelHit collection"  ,
                           _colNameVTX ,
                           std::string("VTXPixelHits") ) ;

  // Output collections
  registerOutputCollection( LCIO::TRACKERHIT,
                           "VTXHitCollection" ,
                           "Name of the vxd TrackerHit output collection"  ,
                           _outColNameVTX ,
                           std::string("VTXTrackerHits") ) ;
}


// =====================================================================
void FPCCDClustering::init() { 

  // usually a good idea to
  printParameters() ;
  std::cout << "printparameter " << std::endl;
  _nRun = 0 ;
  _nEvt = 0 ;

  InitGeometry();

  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_default_seed = _ranSeed;
  
}

// =====================================================================
void FPCCDClustering::InitGeometry() 
{ 
// Save frequently used parameters.

  const gear::VXDParameters &gearVXD = Global::GEAR->getVXDParameters();
  const gear::VXDLayerLayout &layerVXD=gearVXD.getVXDLayerLayout();

  _nLayer = layerVXD.getNLayers() ;
  _geodata.resize(_nLayer);
  _maxLadder = 0;
  
  for(int ly=0;ly<_nLayer;ly++){
    _geodata[ly].nladder = layerVXD.getNLadders(ly);     // Number of ladders in this layer
    if( _maxLadder < _geodata[ly].nladder ) { _maxLadder = _geodata[ly].nladder; }
    _geodata[ly].rmin = layerVXD.getSensitiveDistance(ly); // Distance of sensitive area from IP
    _geodata[ly].dphi = (2*M_PI)/(double)_geodata[ly].nladder;
    _geodata[ly].phi0 = layerVXD.getPhi0(ly);  // phi offset.
    _geodata[ly].sthick = layerVXD.getSensitiveThickness(ly);
    _geodata[ly].sximin = -layerVXD.getSensitiveOffset(ly)
                  		  	-layerVXD.getSensitiveWidth(ly)/2.0;
    _geodata[ly].sximax = -layerVXD.getSensitiveOffset(ly)
                  			  +layerVXD.getSensitiveWidth(ly)/2.0;
    _geodata[ly].hlength = layerVXD.getSensitiveLength(ly);
    _geodata[ly].cosphi.resize( _geodata[ly].nladder ) ;
    _geodata[ly].sinphi.resize( _geodata[ly].nladder ) ;
    _geodata[ly].ladder_incline.resize( _geodata[ly].nladder ) ;
    _geodata[ly].num_xi_pixel = (int)(layerVXD.getSensitiveWidth(ly)/_pixelSize);
    _geodata[ly].num_zeta_pixel = (int)(2*layerVXD.getSensitiveLength(ly)/_pixelSize);
    for(int ld=0;ld<_geodata[ly].nladder;ld++) {
      double phi = _geodata[ly].phi0 + _geodata[ly].dphi*ld;
      _geodata[ly].cosphi[ld] = cos(phi);
      _geodata[ly].sinphi[ld] = sin(phi);

      double incline =  phi + (M_PI/2.0);
      while( incline >  1.0*M_PI || incline < -1.0*M_PI ){ incline > 1.0*M_PI ? incline -= 2.0*M_PI : incline += 2.0*M_PI ;}
      _geodata[ly].ladder_incline[ld] = incline; 
    }
  }
} 


// =====================================================================
void FPCCDClustering::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  std::cout << "processrunheader " << std::endl;
} 

// =====================================================================
void FPCCDClustering::modifyEvent( LCEvent * evt )
{
  
  LCCollectionVec* pHitCol=0;
  try {
    pHitCol = dynamic_cast<LCCollectionVec*>( evt->getCollection( _colNameVTX ) );
    pHitCol->setTransient(true);
  }
  catch(DataNotAvailableException &e){
   if (_debug >= 1) {
        std::cout << "Collection " << _colNameVTX.c_str() << " is unavailable in event " 
                  << _nEvt << std::endl;
   }
  }
  
  if( _debug >= 1 ) {
    std::cout << " Collection =" << _colNameVTX << " nevt=" << _nEvt << std::endl;
    if( pHitCol != 0 ) {
      std::cout << " number of elements is " << pHitCol->getNumberOfElements() << std::endl;
    }
  }
  if( pHitCol != 0 ){    

    
    FPCCDData  theData(_nLayer, _maxLadder);  // prepare object to make pixelhits
    int nhit=theData.unpackPixelHits(*pHitCol);
    if( _debug >= 2 ) { LCTOOLS::dumpEvent( evt ) ;}
    
    if( _debug >=2 ) { std::cout << " nhit=" << nhit << std::endl; }
    if( nhit > 0 ) {  // Output Trackhit, if there are pixel hits
      LCCollectionVec* trkHitVec = new LCCollectionVec( LCIO::TRACKERHIT );
    
      makeTrackerHitVec(theData, *trkHitVec);
      evt->addCollection( trkHitVec, _outColNameVTX ) ;
      std::cout << "dumpEvent Clustering" << std::endl;
      for(int i = pHitCol->getNumberOfElements(); i>0; i--){
        pHitCol->removeElementAt(i-1);
      }
    }
    evt->removeCollection( _colNameVTX );
    theData.clear();

  } // End of process when VXD has hits
  _nEvt ++ ;
}

// ====================================================================
void FPCCDClustering::makeTrackerHitVec(FPCCDData &pHitCol, LCCollectionVec &trkHitVec)
{
  // Convert pixel hits to TrackerHits
  for(int layer=0;layer<_nLayer;layer++){
    int     nladder = _geodata[layer].nladder;
    double   sximin = _geodata[layer].sximin;
    double   sximax = _geodata[layer].sximax;
    double sxiwidth = sximax-sximin;
    double  hlength = _geodata[layer].hlength;

// ----------------------------------------
// Start clustering of hits in each laddder
// ----------------------------------------
    for(int ladder=0;ladder<nladder;ladder++){

      FPCCDLadderHit_t ladderHit;
      PixelHitMap_t::iterator it=pHitCol.itBegin(layer, ladder);      
//(1) Copy all hits into local variables, ladderHit
      while( it != pHitCol.itEnd(layer, ladder) ) {

        FPCCDPixelHit *aHit=(*it).second;
        // Set the noise
        aHit->setEdep( aHit->getEdep() + gsl_ran_gaussian(_rng, _electronNoiseRate/_electronsPerKeV)*1e-6);
        if( _energyDigitization){
          EnergyDigitizer( aHit );
        }
          if( aHit->getEdep() > ( _threshold/_electronsPerKeV ) * 1e-6 ) {  // Do the threshold cut
          
          int xiID=aHit->getXiID();
          int zetaID=aHit->getZetaID();
          FPCCDHitLoc_t hitloc(xiID, zetaID);
          ladderHit.insert(map<FPCCDHitLoc_t, FPCCDPixelHit*>::value_type(hitloc, aHit));
        }
        it++;
      }
// Set the random noise hit
      if( _randomNoise && _energyDigitization){
        unsigned int noiseXi = 0;
        unsigned int noiseZeta = 0;
        double noiseEdep = 0;
        
        for(unsigned int i=0; i< gsl_ran_poisson( _rng, (sxiwidth/_pixelSize)*(2*hlength/_pixelSize)*3.1671*1e-5) ; i++){ // Random noise is generated by 4 sigma plobability.
          
          noiseXi = (unsigned int)gsl_ran_flat( _rng, 0, sxiwidth/_pixelSize );
          noiseZeta = (unsigned int)gsl_ran_flat( _rng, 0, 2*hlength/_pixelSize );          
          FPCCDHitLoc_t noiseHitLoc( noiseXi, noiseZeta );
          FPCCDLadderHit_t::iterator noiseIt=ladderHit.find( noiseHitLoc );
          
          if( noiseIt == ladderHit.end() ){
            noiseEdep = gsl_ran_gaussian_tail( _rng, _threshold, _electronNoiseRate )/(_electronsPerKeV * 1e+6 );
            FPCCDPixelHit* noiseHit = new FPCCDPixelHit(layer, ladder, noiseXi, noiseZeta, noiseEdep, FPCCDPixelHit::kBKG, 0);
            EnergyDigitizer( noiseHit);            
            ladderHit.insert( map<FPCCDHitLoc_t, FPCCDPixelHit*>::value_type( noiseHitLoc, noiseHit) );
          }          
        }
      }
// Do clustering
 
        if( ladderHit.size() > 0 ) {
          FPCCDClusterVec_t clusterVec;
          makeClustersInALadder(layer, ladderHit, clusterVec);
          if( _debug >= 1 ) {
            std::cout << "Layer:" << layer << " ladder:" << ladder << " # of clusters=" << clusterVec.size() << std::endl;
            std::cout << "  # pixels in each clusters are "; 
            for(unsigned int i=0;i<clusterVec.size();i++) {
              std::cout << clusterVec[i]->size() << " " ;
              if ( i > 20 ) {
                std::cout << " ... omitting the rest. " ;
                break;
              }
            }
            std::cout << endl;
          }
          
          // Calulate position from cluster and 
          makeTrackerHit(layer, ladder, clusterVec, trkHitVec);
 
          // Now clean up clusters in clusterVec;
          for(int i=clusterVec.size()-1;i>=0;i--){
            delete clusterVec[i];
          }
        }
        
    } // End of Ladder loop
  } // End of Layer loop
}


//=============================================================================================
void FPCCDClustering::makeTrackerHit(int layer, int ladder, FPCCDClusterVec_t &clusterVec, LCCollectionVec &trkHitVec)
{
  
  CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder( ILDCellIDEncoding::encoder_string , &trkHitVec ) ;
  
  for(unsigned int ic=0;ic<clusterVec.size();ic++) {
    FPCCDCluster_t *cluster=clusterVec[ic];
    double   xiene = 0.0;
    double zetaene = 0.0;
    double  enesum = 0.0;
    int      maxxi = 0;
    int      minxi = 0;
    int    maxzeta = 0;
    int    minzeta = 0;
    FPCCDPixelHit::HitQuality_t trackquality = FPCCDPixelHit::kSingle;

    FPCCDPixelHit   *maxXiHit=(*cluster)[0];
    FPCCDPixelHit   *minXiHit=(*cluster)[0];
    FPCCDPixelHit *maxZetaHit=(*cluster)[0];
    FPCCDPixelHit *minZetaHit=(*cluster)[0];
    
    unsigned int nPix=cluster->size();    
    for(unsigned int i=0;i<nPix;i++) {
      FPCCDPixelHit *aHit=(*cluster)[i];
      trackquality = aHit->getQuality();
      enesum+=aHit->getEdep();
      xiene+=((double)(aHit->getXiID()+0.5))*_pixelSize*aHit->getEdep();
      zetaene+=((double)(aHit->getZetaID()+0.5))*_pixelSize*aHit->getEdep();

      if(i==0){
        maxxi = aHit->getXiID();
        minxi = maxxi;
        maxzeta = aHit->getZetaID();
        minzeta = maxzeta;
      }
      else{
        if(  aHit->getXiID() > maxxi){ maxxi = aHit->getXiID(); maxXiHit=aHit;}
        if(  aHit->getXiID() < minxi){ minxi = aHit->getXiID(); minXiHit=aHit;}
        if(aHit->getZetaID() > maxzeta){ maxzeta = aHit->getZetaID();maxZetaHit=aHit;}
        if(aHit->getZetaID() < minzeta){ minzeta = aHit->getZetaID();minZetaHit=aHit;}
          }
      FPCCDPixelHit::HitQuality_t addedQuality = aHit->getQuality();
      if( trackquality != FPCCDPixelHit::kBKGOverlap){
        if(trackquality == FPCCDPixelHit::kBKG){ addedQuality=FPCCDPixelHit::kBKG ? trackquality=FPCCDPixelHit::kBKG : trackquality=FPCCDPixelHit::kBKGOverlap;}         
      }
      else if(addedQuality==FPCCDPixelHit::kSignalOverlap){trackquality=FPCCDPixelHit::kSignalOverlap;       
      }
      
    }
    
    unsigned int xiWidth = maxxi-minxi + 1;
    unsigned int zetaWidth = maxzeta-minzeta + 1; // cluster shapes info. could be used to reject background.
    unsigned short int tilt=0;
    if( xiWidth==1 || zetaWidth==1 ) tilt=0;
    else{
      if(maxXiHit->getZetaID()-minXiHit->getZetaID()>=0 && maxZetaHit->getXiID()-maxZetaHit->getXiID()>=0) tilt=1;
      else if(maxXiHit->getZetaID()-minXiHit->getZetaID()<=0 && maxZetaHit->getXiID()-maxZetaHit->getXiID()<=0) tilt=2;
      else tilt=0;
    }

    double xiave = xiene/enesum;
    double zetaave = zetaene/enesum;
    double newpos[3];
    double eta0 = _geodata[layer].rmin+0.5*_pixelheight +(layer%2)*(_geodata[layer].sthick-_pixelheight);

    newpos[0] = (xiave+_geodata[layer].sximin)*_geodata[layer].sinphi[ladder]
                + eta0*_geodata[layer].cosphi[ladder];
    newpos[1] =-(xiave+_geodata[layer].sximin)*_geodata[layer].cosphi[ladder]
                + eta0*_geodata[layer].sinphi[ladder];
    newpos[2] = zetaave-_geodata[layer].hlength;
    //store hit variables
 
  
    if ( _new_tracking_system ) {
      TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl ; 
      trkHit->setType( 100 + layer);

      cellid_encoder[ ILDCellIDEncoding::Fields::subdet ] = ILDCellIDEncoding::DetID::VXD ;
      cellid_encoder[ ILDCellIDEncoding::Fields::layer  ] = layer ;
      cellid_encoder[ ILDCellIDEncoding::Fields::module ] = ladder ;
      //int side = 1;
      //SJA:FIXME: for now don't use side   
      //   (*cellid_encoder)[ ILDCellIDEncoding::Fields::side   ] = side ;      
      cellid_encoder[ ILDCellIDEncoding::Fields::side   ] = 0 ;      
     
      cellid_encoder.setCellID( trkHit );
      
      unsigned int cellid1=0;
      cellid1 = ((trackquality << 30) & 0xc0000000) |
        ((tilt << 28) & 0x30000000) |
        ((nPix << 18) & 0x0ffc0000) |
        ((zetaWidth << 9) & 0x0003fe00) |
        ((xiWidth) & 0x000001ff) ;
      
      trkHit->setCellID1(cellid1);
      float u_direction[2] ;
      u_direction[0] = _geodata[layer].ladder_incline[ladder] ;
      u_direction[1] = 0.0 ;

      float v_direction[2] ;
      v_direction[0] = 0.0 ;
      v_direction[1] = 0.0 ;

      trkHit->setU( u_direction ) ;
      trkHit->setV( v_direction ) ;
      
      trkHit->setdU( _pointResoRPhi ) ;
      trkHit->setdV( _pointResoZ ) ;      
      
      trkHit->setPosition( newpos ) ;
      trkHit->setEDep( enesum );
      trkHitVec.addElement( trkHit );
    }
    else{
      TrackerHitImpl* trkHit = new TrackerHitImpl ;
    
      trkHit->setPosition( newpos ) ;
      trkHit->setEDep( enesum );
      trkHit->setType( 101 + layer);
      
      float covMat[TRKHITNCOVMATRIX]={0.,0.,  _pointResoRPhi*_pointResoRPhi,0.,0.,_pointResoZ*_pointResoZ}; // Resolution depends on theta.
      trkHit->setCovMatrix(covMat);
      trkHitVec.addElement( trkHit );
    }
  }  // Repeat pixel hits in this ladder.
}

// =====================================================================
void FPCCDClustering::makeClustersInALadder(int layer, FPCCDLadderHit_t &ladderHit, FPCCDClusterVec_t &clusterVec)
{
//(2) Start pixel hit clustering 
//    All pixels adjucent pixels are clustered, if it has a hit.  Threshold is not considered yet

  int   maxxi = _geodata[layer].num_xi_pixel;
  int maxzeta = _geodata[layer].num_zeta_pixel;

  FPCCDLadderHit_t::iterator ih = ladderHit.begin();
  while( ih != ladderHit.end() ) {
    FPCCDPixelHit *h=(*ih).second;
    (*ih).second=0;
    ih++;
    if( h==0 ) { continue; }
    //  std::cerr << "Found a seed pixel.  Edep=" << h->getEdep() << std::endl;
    //  Found a ssed of a cluster
    
    FPCCDCluster_t *cluster=new FPCCDCluster_t();
    cluster->push_back(h);
    int xi0=h->getXiID();
    int zeta0=h->getZetaID();
    typedef std::pair<int, int> Direction_t;
    std::stack<Direction_t>  direction;   
    for(int i=-1;i<=1;i++) { // Set initial search direction, neighbouring 8 pixels are searched
      for(int j=-1;j<=1;j++) {
	if( i!=0 || j!= 0 ) { direction.push(Direction_t(i,j)); }
      }
    } 
// Now look for all directions and find pixel with hit.
    while( ! direction.empty() ) {
      Direction_t newdir=direction.top();
      int xi=xi0+newdir.first;
      int zeta=zeta0+newdir.second;
//          std::cerr << "searching at (xi,zeta)=" << xi << "  " << zeta << std::endl;
      if( xi < 0 || xi > maxxi || zeta < 0 || zeta > maxzeta ) {
	direction.pop();
	continue;
      }
      FPCCDHitLoc_t newloc(xi, zeta);
      FPCCDLadderHit_t::iterator look=ladderHit.find(newloc);
      if( look == ladderHit.end() ) {
	direction.pop();
	continue;
      }
      FPCCDPixelHit *aHit=(*look).second;
      if( aHit == 0 ) {
  direction.pop();
  continue;
      }
      cluster->push_back( aHit );
//          std::cout << " Found hit .. at xi, zeta=" << xi << " " << zeta << std::endl;
      (*look).second=0;
      for (int ix=-1;ix<=1;ix++) {
        for(int iz=-1;iz<=1;iz++) {
	  if( ix != 0 || iz != 0 ) { 
	     direction.push(Direction_t(newdir.first+ix, newdir.second+iz));
 	  }
        }
      }
    } // Complete search of all nighbours of the seed
    if( cluster->size() > 0 ) {
      clusterVec.push_back(cluster);  
    }
  } 
}


// =====================================================================
void FPCCDClustering::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

  
// =====================================================================
void FPCCDClustering::end(){ 
  
  streamlog_out(MESSAGE) << " end()  " << name() 
                         << " processed " << _nEvt << " events in " << _nRun << " runs "
                         << std::endl ;

}


// =================================================================
// =================================================================
void FPCCDClustering::EnergyDigitizer( FPCCDPixelHit* aHit ){

  int nEle = (int)(aHit->getEdep()*1e+6 * _electronsPerKeV) ;
  int nStep = 1 << _nbitsForEdep ;
  int count = nEle/ _electronsPerStep ; //  int devided by int equales int. 
  
  if(count > nStep ) count = nStep ;

  aHit->setEdep( (count * _electronsPerStep)*1e-6 / _electronsPerKeV) ;
  return ;
}
