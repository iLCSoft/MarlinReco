#include "CurlKillerProcessor.h"
#include <iostream>

//#include "MarlinCED.h"

#include "constants.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gearimpl/FixedPadSizeDiskLayout.h>
//

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace constants ;

// function to generate a unique key for each [r][phi] bin based on a 64bit bitfield
// the most significant 32 bits are used for r 
unsigned long long make_key( unsigned r, unsigned phi ){

  unsigned long long temp = 0xffffffff & r ;

  temp = temp << 32 ;

  return  ( temp )  |  ( 0xffffffff & phi )  ;
} 


// Map to store the enteries for a 2D(r-phi) projection of the Tracker hits 
typedef std::map< unsigned long long  , std::vector<EVENT::TrackerHit*> >  HitMap ; 


CurlKillerProcessor aCurlKillerProcessor ;


CurlKillerProcessor::CurlKillerProcessor() : Processor("CurlKillerProcessor") {
  
  // modify processor description
  _description = "CurlKillerProcessor: Using a 2D(r-phi) histogram, hits from patterns (curlers) traversing the TPC in Z whilst retaining constant r-phi are removed from a new TrackerHit collection " ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "InputCollectionName" , 
			      "Name of the TrackerHit collection"  ,
			      _inputColName ,
			      std::string("TPCTrackerHits") ) ;
  registerProcessorParameter( "CutCollectionName" , 
			      "Name of the cut away TrackerHit collection"  ,
			      _cutColName ,
			      std::string("cutTPCTrackeHits") ) ;
  registerProcessorParameter( "RemainingCollectionName" , 
			      "Name of the remaining TrackerHit collection"  ,
			      _remainingColName ,
			      std::string("remainingTPCTrackerHits") ) ;
  registerProcessorParameter( "BinSize" , 
			      "Bin size in square root of pad multiples"  ,
			      _binSize ,
			      int(2) ) ;
  registerProcessorParameter( "MultiplicityCut" , 
			      "Cut for the number of hits allowed in one bin"  ,
			      _multiplicityCut,
			      int(4) ) ;
  registerProcessorParameter( "PadHeight" , 
			      "TPC PadHeight"  ,
			      _padHeight,
			      float(6.2) ) ;
  registerProcessorParameter( "PadWidth" , 
			      "TPC PadWidth"  ,
			      _padWidth,
			      float(2.2) ) ;
}


void CurlKillerProcessor::init() { 

  // usually a good idea to
  printParameters() ;

  //  MarlinCED::init(this) ;


  
  _nRun = 0 ;
  _nEvt = 0 ;


}

void CurlKillerProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void CurlKillerProcessor::processEvent( LCEvent * evt ) { 

    
  // this gets called for every event 
  // usually the working horse ...

  static bool firstEvent = true ;  

  if(firstEvent==true) cout << "CurlKillerProcessor called for first event" << endl;

  firstEvent = false ;

  // Reset drawing buffer and START drawing collection
  //    MarlinCED::newEvent(this) ;
  
  LCCollection* THcol = 0 ;

  try{
    THcol = evt->getCollection( _inputColName ) ;
  }

  catch(DataNotAvailableException &e){
  }
  
  if( THcol != 0 ){
    
    int n_THits = THcol->getNumberOfElements()  ;

    
    const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
    const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

    const gear::DoubleVec & planeExt = padLayout.getPlaneExtent() ;

    double gearRMin = planeExt[0] ;
    double gearRMax = planeExt[1] ;


    // FIXME:SJA this gives strange values for the padwidth
//     double binHeight = padLayout.getRowHeight(1) * (double) _binSize ;
//     double binWidth = padLayout.getPadWidth(1)  * (double) _binSize ;
    double binHeight = _padHeight * (double) _binSize ;
    double binWidth = _padWidth * (double) _binSize ;
    
    // create bined pad layout
    const gear::FixedPadSizeDiskLayout padsAsBins(gearRMin, gearRMax, binHeight, binWidth) ;                                       
    // create hit map
    HitMap hitMap ; 

    LCCollectionVec* cutCol = new LCCollectionVec( LCIO::TRACKERHIT ) ;
    cutCol->setSubset() ; 

    LCCollectionVec* remainingCol = new LCCollectionVec( LCIO::TRACKERHIT ) ;
    remainingCol->setSubset() ; 

    for(int i=0; i< n_THits; i++){
      
      TrackerHit* THit = dynamic_cast<TrackerHit*>( THcol->getElementAt( i ) ) ;
      
      double *pos;
      pos = (double*) THit->getPosition(); 

      double rad = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);      
      
      //get phi of current hit
      float phi = atan2(pos[1],pos[0]);
      if (phi<0.) phi=phi+twopi;     
      
      int padIndex = padsAsBins.getNearestPad(rad,phi);
      unsigned int iRowHit = padsAsBins.getRowNumber(padIndex);      
      unsigned int iPhiHit = padsAsBins.getPadNumber(padIndex);
      
      // enter hit into hitMap
      hitMap[  make_key( iRowHit, iPhiHit ) ].push_back(  THit ) ;      

    }

    //loop over hitmap and fill both collections of cut and remaining hits

    for( HitMap::iterator it = hitMap.begin() ;it != hitMap.end() ;  ++it ) {
      
      const std::vector<EVENT::TrackerHit*>& v = it->second ;
         
      if(   v.size() >=  (unsigned) _multiplicityCut ) {

        for( unsigned i = 0 ; i < v.size() ; i++){
          cutCol->addElement(  v[i] ) ;
        } 

//                 int color = 0x88ff88 ;
//                 int layer = 6 ;
//                 int marker = 2 ;
//                 int size = 1 ;
        
                //                MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer) ;  
        
      }

      if(   v.size() < (unsigned) _multiplicityCut  ) {
 
        for( unsigned i = 0 ; i < v.size() ; i++){
          remainingCol->addElement(  v[i] ) ;
        }            

//                  int color = 0x88ffff ;
//                  int layer = 7 ;
//                  int marker = 2 ;
//                  int size = 1 ;
        
                 //                 MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer) ; 
      }
      
    }   
    
    evt->addCollection( remainingCol , _remainingColName ) ;
    evt->addCollection( cutCol , _cutColName ) ;
    

  }
  
  //++++++++++++++++++++++++++++++++++++
  //    MarlinCED::draw(this) ;
  //++++++++++++++++++++++++++++++++++++

  
  _nEvt ++ ;
  
}



void CurlKillerProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void CurlKillerProcessor::end(){ 
  
  //   std::cout << "CurlKillerProcessor::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;

}

