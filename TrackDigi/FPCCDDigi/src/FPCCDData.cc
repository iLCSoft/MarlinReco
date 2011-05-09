/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "FPCCDPixelHit.h"
#include "FPCCDData.h"

#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/LCTOOLS.h>
#include <iostream>

// =====================================================================
FPCCDData::FPCCDData(int maxlayer, int maxladder): _maxlayer(maxlayer), _maxladder(maxladder)
{
  _pxHits.resize(_maxlayer);
  for(int i=0;i<_maxlayer;i++){
    _pxHits[i].resize(_maxladder);
  }
}

// ===================================================================
void FPCCDData::clear()
{
  for(int layer=0;layer<_maxlayer;layer++){
    for(int ladder=0;ladder<_maxladder;ladder++){
            PixelHitMap_t::iterator it = _pxHits[layer][ladder].begin();
      while( it != _pxHits[layer][ladder].end()){
        delete (*it).second;
        it++ ;
      }
      _pxHits[layer][ladder].clear();
    }
  }
}

// ===================================================================
void FPCCDData::dump()
{
  std::cout << "FPCCDData::dump() prints out all Pixel hits" << std::endl;
  for(int layer=0;layer<_maxlayer;layer++) {
    for(int ladder=0;ladder<_maxladder;ladder++) {
      if( _pxHits[layer][ladder].size() > 0 ) {
        PixelHitMap_t::iterator it=_pxHits[layer][ladder].begin();
        while( it != _pxHits[layer][ladder].end() ) {
          (*it).second->print();
          it++;
        }
      }
    }
  }
}

// ====================================================================
void FPCCDData::addPixelHit(FPCCDPixelHit &aHit, bool isSignal)
{
  int layer=aHit.getLayerID();
  int ladder=aHit.getLadderID();

  unsigned int hitid=aHit.encodeCellWord();
  FPCCDPixelHit::HitQuality_t addedQuality;
  { isSignal ? addedQuality = FPCCDPixelHit::kSingle : addedQuality = FPCCDPixelHit::kBKG ; }  

  PixelHitMap_t::iterator it=_pxHits[layer][ladder].find(hitid);
  // Append ADC if same pixelID exist in the map
  if( it != _pxHits[layer][ladder].end() ) {
    FPCCDPixelHit *findhit=(*it).second;  
    //    addedQuality = aHit.getQuality();
    findhit->addPixelHit(aHit, addedQuality);
     }
  // Add new hit if not exist in the map
  else {
    aHit.setQuality(addedQuality);
    _pxHits[layer][ladder].insert(PixelHitMap_t::value_type(hitid, new FPCCDPixelHit(aHit)));
  }
}

// ===================================================================
void FPCCDData::packPixelHits(EVENT::LCCollection &colvec)
{
  // This function writes all CCD data into hitvec.
  // _pxHits are cleared after copying to save space
  for(int layer=0;layer<_maxlayer;layer++) {
    for(int ladder=0;ladder<_maxladder;ladder++) {
      PixelHitMap_t::iterator it=_pxHits[layer][ladder].begin();
      if( it != _pxHits[layer][ladder].end() ) {
        IMPL::LCGenericObjectImpl *out=new IMPL::LCGenericObjectImpl();
        unsigned int index=0;
        // store layer number, ladder number in the first word
        // Most significant 4 bits (left 4 bits) should be kept unused 
        // for future use as a data format type
        int word0=( (layer & 0x000000FF ) << 8 | ( ladder & 0x000000FF ) );
        out->setIntVal(index, word0 );
        index++;
        while( it !=_pxHits[layer][ladder].end() ) {
          FPCCDPixelHit *aHit=(*it).second;
          // First word, from left to right,
          // MSB =0 to indicate hitID word
          //     next 2 bit for quality
          //     next 13 + 16 bit for hit id ( xi and zeta )
          unsigned int hitid=(unsigned int)aHit->encodeCellWord();
          unsigned int quality=(unsigned int)aHit->getQuality();
          int hitwd= ( ( (unsigned int)quality << 29 & 0x60000000 ) |
                       ( (hitid & 0x1FFFFFFF ) ) );
          out->setIntVal(index++, hitwd );
          // 2nd word is edep
          union intfloat { float edep; int iedep; } edepout;

          edepout.edep=aHit->getEdep();
          out->setIntVal(index++, edepout.iedep);
          it++;
        } // moving data in aHit to LCGenericObjectImpl
        colvec.addElement(out);  // add one element
      }
    }
  }
  clear();
}

// =====================================================
int FPCCDData::unpackPixelHits(EVENT::LCCollection &col)
{
// Convert pixelhit data in col to _pxHits 
// returns total number of PixelHits converted.

  int nhits=0;
  int nElements=col.getNumberOfElements();

//   UTIL::LCTOOLS::printLCGenericObjects(col);

  clear();
  for(int ie=0;ie<nElements;ie++){
    int ig=0;
    EVENT::LCGenericObject *obj=dynamic_cast<EVENT::LCGenericObject*>(col.getElementAt(ie));
    unsigned int iw0=obj->getIntVal(ig);
    int layer=( ( iw0>>8 ) & 0x000000FF );
    int ladder= ( iw0 & 0x000000FF ) ;
    ig++;
    while( ig < obj->getNInt() ) {
      unsigned int iw=obj->getIntVal(ig);
      FPCCDPixelHit *aHit=new FPCCDPixelHit(layer, ladder);
      // Converting first word
      aHit->decodeCellWord(iw);
      unsigned int qwd=(iw & 0x60000000 ) >> 29 ;
      switch ( qwd ) {
      case 0 : aHit->setQuality(FPCCDPixelHit::kSingle); break;
      case 1 : aHit->setQuality(FPCCDPixelHit::kSignalOverlap) ; break ;
      case 2 : aHit->setQuality(FPCCDPixelHit::kBKGOverlap) ; break ;
      case 3 : aHit->setQuality((FPCCDPixelHit::HitQuality_t)(FPCCDPixelHit::kSignalOverlap|FPCCDPixelHit::kBKGOverlap)) ; break ;
      }
      ig++;
      // converting second word ; 
      union intfloat { float edep ; int iedep ; } edepin ;
      edepin.iedep=obj->getIntVal(ig);
      aHit->setEdep(edepin.edep);
      ig++;
      unsigned int hitid=aHit->encodeCellWord();
      _pxHits[layer][ladder].insert(PixelHitMap_t::value_type(hitid, aHit));
      nhits++;
//      aHit->print();
    }
  }
  return nhits;
}

// =====================================================
void FPCCDData::Add(FPCCDData &bgHit)
{
  FPCCDPixelHit *aHit;
  
  for(int i=0; i<_maxlayer; i++){  
    for(int j=0; j<_maxladder; j++){
      
      
      PixelHitMap_t::iterator it= bgHit.itBegin(i, j);
      while( it != bgHit.itEnd(i, j)){
        
        aHit=(*it).second;
        
        FPCCDData::addPixelHit( *aHit, false);
        
        it++;
      }
      
    }
  }
}
