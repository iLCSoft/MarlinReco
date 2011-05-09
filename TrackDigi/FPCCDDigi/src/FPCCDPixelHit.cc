/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "FPCCDPixelHit.h"
#include <iostream>
#include <assert.h>
#include <EVENT/MCParticle.h>

#define SHIFT_LADDER 0  // 5 bits
#define SHIFT_LAYER  5  // 3 bits
#define SHIFT_ZETA   0  // 16 bits
#define SHIFT_XI    16  // 13 bits 

#define MASK_LADDER  0x0000001F
#define MASK_LAYER   0x000000E0
#define MASK_ZETA    0x0000FFFF
#define MASK_XI      0x1FFF0000

#define FLAG_ADCWORD 0x80000000
#define FLAG_MCPART  0xC0000000
#define FPCCD_ADC_MAX      0x00003FFF
#define FPCCD_QUALITY      0x0000C000
#define MASK_QUALITY 0xC000
#define SHIFT_QUALITY 14

#define MCPARTID_MAX  16384 // 2^14      

// =====================================================================
FPCCDPixelHit::FPCCDPixelHit(unsigned short int layerID, unsigned short int ladderID,
		unsigned short int xiID,    unsigned short int zetaID,
		float edep,     HitQuality_t quality,
		EVENT::MCParticle *mc)
{
  // check validity of input variables
  assert( layerID < 8 );
  assert( ladderID < 32 );
  assert( xiID < 8192 );
  
  _layerID=layerID;
  _ladderID=ladderID;
  _xiID=xiID;
  _zetaID=zetaID;
  _edep=edep;
  _quality=quality;

  if( mc != 0 ) { _MCParticleVec.push_back(mc); }

}

// =====================================================================
void FPCCDPixelHit::print()
{
  std::cout << "ly=" << _layerID << " ld=" << _ladderID 
	<< " zeta=" << _zetaID << " xi=" << _xiID 
       << " edep=" << _edep << " qual.=" << _quality;
  std::cout << std::endl;

}

// =====================================================================
void FPCCDPixelHit::addPixelHit(FPCCDPixelHit &aHit, HitQuality_t addedQuality)
{

  _edep += aHit.getEdep();

  if( _quality == kSingle || _quality == kSignalOverlap ){
    if( addedQuality == kSingle || addedQuality == kSignalOverlap ) { _quality=kSignalOverlap; }
    else { _quality=kBKGOverlap; }
  }
  else if( _quality == kBKG){
    if( addedQuality == kSingle || addedQuality == kSignalOverlap ) { _quality = kBKGOverlap;}
  }
//   if( aHit.getMCParticleVec().size() > 0 ) {
    
//     for(int i=0;i<aHit.getMCParticleVec().size();i++) {
//       _MCParticleVec.push_back(aHit.getMCParticleVec()[i]);
//     }
//  }
}

// =====================================================================
unsigned int FPCCDPixelHit::encodeCellWord()
{
  assert(_xiID < 8192);
  unsigned int ic=( ( ((unsigned int)_xiID)<<SHIFT_XI & MASK_XI ) |
                    ( ((unsigned int)_zetaID)<<SHIFT_ZETA & MASK_ZETA )  );
  return ic;

}

// ====================================================================
void FPCCDPixelHit::decodeCellWord(unsigned int iword)
{
  _zetaID = (unsigned short int) ( (iword & MASK_ZETA ) >> SHIFT_ZETA );
  _xiID   = (unsigned short int) ( (iword & MASK_XI ) >> SHIFT_XI ) ;
}
