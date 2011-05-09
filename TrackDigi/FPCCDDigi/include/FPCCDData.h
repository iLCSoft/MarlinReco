#ifndef FPCCDData_h
#define FPCCDData_h 1

#include <vector>
#include <map>

/** ======= FPCCDData ========== <br>
 * A class to store all FPCCD's data. <br>
 * Used by FPCCDDigitizer and FPCCDClustering 
 *
 * packPixelHits(...) copies PixelHit data in _pxHits(PixelDataBuf_t) 
 * to LCGenericObject of name VTXPixelHits as integer data and 
 * unpackPixelHits(...) copies LCGenericObject to _pxHits.
 *
 * Format of LCGenericObject is as follows:
 *  PixelHits in one ladder are packed in one element.
 *  - First 32 bit word contains layerID and ladderID of the element.
 *    Counting bit position from left to right, layerID and ladderID
 *    information are stored as follows. 
 *      bit 31-16: reserved for future use
 *      bit 15-8 : layerID
 *      bit 7-0  : ladderID
 *  - Last words are for pixel hits in a ladder. Two word is used 
 *    for each pixel hit. 
 *    -- First word
 *       bit 31: 0 
 *       bit 30-29 : hit quality bit
 *               0 = hit is created by a single SimTrackerHit of 
 *                  a signal event.
 *               1 = hit is created by multiple SimTrackerHits of 
 *                  a signal event
 *               2 = hit is created by a sigle SimTrackerHit of 
 *                  a signal event and hits by background events
 *               3 = hit is created by multiple SimTrackerHits 
 *                  of a signal event and background events. 
 *       bit 28-16 : xiID 
 *       bit 15-0  : zetaID
 *    -- Second word
 *       dEdx value given by SimTrackerHit object.
 *       Float value is stored with a help of union.
 *
 * <br>
 * @author Akiya Miyamoto, KEK: 2010-04-19
 * 
 */

namespace IMPL {  class LCCollectionVec; }
namespace EVENT {  class LCCollection; }
class FPCCDPixelHit;

typedef std::map<unsigned int, FPCCDPixelHit*> PixelHitMap_t;
typedef std::vector< std::vector<PixelHitMap_t> > PixelDataBuf_t;


// =================================================================
class FPCCDData {
 protected:
  PixelDataBuf_t _pxHits; // Hit map for each layer/ladder
  int _maxlayer;
  int _maxladder;

 public:
  FPCCDData(int max_layer, int max_ladder);

  // Add
  //  if same hit exists, add ADC values
  //  if new hit, create a new Pixel hit and store in map
  void addPixelHit(FPCCDPixelHit &aHit, bool isSignal);

  // Clear
  void clear();

  // Save _pxHits info in LCCollectionVec
  void packPixelHits(EVENT::LCCollection &col);
  
  // Copy pixelhit info in _pxHits;
  int unpackPixelHits(EVENT::LCCollection &col);

  //
  void Add(FPCCDData &bkgHit);

  void Add(FPCCDData &bkgHit, int layer, int ladder);
  // iterators to get pixel hit.
  PixelHitMap_t::iterator itBegin(int layer, int ladder){ 
	return _pxHits[layer][ladder].begin();
  }
  PixelHitMap_t::iterator itEnd(int layer, int ladder){
	return _pxHits[layer][ladder].end();
  }


  // dump
  void dump();

};

#endif



