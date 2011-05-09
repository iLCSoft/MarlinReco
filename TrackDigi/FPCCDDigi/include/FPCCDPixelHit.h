#ifndef FPCCDPixelHit_h
#define FPCCDPixelHit_h 1

#include <vector>

/** ======= FPCCDPixelhit ========== <br>
 * A class represents one FPCCD pixel hit <br>
 * PixelHit ID is given by ( layer, ladder, xi, zeta ). xi and zeta are 
 * pixel address in a coordinate system local to a ladder. zeta is along 
 * Z axis, xi is in the ladder plain and eta is perpendicular to the plain 
 * of ladder. 
 *
 *  cellID convension for FPCCD
 *  layer ID is from 0 to 5 ( for 6 layers configuration )
 *  ladder ID is from 0 to ( number of ladders - 1 ), depending on the layer ID.
 *  Pixel IDs start from 0.
 *  (xi, zeta) is a local coordinate attached to a ladder.
 *  xi is coordinate along shorter axix of ladder
 *  zeta is coordinate along longer axis of ladder.
 *  the direction of xi axis is same as X axis of laboratory coordinate when 
 * a ladder is at phi=90 degree.
 *  the direction of zeta axis is same as Z axis of laboratory coordinate.
 *  the eta axis is perpendicular to a ladder.  when phi=90 degree, 
 * it matched Y axis.
 *  The origin of local coordinate system, (xi, eta, zeta) is placed at 
 * the intersection of the line from IP and perpendicular to the ladder.
 *
 * See also comment in FPCCDDigitizer::encodeFPCCDID(...)
 *
 * <br>
 * @author Akiya Miyamoto, KEK: 2010-04-19
 * 
 */

namespace EVENT {
  class MCParticle;
}

// =================================================================
class FPCCDPixelHit 
{
 public:
  typedef enum { kSingle=0, kSignalOverlap=0x01, kBKGOverlap=0x02, kBKG=0x03} HitQuality_t;

 protected:
  unsigned short int _layerID;
  unsigned short int _ladderID;
  unsigned short int _xiID;
  unsigned short int _zetaID;
  float              _edep;
  HitQuality_t _quality;
  std::vector<EVENT::MCParticle*> _MCParticleVec;

 public:
  FPCCDPixelHit(unsigned short int layerID=0, unsigned short int ladderID=0,
		unsigned short int xiID=0,    unsigned short int zetaID=0,
		float edep=0.0,     HitQuality_t quality=kSingle,
		EVENT::MCParticle *mc=0);

  void setLayerID(int layerid){ _layerID=layerid; }
  void setLadderID(int ladderid){ _ladderID=ladderid;}
  void setXiID(int xiid){ _xiID=xiid; }
  void setZetaID(int zetaid) { _zetaID=zetaid; }
  void setEdep(float edep){ _edep=edep; }
  void setQuality(HitQuality_t quality){ _quality=quality; }
  
  int getLayerID(){ return _layerID; }
  int getLadderID(){ return _ladderID; }
  int getXiID(){ return _xiID; }
  int getZetaID(){ return _zetaID; }
  float getEdep(){ return _edep; }
  HitQuality_t getQuality(){ return _quality; }
  std::vector<EVENT::MCParticle*> getMCParticleVec(){ return _MCParticleVec; }
  EVENT::MCParticle *getMCParticle(int index){ return _MCParticleVec[index]; }
  int getNMCParticles(){ return _MCParticleVec.size(); }

  // add pixel hit
  void addPixelHit(FPCCDPixelHit &aHit, HitQuality_t addedQuality);
  
  void print();
  unsigned int encodeCellWord();
  void decodeCellWord(unsigned int word);

};

#endif
