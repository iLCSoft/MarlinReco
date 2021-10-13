#ifndef VXDGeometry_h
#define VXDGeometry_h 1

//#include "lcio.h"
#include <map>
#include <vector>
//#include <gsl/gsl_rng.h>

#include <gear/GearMgr.h>
#include "CLHEP/Vector/TwoVector.h"

/** Helper struct for VXD ladder geometry */
struct VXDLadder{
  double phi{};   // phi of ladder - rotation araound z-axis
  gear::Vector3D trans{}; // translation after rotation
//   CLHEP::Hep2Vector p0 ;  // 'left' end of ladder in r-phi
//   CLHEP::Hep2Vector p1 ;  // 'right' end of ladder in r-phi
//   CLHEP::Hep2Vector u  ;  // unit vector along ladder in r-phi
};
typedef std::vector< std::vector< VXDLadder > > VXDLadders ;

/** Helper struct for VXD layer geometry */
struct VXDLayer{
  double rMin{};
  double rMax{};
  double length{};
  double width{};
  double thickness{};
  double gap{};
  double ladderArea{};
  int nLadders{};
};

typedef std::vector< VXDLayer >  VXDLayers ;



/** ======= VXDGeometry ========== <br>
 * Helper class for VXD geomtry transformations: from lab frame to ladder frame  
 * coordinates and inverse. 
 * The ladder reference system has its origin at the middle of the sensitive box in x and y, 
 * where the x-coordinate runs along the thickness of the ladder and t he y coordinate runs along its 
 * width. The z coordinate is the same as in the lab frame. 
 *  
 * <br>
 * @version $Id$
 * @author F.Gaede, DESY
 */
class VXDGeometry  {
  
public:

  VXDGeometry(const VXDGeometry&) = delete;
  VXDGeometry& operator=(const VXDGeometry&) = delete;  

  VXDGeometry(gear::GearMgr* gearMgr) ;

  
  /** Return the pair (layerID, ladderID) for the given position, 
   *  (-1,-1) if not in sensitive volume (in the given layer).
   */
  std::pair<int,int> getLadderID( gear::Vector3D labPos, int layerID=-1 ) ;
  
  /** Convert a position in the lab frame to local ladder coordinates 
   *  (x_ladder==0 is the middle of the sensitive).
   */
  gear::Vector3D lab2LadderPos( gear::Vector3D labPos, int layerID, int ladderID) ;
  
  /** Convert a position in local ladder coordinates 
   *  (x_ladder==0 is the middle of the sensitive) to the lab frame.
   */
  gear::Vector3D ladder2LabPos( gear::Vector3D ladderPos, int layerID, int ladderID) ;

  /** Convert a direction in the lab frame to local ladder coordinates 
   *  (x_ladder==0 is the middle of the sensitive). 
   */
  gear::Vector3D lab2LadderDir( gear::Vector3D labPos, int layerID, int ladderID) ;
  
  /** Convert a direction in local ladder coordinates 
   *  (x_ladder==0 is the middle of the sensitive) to the lab frame.
   */
  gear::Vector3D ladder2LabDir( gear::Vector3D ladderPos, int layerID, int ladderID) ;

  void test() ;


protected:

  void init() ;
  VXDGeometry(){}

  gear::GearMgr* _gearMgr{};
  VXDLadders _vxdLadders{};
  VXDLayers  _vxdLayers{};

} ;

#endif
