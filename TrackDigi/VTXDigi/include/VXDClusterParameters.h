#ifndef VXDClusterParameters_h
#define VXDClusterParameters_h 1

#include "LCRTRelations.h"
#include "gear/GEAR.h"
#include "lcio.h"

/** ======= VXDClusterParameters ========== <br>
 * Holds cluster parameters for a VXD hit - to be attached at runtime to the SimTrackerHit.
 * Cluster is defined by two  (principal component) axis and the corresponding extensions (eigenvalues).
 *
 * @version $Id$
 * @author F.Gaede, DESY
 */

class VXDClusterParameters {

public:
  /** The standard constructor: takes the position of the hit in ladder coordinates
   *  and the cluster axes in  ladder coordinates. Also store the layerId and ladderId for saving cpu
   *  time later.
   */
  VXDClusterParameters(const gear::Vector3D& localPos, const gear::Vector3D& a, const gear::Vector3D& b, int lay,
                       int lad);

  gear::Vector3D getClusterPosition() { return _localPos; }
  gear::Vector3D getClusterAxisA() { return _cluAxisA; }
  gear::Vector3D getClusterAxisB() { return _cluAxisB; }

  int getLayerId() { return _layerId; }
  int getLadderId() { return _ladderId; }

  /** True, if the 2D projection to the ladder surface of the given position is within the ellipse
   *  defined by the two cluster axes.
   */
  bool isPointInClusterEllipse(const gear::Vector3D& pos);

protected:
  VXDClusterParameters() : _localPos(0, 0, 0), _cluAxisA(0, 0, 0), _cluAxisB(0, 0, 0), _layerId(-1), _ladderId(-1) {}

  gear::Vector3D _localPos{}; // position of hit on ladder
  gear::Vector3D _cluAxisA{};
  gear::Vector3D _cluAxisB{};
  int _layerId{};
  int _ladderId{};
};

/** Allows to use VXDClusterParameters as runtime extension object.
 */
struct ClusterParams : public lcio::LCOwnedExtension<ClusterParams, VXDClusterParameters> {};

#endif
