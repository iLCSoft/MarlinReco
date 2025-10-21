/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "VXDClusterParameters.h"
#include "VXDGeometry.h"

VXDClusterParameters::VXDClusterParameters(const gear::Vector3D& localPos, const gear::Vector3D& a,
                                           const gear::Vector3D& b, int lay, int lad)
    : _localPos(localPos), _cluAxisA(a), _cluAxisB(b), _layerId(lay), _ladderId(lad) {

  // make sure a is the larger axes

  if (_cluAxisA.r() < _cluAxisB.r()) {

    gear::Vector3D tmp = _cluAxisA;
    _cluAxisA = _cluAxisB;
    _cluAxisB = tmp;
  }
}

bool VXDClusterParameters::isPointInClusterEllipse(const gear::Vector3D& pos) {

  // compute focal points of ellipse:

  double a = _cluAxisA.r();
  double b = _cluAxisB.r();
  double epsilon = sqrt(a * a - b * b);

  gear::Vector3D f0 = _localPos + epsilon / a * _cluAxisA;
  gear::Vector3D f1 = _localPos - epsilon / a * _cluAxisA;

  // distance of points from focal points
  gear::Vector3D d0 = f0 - pos;
  gear::Vector3D d1 = f1 - pos;

  // ignore x-coordinate, i.e. project to ladder surface
  d0[0] = 0.0;
  d1[0] = 0.0;

  double dSum = d0.r() + d1.r();

  return dSum <= 2. * a;
}
