

#include <iostream>
#include "voxel.h"

using namespace std;

Voxel_tpc::Voxel_tpc(){
}
Voxel_tpc::Voxel_tpc(int row, int phi, int z, double pos[3], double posRPhi[2], double dedx, double RPhiRes, double ZRes)
{
  _row_index = row;
  _phi_index = phi;
  _z_index = z;
  _coord.setX(pos[0]);
  _coord.setY(pos[1]);
  _coord.setZ(pos[2]);
  _dE_dx = dedx;
  _rPhiRes = RPhiRes;
  _zRes = ZRes;
  _isMerged = false;
  _isClusterHit = false;
}

Voxel_tpc::Voxel_tpc(int row, int phi, int z, CLHEP::Hep3Vector coord, double dedx, double RPhiRes, double ZRes)
{
  _row_index = row;
  _phi_index = phi;
  _z_index = z;
  _coord=coord;
  _dE_dx = dedx;
  _rPhiRes = RPhiRes;
  _zRes = ZRes;
  _isMerged = false;
  _isClusterHit = false;
}

Voxel_tpc::~Voxel_tpc()
{
}

//bool Voxel_tpc::compare_phi( Voxel_tpc * & a, Voxel_tpc * & b){

//  return ( a->getPhiIndex() < b->getPhiIndex());
  
//}

int Voxel_tpc::clusterFind(vector <Voxel_tpc*>* hitList){
  
  if(!this->IsClusterHit()){
    hitList->push_back(this);
    this->setIsClusterHit();
    for(int i=0; i<this->getNumberOfAdjacent();++i){
      int clusterSize = getAdjacent(i)->clusterFind(hitList);
    }
  }

  return hitList->size();
}

