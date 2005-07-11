#include <iostream>
#include "voxel.h"

using namespace std;

Voxel_tpc::Voxel_tpc(){
}

Voxel_tpc::Voxel_tpc(int row, int phi, int z, double pos[3], double dedx) 
{
  row_index = row;
  phi_index = phi;
  z_index = z;
  xyz[0] = pos[0];
  xyz[1] = pos[1];
  xyz[2] = pos[2];
  dE_dx = dedx;
}

Voxel_tpc::~Voxel_tpc()
{
}

void Voxel_tpc::setAdjacent(Voxel_tpc * p_voxel)
{
  adjacent_voxel.push_back(p_voxel);
}

Voxel_tpc * Voxel_tpc::getFirstAdjacent()
{
  return *(adjacent_voxel.begin());
}

int  Voxel_tpc::getNumberOfAdjacent()
{
  return adjacent_voxel.size();
}

//bool Voxel_tpc::compare_phi( Voxel_tpc * & a, Voxel_tpc * & b){

//  return ( a->getPhiIndex() < b->getPhiIndex());
  
//}
