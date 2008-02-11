// A header file which defines a voxel class for the TPC
#include <vector>

using namespace std;

class Voxel_tpc{

 public:
  Voxel_tpc();  
  // the intialation in the constructor here would be preferable though I don't know how to intialise
  // the array xyz[3] here with pos[3], for the mean time the constructor will be put in the .cc file
  //  Voxel_tpc(int row, int phi, int z, double pos[3]) : row_index(row), phi_index(phi), z_index(z){}
  Voxel_tpc(int row, int phi, int z, double pos[3], double posRPhi[2], double dedx, double rPhiRes, double zRes);
  ~Voxel_tpc();
  void setAdjacent(Voxel_tpc *); 
  int getRowIndex() {return row_index;};
  int getPhiIndex() {return phi_index;};
  int getZIndex() {return z_index;};
  Voxel_tpc * getFirstAdjacent();
  int getNumberOfAdjacent(); 
  double getX() {return xyz[0];};
  double getY() {return xyz[1];};
  double getZ() {return xyz[2];};
  double getR() {return rphi[0];};
  double getPhi() {return rphi[1];};
  double getdEdx() {return dE_dx;};
  double getRPhiRes() {return rPhiRes;};
  double getZRes() {return zRes;};
  //  bool compare_phi( Voxel_tpc * & a, Voxel_tpc * & b);



 private:
  int row_index; 
  int phi_index;
  int z_index;
  vector <Voxel_tpc *> adjacent_voxel;
  double xyz[3];   
  double rphi[2];   
  double dE_dx;
  double rPhiRes;
  double zRes;
};
