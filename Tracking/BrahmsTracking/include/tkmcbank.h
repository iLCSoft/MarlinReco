#ifndef Tk_MC_Bank_h
#define Tk_MC_Bank_h 1

#include<vector>
#include <cfortran.h>



class Tk_MC_Bank
{

 public:
  Tk_MC_Bank();
  ~Tk_MC_Bank();
  void clear();
  // this class uses PTRTYPE which is defined in cpointer.h as of type long
  void add_mc(float, float, float, float, float, float, float, int, int, int, int);
  void remove_mc(int);
  void setPX(float PX, int mc){mc_bank[mc].px = PX;};
  void setPY(float PY, int mc){mc_bank[mc].py = PY;};
  void setPZ(float PZ, int mc){mc_bank[mc].pz = PZ;};
  void setEnergy(float E, int mc){mc_bank[mc].energy = E;};
  void setX(float X, int mc){mc_bank[mc].x = X;};
  void setY(float Y, int mc){mc_bank[mc].y = Y;};
  void setZ(float Z, int mc){mc_bank[mc].z = Z;};
  void setGEANT_PID(int GEANT_PID, int mc){mc_bank[mc].geant_pid = GEANT_PID;};
  void setMCTRACK(int MCTRACK, int mc){mc_bank[mc].mctrack = MCTRACK;};
  void setNHITS(int NHITS, int mc){mc_bank[mc].nhits = NHITS;};
  void setHEPEVT_NUM(int HEPEVT_NUM, int mc){mc_bank[mc].hepevt_num = HEPEVT_NUM;};


  int   size(){return mc_bank.size();};
  float getPX(int i){return mc_bank[i].px;};
  float getPY(int i){return mc_bank[i].py;};
  float getPZ(int i){return mc_bank[i].pz;};
  float getEnergy(int i){return mc_bank[i].energy;};
  float getX(int i){return mc_bank[i].x;};
  float getY(int i){return mc_bank[i].y;};
  float getZ(int i){return mc_bank[i].z;};
  int   getGEANT_PID(int i){return mc_bank[i].geant_pid;};
  int   getMCTRACK(int i){return mc_bank[i].mctrack;};
  int   getNHITS(int i){return mc_bank[i].nhits;};
  int   getHEPEVT_NUM(int i){return mc_bank[i].hepevt_num;};

 private:

struct tk_mc 
{
  float    px;
  float    py;
  float    pz;
  float    energy;
  float    x;
  float    y;
  float    z;
  int geant_pid;
  int mctrack;
  int nhits;
  int hepevt_num;
};

 std::vector <tk_mc> mc_bank;


};

// Global pointer to the Tk_MC_Bank structure which is defined in Tk_MC_Bank.cc 
extern Tk_MC_Bank * TkMCBank;


#endif
