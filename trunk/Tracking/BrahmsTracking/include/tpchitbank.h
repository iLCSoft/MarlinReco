// A header file which defines the basics of a tpc

#ifndef TPC_Hit_Bank_h
#define TPC_Hit_Bank_h 1

#include<vector>
#include <cfortran.h>

using namespace std;



class TPC_Hit_Bank
{

 public:
  TPC_Hit_Bank();
  ~TPC_Hit_Bank();
  //  void add_hit(float, float, float, float, int, int, int, int, int, float, float);
  void clear();
  void add_hit(float, float, float, float, int, float, float, int);
  void remove_hit(int);
  void setX(float X, int hit){hit_bank[hit].x = X;};
  void setY(float Y, int hit){hit_bank[hit].y = Y;};
  void setZ(float Z, int hit){hit_bank[hit].z = Z;};
  void setEnergy(float E, int hit){hit_bank[hit].energy = E;};
  void setSubdetectorID(int SubID, int hit){hit_bank[hit].subdetector_ID = SubID;};
  //  void setPntToFirstExclusion(int NEx, int hit){hit_bank[hit].pointer_to_first_exclusion = NEx;};
  //  void getNExclusion(int PntToEx, int hit){hit_bank[hit].number_of_exclusions = PntToEx;};
  //  void setResolutionCode(int ResC, int hit){hit_bank[hit].resolution_code = ResC;};
  void setResolution1(float Res1, int hit){hit_bank[hit].resolution_1 = Res1;};
  void setResolution2(float Res2, int hit){hit_bank[hit].resolution_2 = Res2;};
  void setTrackID(int TrkID, int hit){hit_bank[hit].track_ID = TrkID;};


  int   size(){return hit_bank.size();};
  float getX(int i){return hit_bank[i].x;};
  float getY(int i){return hit_bank[i].y;};
  float getZ(int i){return hit_bank[i].z;};
  float getEnergy(int i){return hit_bank[i].energy;};
  int   getSubdetectorID(int i){return hit_bank[i].subdetector_ID;};
  //  int   getPntToFirstExclusion(int i){return hit_bank[i].pointer_to_first_exclusion;};
  //  int   getNExclusion(int i){return hit_bank[i].number_of_exclusions;};
  //  int   getResolutionCode(int i){return hit_bank[i].resolution_code;};
  float getResolution1(int i){return hit_bank[i].resolution_1;};
  float getResolution2(int i){return  hit_bank[i].resolution_2;};
  int getTrackID(int i){return hit_bank[i].track_ID;};

 private:

struct TPC_hit 
{
  float    x;
  float    y;
  float    z;
  float    energy;
  int subdetector_ID;
  //  int pointer_to_first_exclusion;
  //  int number_of_exclusions;
  //  int resolution_code;
  float resolution_1;
  float resolution_2;
  int track_ID;
};

 vector <TPC_hit> hit_bank;


};

// Global pointer to the TPC_Hit_Bank structure which is defined in TPC_Hit_Bank.cc 
extern TPC_Hit_Bank * TPCHitBank;

// Common Block used to store the number of hits in the hit bank
//typedef struct { int ntphits; } CNTPC_DEF;
//#define  CNTPC COMMON_BLOCK( CNTPC,cntpc)
//COMMON_BLOCK_DEF(CNTPC_DEF,CNTPC);
//extern CNTPC_DEF CNTPC;


#endif
