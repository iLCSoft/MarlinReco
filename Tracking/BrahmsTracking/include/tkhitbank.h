
#ifndef Tk_Hit_Bank_h
#define Tk_Hit_Bank_h 1

#include<vector>
#include<cfortran.h>
#include<map>
#include<string>
//stl exception handler
#include<stdexcept>

class Tk_Hit_Bank
{

 public:
  Tk_Hit_Bank();
  ~Tk_Hit_Bank();
  void clear();
  void add_hit(float, float, float, float, int, int, int, int, int, float, float);
  //void add_hit(float, float, float, float, int, float, float, int);
  //  void remove_hit(int);
  void setX(float X, int hit){hit_bank[hit].x = X;};
  void setY(float Y, int hit){hit_bank[hit].y = Y;};
  void setZ(float Z, int hit){hit_bank[hit].z = Z;};
  void setEnergy(float E, int hit){hit_bank[hit].energy = E;};
  void setSubdetectorID(int SubID, int hit){hit_bank[hit].subdetector_ID = SubID;};
  void setTrackID(int TrkID, int hit){hit_bank[hit].track_ID = TrkID;};
  void setPntToFirstExclusion(int NEx, int hit){hit_bank[hit].pointer_to_first_exclusion = NEx;};
  void setNExclusion(int PntToEx, int hit){hit_bank[hit].number_of_exclusions = PntToEx;};
  void setResolutionCode(int ResC, int hit){hit_bank[hit].resolution_code = ResC;};
  void setResolution1(float Res1, int hit){hit_bank[hit].resolution_1 = Res1;};
  void setResolution2(float Res2, int hit){hit_bank[hit].resolution_2 = Res2;};



  int   size(){return hit_bank.size();};
  float getX(int i){return hit_bank[i].x;};
  float getY(int i){return hit_bank[i].y;};
  float getZ(int i){return hit_bank[i].z;};
  float getEnergy(int i){return hit_bank[i].energy;};
  int   getSubdetectorID(int i){return hit_bank[i].subdetector_ID;};;
  int   getTrackID(int i){return hit_bank[i].track_ID;};
  int   getPntToFirstExclusion(int i){return hit_bank[i].pointer_to_first_exclusion;};
  int   getNExclusion(int i){return hit_bank[i].number_of_exclusions;};
  int   getResolutionCode(int i){return hit_bank[i].resolution_code;};
  float getResolution1(int i){return hit_bank[i].resolution_1;};
  float getResolution2(int i){return  hit_bank[i].resolution_2;};


  //  omega>0 ? +1 : -1;

  // These methods have no protection against using the wrong key e.g. tpc instead of TPC ...

  void setFirstHitIndex(std::string subdet)
    { if(firstHitEntry.find(subdet) == firstHitEntry.end()) {firstHitEntry[subdet]=hit_bank.size();}
    else { throw std::runtime_error("subdetector does not exit") ;}};

  unsigned int getFirstHitIndex(std::string subdet)
    { if(firstHitEntry.find(subdet) != firstHitEntry.end()) {return firstHitEntry[subdet];}
    else { throw std::runtime_error("subdetector does not exit") ;}}; 
  
  void setLastHitIndex(std::string subdet)
   { if(lastHitEntry.find(subdet) == lastHitEntry.end()) {lastHitEntry[subdet]=hit_bank.size()-1;}
    else { throw std::runtime_error("subdetector does not exit") ;}};

  unsigned int getLastHitIndex(std::string subdet)
  { if(lastHitEntry.find(subdet) != lastHitEntry.end()) {return  lastHitEntry[subdet];}
  else { throw std::runtime_error("subdetector does not exit") ;}};  
  
  int getNumOfSubDetHits(std::string subdet)
    { if(firstHitEntry.find(subdet) == firstHitEntry.end()) {return 0;} 
    else {return 1 + lastHitEntry[subdet] - firstHitEntry[subdet];}};
  

 private:

struct tk_hit 
{
  float    x;
  float    y;
  float    z;
  float    energy;
  int subdetector_ID;
  int track_ID;
  int pointer_to_first_exclusion;
  int number_of_exclusions;
  int resolution_code;
  float resolution_1;
  float resolution_2;

};

 std::vector <tk_hit> hit_bank;
 
 // stores the postion of the first hit from a subdetector
 std::map <const std::string , unsigned int> firstHitEntry;
 // stores the postion of the last hit from a subdetector
 std::map <const std::string , unsigned int> lastHitEntry;

};

// Global pointer to the Tk_Hit_Bank structure which is defined in Tk_Hit_Bank.cc 
extern Tk_Hit_Bank * TkHitBank;


#endif
