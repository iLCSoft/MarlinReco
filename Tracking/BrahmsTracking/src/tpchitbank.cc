#include <iostream>
#include <vector>
//stl exception handler
#include<stdexcept>

#include "tpchitbank.h"

using namespace std;

// Global pointer to the TPC_Hit_Bank structure
TPC_Hit_Bank * TPCHitBank;

// Declaration of Common Block used to store the number of hits in the hit bank defined in TPChitbank.h
CNTPC_DEF CNTPC;

TPC_Hit_Bank::TPC_Hit_Bank()
{
}

TPC_Hit_Bank::~TPC_Hit_Bank()
{
}

//void TPC_Hit_Bank::add_hit(float X, float Y, float Z, float E, int SubID, int TrkID, int PntToEx, int NEx, int ResC, float Res1, float Res2)
void TPC_Hit_Bank::add_hit(float X, float Y, float Z, float E, int SubID, float Res1, float Res2, int TrkID)
{

  TPC_hit ahit;

  ahit.x = X;
  ahit.y = Y;
  ahit.z = Z;
  ahit.energy = E;
  ahit.subdetector_ID = SubID;
  //  ahit.pointer_to_first_exclusion = PntToEx;
  //  ahit.number_of_exclusions = NEx;
  //  ahit.resolution_code = ResC;
  ahit.resolution_1 = Res1;
  ahit.resolution_2 = Res2;
  ahit.track_ID = TrkID;

  hit_bank.push_back(ahit);

}

// due to the fact that the bank will be accessed by integer value of the hit number this is probaly dangerous
// so leave it out for now

//void TPC_Hit_Bank::remove_hit(int hit)
//{
//  hit_bank.erase(hit_bank.begin()+hit);
//}

float readtpchitscpp(int attribute, int hit)
{

  if(hit>TPCHitBank->size()) return 0.;

  hit = hit - 1;

  switch (attribute) {
  case 1: 
    return TPCHitBank->getX(hit);
    break;
  case 2: 
    return TPCHitBank->getY(hit);
     break;
  case 3: 
    return TPCHitBank->getZ(hit);
    break;
  case 4: 
    return TPCHitBank->getEnergy(hit);
    break;
  case 5: 
    return TPCHitBank->getSubdetectorID(hit);
    break;
  case 6: 
    return TPCHitBank->getResolution1(hit);
    break;
  case 7: 
    return TPCHitBank->getResolution2(hit);
    break;
  case 8: 
    return TPCHitBank->getTrackID(hit);
    break;
  default: 
    throw runtime_error("hit attribute not valid");
  } 
}

int writetpccpp(float value, int attribute, int hit)
{

  hit = hit - 1;

  if(hit>TPCHitBank->size()) return 1;

  switch (attribute) {
  case 1: 
    TPCHitBank->setX(value, hit);
    return 0;
    break;
  case 2: 
    TPCHitBank->setY(value, hit);
    return 0;
    break;
  case 3: 
    TPCHitBank->setZ(value, hit);
    return 0;
    break;
  case 4: 
    TPCHitBank->setEnergy(value, hit);
    return 0;
    break;
  case 5: 
    TPCHitBank->setSubdetectorID(((int)value),hit);
    return 0;
    break;
  case 6: 
    TPCHitBank->setResolution1(value,hit);
    return 0;
    break;
  case 7: 
    TPCHitBank->setResolution2(value,hit);
    return 0; 
    break;
  case 8: 
    TPCHitBank->setTrackID(((int)value),hit);
    return 0;
    break;
  default: 
    throw runtime_error("hit attribute not valid");
  } 

}
