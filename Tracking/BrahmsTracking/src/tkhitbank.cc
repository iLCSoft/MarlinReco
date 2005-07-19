#include <iostream>
#include <vector>
#include <string>
//stl exception handler
#include<stdexcept>

#include "tkhitbank.h"

using namespace std;

// Global pointer to the Tk_Hit_Bank structure
Tk_Hit_Bank * TkHitBank;


Tk_Hit_Bank::Tk_Hit_Bank()
{
}

Tk_Hit_Bank::~Tk_Hit_Bank()
{
}

void Tk_Hit_Bank::add_hit(float X, float Y, float Z, float E, int SubID,  int TrkID, int PntToEx, int NEx, int ResC, float Res1, float Res2)
  //void Tk_Hit_Bank::add_hit(float X, float Y, float Z, float E, int SubID, float Res1, float Res2, int TrkID)
{

  tk_hit ahit;

  ahit.x = X;
  ahit.y = Y;
  ahit.z = Z;
  ahit.energy = E;
  ahit.subdetector_ID = SubID;
  ahit.track_ID = TrkID;
  ahit.pointer_to_first_exclusion = PntToEx;
  ahit.number_of_exclusions = NEx;
  ahit.resolution_code = ResC;
  ahit.resolution_1 = Res1;
  ahit.resolution_2 = Res2;


  hit_bank.push_back(ahit);

}

// due to the fact that the bank will be accessed by integer value of the hit number this is probaly dangerous
// so leave it out for now

//void Tk_Hit_Bank::remove_hit(int hit)
//{
//  hit_bank.erase(hit_bank.begin()+hit);
//}

int subdetfirsthitindex(string subdet)
{
  return TkHitBank->getFirstHitIndex( subdet )+1;
}

int numofsubdethits(string subdet)
{
  return TkHitBank->getNumOfSubDetHits( subdet );
}

float readtkhitscpp(int attribute, int hit)
{

  if(hit>TkHitBank->size()) return 0.;

  hit = hit - 1;

  switch (attribute) {
  case 1: 
    //    std::cout << "getX = " << TkHitBank->getX(hit) << std::endl;
    return TkHitBank->getX(hit);
    break;
  case 2: 
    //    std::cout << "getY = " << TkHitBank->getY(hit) << std::endl;
    return TkHitBank->getY(hit);
     break;
  case 3: 
    //    std::cout << "getZ = " << TkHitBank->getZ(hit) << std::endl;
    return TkHitBank->getZ(hit);
    break;
  case 4: 
    //    std::cout << "getEnergy = " << TkHitBank->getEnergy(hit) << std::endl;
    return TkHitBank->getEnergy(hit);
    break;
  case 5: 
    //    std::cout << "getSubdetectorID = " << TkHitBank->getSubdetectorID(hit) << std::endl;
    return TkHitBank->getSubdetectorID(hit)+0.5;
    break;
  case 6: 
    //    std::cout << "getTrackID = " << TkHitBank->getTrackID(hit)<< std::endl;
    return TkHitBank->getTrackID(hit)+0.5;
    break;
  case 7: 
    //    std::cout << "getPntToFirstExclusion = " << TkHitBank->getPntToFirstExclusion(hit) << std::endl;
    return TkHitBank->getPntToFirstExclusion(hit)+0.5;
    break;
  case 8: 
    //    std::cout << "TkHitBank->getNExclusion = " << TkHitBank->getNExclusion(hit)<< std::endl;
    return TkHitBank->getNExclusion(hit)+0.5;
    break;
  case 9: 
    //    std::cout << "getResolutionCode = " << TkHitBank->getResolutionCode(hit) << std::endl;
    return TkHitBank->getResolutionCode(hit)+0.5;
    break;
  case 10: 
    //    std::cout << "getResolution1 = " << TkHitBank->getResolution1(hit)<< std::endl;
    return TkHitBank->getResolution1(hit);
    break;
  case 11: 
    //    std::cout << "getResolution2 = " << TkHitBank->getResolution2(hit) << std::endl;
    return TkHitBank->getResolution2(hit);
    break;
  default: 
    throw runtime_error("hit attribute not valid");
  } 
}

int writetkhitcpp(float value, int attribute, int hit)
{

  hit = hit - 1;

  if(hit>TkHitBank->size()) return 1;

  switch (attribute) {
  case 1: 
    TkHitBank->setX(value, hit);
    return 0;
    break;
  case 2: 
    TkHitBank->setY(value, hit);
    return 0;
    break;
  case 3: 
    TkHitBank->setZ(value, hit);
    return 0;
    break;
  case 4: 
    TkHitBank->setEnergy(value, hit);
    return 0;
    break;
  case 5: 
    TkHitBank->setSubdetectorID(((int)value),hit);
    return 0;
    break;
  case 6: 
    TkHitBank->setTrackID(((int)value),hit);
    return 0;
    break;
  case 7: 
    TkHitBank->setPntToFirstExclusion(((int)value),hit);
    return 0;
    break;
  case 8: 
    TkHitBank->setNExclusion(((int)value),hit);
    return 0;
    break;
  case 9: 
    TkHitBank->setResolutionCode(((int)value),hit);
    return 0;
    break;
  case 10: 
    TkHitBank->setResolution1(value,hit);
    return 0;
    break;
  case 11: 
    TkHitBank->setResolution2(value,hit);
    return 0;
    break;
  default: 
    throw runtime_error("hit attribute not valid");
  } 

}
