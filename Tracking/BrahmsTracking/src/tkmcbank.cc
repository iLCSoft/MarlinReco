#include <iostream>
#include <vector>
//stl exception handler
#include<stdexcept>

#include "tkmcbank.h"

using namespace std;

// Global pointer to the Tk_MC_Bank structure
Tk_MC_Bank * TkMCBank;


Tk_MC_Bank::Tk_MC_Bank()
{
}

Tk_MC_Bank::~Tk_MC_Bank()
{
}



void Tk_MC_Bank::add_mc(float PX, float PY, float PZ, float E, float X, float Y, float Z, int GEANT_PID, int MCTRACK, int NHITS, int HEPEVT_NUM)

{

  tk_mc amc;

  amc.px = PX;
  amc.py = PY;
  amc.pz = PZ;
  amc.energy = E;
  amc.x = X;
  amc.y = Y;
  amc.z = Z;
  amc.geant_pid = GEANT_PID; 
  amc.mctrack = MCTRACK;
  amc.nhits = NHITS;
  amc.hepevt_num = HEPEVT_NUM;

  mc_bank.push_back(amc);

}

// due to the fact that the bank will be accessed by integer value of the mc number this is probaly dangerous
// so leave it out for now

//void Tk_MC_Bank::remove_mc(int mc)
//{
//  mc_bank.erase(mc_bank.begin()+mc);
//}

float readtkmccpp(int attribute, int mc)
{

  if(mc>TkMCBank->size()) return 0.;

  mc = mc - 1;

  switch (attribute) {
  case 1: 
    //    std::cout << "getX = " << TkMCBank->getPX(mc) << std::endl;
    return TkMCBank->getPX(mc);
    break;
  case 2: 
    //    std::cout << "getY = " << TkMCBank->getPY(mc) << std::endl;
    return TkMCBank->getPY(mc);
     break;
  case 3: 
    //    std::cout << "getPZ = " << TkMCBank->getPZ(mc) << std::endl;
    return TkMCBank->getPZ(mc);
    break;
  case 4: 
    //    std::cout << "getEnergy = " << TkMCBank->getEnergy(mc) << std::endl;
    return TkMCBank->getEnergy(mc);
    break;
  case 5: 
    //    std::cout << "getX = " << TkMCBank->getPX(mc) << std::endl;
    return TkMCBank->getPX(mc);
    break;
  case 6: 
    //    std::cout << "getY = " << TkMCBank->getPY(mc) << std::endl;
    return TkMCBank->getPY(mc);
     break;
  case 7: 
    //    std::cout << "getPZ = " << TkMCBank->getPZ(mc) << std::endl;
    return TkMCBank->getPZ(mc);
    break;
  case 8: 
    //    std::cout << "getGEANT_PID = " << TkMCBank->getGEANT_PID(mc) << std::endl;
    return TkMCBank->getGEANT_PID(mc)+0.5;
    break;
  case 9: 
    //    std::cout << "getGEANT_TRACK_NUMBER = " << TkMCBank->getGEANT_TRACK_NUMBER(mc)<< std::endl;
    return TkMCBank->getMCTRACK(mc)+0.5;
    break;
  case 10: 
    //    std::cout << "getNHITS = " << TkMCBank->getNHITS(mc)<< std::endl;
    return TkMCBank->getNHITS(mc);
    break;
  case 11: 
    //    std::cout << "getHEPEVT_NUM = " << TkMCBank->getHEPEVT_NUM(mc) << std::endl;
    return TkMCBank->getHEPEVT_NUM(mc);
    break;
  default: 
    throw runtime_error("mc attribute not valid");
  } 
}

int writetkmccpp(float value, int attribute, int mc)
{

  mc = mc - 1;

  if(mc>TkMCBank->size()) return 1;

  switch (attribute) {
  case 1: 
    TkMCBank->setPX(value, mc);
    return 0;
    break;
  case 2: 
    TkMCBank->setPY(value, mc);
    return 0;
    break;
  case 3: 
    TkMCBank->setPZ(value, mc);
    return 0;
    break;
  case 4: 
    TkMCBank->setEnergy(value, mc);
    return 0;
    break;
  case 5: 
    TkMCBank->setX(value, mc);
    return 0;
    break;
  case 6: 
    TkMCBank->setY(value, mc);
    return 0;
    break;
  case 7: 
    TkMCBank->setZ(value, mc);
    return 0;
    break;
  case 8: 
    TkMCBank->setGEANT_PID(((int)value),mc);
    return 0;
    break;
  case 9: 
    TkMCBank->setMCTRACK(((int)value),mc);
    return 0;
    break;
  case 10: 
    TkMCBank->setNHITS(((int)value),mc);
    return 0;
    break;
  case 11: 
    TkMCBank->setHEPEVT_NUM(((int)value),mc);
    return 0;
    break;
  default: 
    throw runtime_error("mc attribute not valid");
  } 

}
