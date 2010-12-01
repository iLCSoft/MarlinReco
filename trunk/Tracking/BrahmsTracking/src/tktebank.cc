#include <iostream>
#include <vector>
//stl exception handler
#include<stdexcept>

#include "tktebank.h"

using namespace std;

// Global pointer to the Tk_Te_Bank structure
Tk_Te_Bank * TkTeBank;


Tk_Te_Bank::Tk_Te_Bank()
{
}

Tk_Te_Bank::~Tk_Te_Bank()
{
}

void Tk_Te_Bank::clear(){
  te_bank.clear();
  itedat_bank.clear();
}

void Tk_Te_Bank::add_te(int subid,int submod,int unused,int MesrCode,int PnteTE,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float cov1,float cov2,float cov3,float cov4,float cov5,float cov6,float cov7,float cov8,float cov9,float cov10,float cov11,float cov12,float cov13,float cov14,float cov15)
{

  tk_te ate;
  ate.subdetector_ID = subid;
  ate.submodule = submod;
  ate.unused = unused;
  ate.measurement_code = MesrCode;
  ate.pointer_to_end_of_TE = PnteTE;
  ate.charge = Q;
  ate.ndf = ndf;
  ate.chi2 = chi2;
  ate.length = L;
  ate.coord1_of_ref_point = cord1;
  ate.coord2_of_ref_point = cord2;
  ate.coord3_of_ref_point = cord3;
  ate.theta = theta;
  ate.phi = phi;
  ate.invp = invp;
  ate.de_dx = dedx;
  ate.covmatrix1 = cov1;
  ate.covmatrix2 = cov2;
  ate.covmatrix3 = cov3;
  ate.covmatrix4 = cov4;
  ate.covmatrix5 = cov5;
  ate.covmatrix6 = cov6;
  ate.covmatrix7 = cov7;
  ate.covmatrix8 = cov8;
  ate.covmatrix9 = cov9;
  ate.covmatrix10 = cov10;
  ate.covmatrix11 = cov11;
  ate.covmatrix12 = cov12;
  ate.covmatrix13 = cov13;
  ate.covmatrix14 = cov14;
  ate.covmatrix15 = cov15;


  te_bank.push_back(ate);

  unsigned int size = itedat_bank.size() ;
  itedat_bank.resize(size+1) ;

}

// due to the fact that the bank will be accessed by integer value of the te number this is probaly dangerous
// so leave it out for now

//void Tk_Te_Bank::remove_te(int te)

//{
//  te_bank.erase(te_bank.begin()+te);
//}

//int tkmktecpp(int subid,int submod,int unused,int MesrCode,int PnteTE,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float cov1,float cov2,float cov3,float cov4,float cov5,float cov6,float cov7,float cov8,float cov9,float cov10,float cov11,float cov12,float cov13,float cov14,float cov15)
//{
// 
//  TkTeBank->add_te(subid, submod, unused, MesrCode, PnteTE, Q, ndf, chi2, L, cord1, cord2, cord3, theta, phi, invp, dedx, cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8, cov9, cov10, cov11, cov12, cov13, cov14, cov15);
// 
//  return 0;
//}

int tkmktecpp(int subid,int submod,int unused,int MesrCode,int PnteTE,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float* cov)
{

  float cov1 = cov[0]; 
  float cov2 = cov[1];
  float cov3 = cov[2];
  float cov4 = cov[3];
  float cov5 = cov[4];
  float cov6 = cov[5];
  float cov7 = cov[6];
  float cov8 = cov[7];
  float cov9 = cov[8];
  float cov10 = cov[9];
  float cov11 = cov[10];
  float cov12 = cov[11];
  float cov13 = cov[12];
  float cov14 = cov[13];
  float cov15 = cov[14];

//   cout << "cov1 = " << cov1 << endl; 
//   cout << "cov2 = " << cov2 << endl; 
//   cout << "cov3 = " << cov3 << endl; 
//   cout << "cov4 = " << cov4 << endl; 
//   cout << "cov5 = " << cov5 << endl; 
//   cout << "cov6 = " << cov6 << endl; 
//   cout << "cov15 = " << cov15 << endl; 

  TkTeBank->add_te(subid, submod, unused, MesrCode, PnteTE, Q, ndf, chi2, L, cord1, cord2, cord3, theta, phi, invp, dedx, cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8, cov9, cov10, cov11, cov12, cov13, cov14, cov15);
  return 0;
}

float rreadtktecpp(int attribute, int te)
{

  te = te - 1;

  if(te>TkTeBank->size()) return 0.;

  switch (attribute) {
  case 1: 
    return TkTeBank->getSubdetector_ID(te);
    break;
  case 2: 
    return TkTeBank->getSubmodule(te);
     break;
  case 3: 
    return TkTeBank->getUnused(te);
    break;
  case 4: 
    return TkTeBank->getMeasurement_code(te);
    break;
  case 5: 
    return TkTeBank->getPointer_to_end_of_TE(te);
    break;
  case 6: 
    return TkTeBank->getCharge(te);
    break;
  case 7: 
    return TkTeBank->getNdf(te);
    break;
  case 8: 
    return TkTeBank->getChi2(te);
    break;
  case 9: 
    return TkTeBank->getLength(te);
    break;
  case 10: 
    return TkTeBank->getCoord1_of_ref_point(te);
    break;
  case 11: 
    return TkTeBank->getCoord2_of_ref_point(te);
    break;
  case 12: 
    return TkTeBank->getCoord3_of_ref_point(te);
    break;
  case 13: 
    return TkTeBank->getTheta(te);
    break;
  case 14: 
    return TkTeBank->getPhi(te);
    break;
  case 15: 
    return TkTeBank->getInvp(te);
    break;
  case 16: 
    return TkTeBank->getDe_dx(te);
    break;
  case 17: 
    return TkTeBank->getCovmatrix1(te);
    break;
  case 18: 
    return TkTeBank->getCovmatrix2(te);
    break;
  case 19: 
    return TkTeBank->getCovmatrix3(te);
    break;
  case 20: 
    return TkTeBank->getCovmatrix4(te);
    break;
  case 21: 
    return TkTeBank->getCovmatrix5(te);
    break;
  case 22: 
    return TkTeBank->getCovmatrix6(te);
    break;
  case 23: 
    return TkTeBank->getCovmatrix7(te);
    break;
  case 24: 
    return TkTeBank->getCovmatrix8(te);
    break;
  case 25: 
    return TkTeBank->getCovmatrix9(te);
    break;
  case 26: 
    return TkTeBank->getCovmatrix10(te);
    break;
  case 27: 
    return TkTeBank->getCovmatrix11(te);
    break;
  case 28: 
    return TkTeBank->getCovmatrix12(te);
    break;
  case 29: 
    return TkTeBank->getCovmatrix13(te);
    break;
  case 30: 
    return TkTeBank->getCovmatrix14(te);
    break;
  case 31: 
    return TkTeBank->getCovmatrix15(te);
    break;
  default: 
    throw runtime_error("te attribute not valid");
  } 
}

int ireadtktecpp(int attribute, int te)
{

  te = te - 1;

  if(te>TkTeBank->size()) return 0;

  switch (attribute) {
  case 1: 
    return (int)TkTeBank->getSubdetector_ID(te);
    break;
  case 2: 
    return (int)TkTeBank->getSubmodule(te);
     break;
  case 3: 
    return (int)TkTeBank->getUnused(te);
    break;
  case 4: 
    return (int)TkTeBank->getMeasurement_code(te);
    break;
  case 5: 
    return (int)TkTeBank->getPointer_to_end_of_TE(te);
    break;
  case 6: 
    return (int)TkTeBank->getCharge(te);
    break;
  case 7: 
    return (int)TkTeBank->getNdf(te);
    break;
  case 8: 
    return (int)TkTeBank->getChi2(te);
    break;
  case 9: 
    return (int)TkTeBank->getLength(te);
    break;
  case 10: 
    return (int)TkTeBank->getCoord1_of_ref_point(te);
    break;
  case 11: 
    return (int)TkTeBank->getCoord2_of_ref_point(te);
    break;
  case 12: 
    return (int)TkTeBank->getCoord3_of_ref_point(te);
    break;
  case 13: 
    return (int)TkTeBank->getTheta(te);
    break;
  case 14: 
    return (int)TkTeBank->getPhi(te);
    break;
  case 15: 
    return (int)TkTeBank->getInvp(te);
    break;
  case 16: 
    return (int)TkTeBank->getDe_dx(te);
    break;
  case 17: 
    return (int)TkTeBank->getCovmatrix1(te);
    break;
  case 18: 
    return (int)TkTeBank->getCovmatrix2(te);
    break;
  case 19: 
    return (int)TkTeBank->getCovmatrix3(te);
    break;
  case 20: 
    return (int)TkTeBank->getCovmatrix4(te);
    break;
  case 21: 
    return (int)TkTeBank->getCovmatrix5(te);
    break;
  case 22: 
    return (int)TkTeBank->getCovmatrix6(te);
    break;
  case 23: 
    return (int)TkTeBank->getCovmatrix7(te);
    break;
  case 24: 
    return (int)TkTeBank->getCovmatrix8(te);
    break;
  case 25: 
    return (int)TkTeBank->getCovmatrix9(te);
    break;
  case 26: 
    return (int)TkTeBank->getCovmatrix10(te);
    break;
  case 27: 
    return (int)TkTeBank->getCovmatrix11(te);
    break;
  case 28: 
    return (int)TkTeBank->getCovmatrix12(te);
    break;
  case 29: 
    return (int)TkTeBank->getCovmatrix13(te);
    break;
  case 30: 
    return (int)TkTeBank->getCovmatrix14(te);
    break;
  case 31: 
    return (int)TkTeBank->getCovmatrix15(te);
    break;
  default: 
    throw runtime_error("te attribute not valid");
  } 
}


int writetktecpp(float value, int attribute, int te){

  if(te>TkTeBank->size()) return 1;

  te = te - 1;

  switch (attribute) {
  case 1: 
    TkTeBank->setSubdetector_ID((int)value,te);
    return 0;
    break;
  case 2: 
    TkTeBank->setSubmodule((int)value,te);
    return 0;
    break;
  case 3: 
    TkTeBank->setUnused((int)value,te);
    return 0;
    break;
  case 4: 
    TkTeBank->setMeasurement_code((int)value,te);
    return 0;
    break;
  case 5: 
    TkTeBank->setPointer_to_end_of_TE((int)value,te);
    return 0;
    break;
  case 6: 
    TkTeBank->setCharge((int)value,te);
    return 0;
    break;
  case 7: 
    TkTeBank->setNdf((int)value,te);
    return 0;
    break;
  case 8: 
    TkTeBank->setChi2(value,te);
    return 0;
    break;
  case 9: 
    TkTeBank->setLength(value,te);
    return 0;
    break;
  case 10: 
    TkTeBank->setCoord1_of_ref_point(value,te);
    return 0;
    break;
  case 11: 
    TkTeBank->setCoord2_of_ref_point(value,te);
    return 0;
    break;
  case 12: 
    TkTeBank->setCoord3_of_ref_point(value,te);
    return 0;
    break;
  case 13: 
    TkTeBank->setTheta(value,te);
    return 0;
    break;
  case 14: 
    TkTeBank->setPhi(value,te);
    return 0;
    break;
  case 15: 
    TkTeBank->setInvp(value,te);
    return 0;
    break;
  case 16: 
    TkTeBank->setDe_dx(value,te);
    return 0;
    break;
  case 17: 
    TkTeBank->setCovmatrix1(value,te);
    return 0;
    break;
  case 18: 
    TkTeBank->setCovmatrix2(value,te);
    return 0;
    break;
  case 19: 
    TkTeBank->setCovmatrix3(value,te);
    return 0;
    break; 
  case 20: 
    TkTeBank->setCovmatrix4(value,te);
    return 0;
    break;
  case 21: 
    TkTeBank->setCovmatrix5(value,te);
    return 0;
    break;
  case 22: 
    TkTeBank->setCovmatrix6(value,te);
    return 0;
    break;
  case 23: 
    TkTeBank->setCovmatrix7(value,te);
    return 0;
    break;
  case 24: 
    TkTeBank->setCovmatrix8(value,te);
    return 0;
    break;
  case 25: 
    TkTeBank->setCovmatrix9(value,te);
    return 0;
    break;
  case 26: 
    TkTeBank->setCovmatrix10(value,te);
    return 0;
    break;
  case 27: 
    TkTeBank->setCovmatrix11(value,te);
    return 0;
    break;
  case 28: 
    TkTeBank->setCovmatrix12(value,te);
    return 0;
    break;
  case 29: 
    TkTeBank->setCovmatrix13(value,te);
    return 0;
    break;
  case 30: 
    TkTeBank->setCovmatrix14(value,te);
    return 0;
    break;
  case 31: 
    TkTeBank->setCovmatrix15(value,te);
    return 0;
    break;
  default: 
    cout << "attribute = " << attribute << endl;  
    throw runtime_error("te attribute not valid");
  } 
  
}

int addhittktecpp(int hit, int te) 
{
  te = te - 1;
  hit = hit - 1;

  TkTeBank->addHit(hit,te);
  return 0;
}

int writetkitedatcpp(int value, int attribute, int te){
  te = te - 1;
  
  switch (attribute) {
  case 1: 
    TkTeBank->setPosOfFirstHitInHitList(value,te);
    return 0;
    break;
  case 2: 
    TkTeBank->setNumOfHits(value,te);
    return 0;
    break;
  case 3: 
    TkTeBank->setPointrToFirstExclusion(value,te);
    return 0;
    break;
  case 4: 
    TkTeBank->setNumOfExclusions(value,te);
    return 0;
    break;
  case 5: 
    TkTeBank->setTrackNo(value,te);
    return 0;
    break;
  default: 
    std::cout << "attribute = " << attribute <<  std::endl;
    throw runtime_error("writetkteitedatcpp: te attribute not valid");
  } 

}

int readtkitedatcpp(int attribute, int te)
{

  te = te - 1;

  if(te>TkTeBank->size()) return 0;

  switch (attribute) {
  case 1: 
    return TkTeBank->getPosOfFirstHitInHitList(te);
    break;
  case 2: 
    return TkTeBank->getNumOfHits(te);
    break;
  case 3: 
    return TkTeBank->getPointrToFirstExclusion(te);
    break;
  case 4: 
    return TkTeBank->getNumOfExclusions(te);
    break;
  case 5: 
    return TkTeBank->gettrackNo(te);
    break;
  default:
    std::cout << "attribute = " << attribute << std::endl ;
    std::cout << "te = " << te << std::endl ;
    throw runtime_error("te attribute not valid");
  } 

}
