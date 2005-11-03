#include <iostream>
#include <vector>
//stl exception handler
#include<stdexcept>

#include "tktkbank.h"

using namespace std;

// Global pointer to the Tk_Tk_Bank structure
Tk_Tk_Bank * TkTkBank;


Tk_Tk_Bank::Tk_Tk_Bank()
{
}

Tk_Tk_Bank::~Tk_Tk_Bank()
{
}

void Tk_Tk_Bank::add_tk(int modid,int subdetbits,int MesrCode,int tracktype,int numtes,int Charge,int unused,int ndf,float chi2,float L,float xstart,float ystart,float zstart,float xend,float yend, float zend,float cord1,float cord2,float cord3,float theta,float phi,float invp,float cov1,float cov2,float cov3,float cov4,float cov5,float cov6,float cov7,float cov8,float cov9,float cov10,float cov11,float cov12,float cov13,float cov14,float cov15)
{

  tk_tk atk;
  atk.modid = modid;
  atk.subdetbits = subdetbits;
  atk.measurement_code = MesrCode;
  atk.tracktype = tracktype;
  atk.numtes = numtes;
  atk.charge = Charge;
  atk.unused = unused;
  atk.ndf = ndf;
  atk.chi2 = chi2;
  atk.length = L;
  atk.xstart = xstart ;
  atk.ystart = ystart ;
  atk.zstart = zstart ;
  atk.xend = xend ;
  atk.yend = yend ;
  atk.zend = zend ;
  atk.coord1_of_ref_point = cord1;
  atk.coord2_of_ref_point = cord2;
  atk.coord3_of_ref_point = cord3;
  atk.theta = theta;
  atk.phi = phi;
  atk.invp = invp;
  atk.covmatrix1 = cov1;
  atk.covmatrix2 = cov2;
  atk.covmatrix3 = cov3;
  atk.covmatrix4 = cov4;
  atk.covmatrix5 = cov5;
  atk.covmatrix6 = cov6;
  atk.covmatrix7 = cov7;
  atk.covmatrix8 = cov8;
  atk.covmatrix9 = cov9;
  atk.covmatrix10 = cov10;
  atk.covmatrix11 = cov11;
  atk.covmatrix12 = cov12;
  atk.covmatrix13 = cov13;
  atk.covmatrix14 = cov14;
  atk.covmatrix15 = cov15;


  tk_bank.push_back(atk);

  unsigned int size = itkdat_bank.size() ;
  itkdat_bank.resize(size+1) ;

}

// due to the fact that the bank will be accessed by integer value of the tk number this is probaly dangerous
// so leave it out for now

//void Tk_Tk_Bank::remove_tk(int tk)

//{
//  tk_bank.erase(tk_bank.begin()+tk);
//}

//int tkmktkcpp(int subid,int submod,int unused,int MesrCode,int PntkTK,int Q,int ndf,float chi2,float L,float cord1,float cord2,float cord3,float theta,float phi,float invp,float dedx,float cov1,float cov2,float cov3,float cov4,float cov5,float cov6,float cov7,float cov8,float cov9,float cov10,float cov11,float cov12,float cov13,float cov14,float cov15)
//{
// 
//  TkTkBank->add_tk(subid, submod, unused, MesrCode, PnteTK, Q, ndf, chi2, L, cord1, cord2, cord3, theta, phi, invp, dedx, cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8, cov9, cov10, cov11, cov12, cov13, cov14, cov15);
// 
//  return 0;
//}

int tkmktkcpp(int modid,int subdetbits,int MesrCode,int tracktype, int numtes,int Charge,int unused,int ndf,float chi2,float L,float xstart, float ystart, float zstart, float xend, float yend, float zend, float cord1,float cord2,float cord3,float theta,float phi,float invp,float* cov)
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

  TkTkBank->add_tk(modid, subdetbits, MesrCode, tracktype, numtes, Charge, unused, ndf, chi2, L, xstart, ystart, zstart, xend, yend, zend, cord1, cord2, cord3, theta, phi, invp, cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8, cov9, cov10, cov11, cov12, cov13, cov14, cov15);
  return 0;
}

float rreadtktkcpp(int attribute, int tk)
{

  tk = tk - 1;

  if(tk>TkTkBank->size()) return 0.;

  switch (attribute) {
  case 1: 
    return TkTkBank->getMod_ID(tk);
    break;
  case 2: 
    return TkTkBank->getSubdetbits(tk);
     break;
  case 3: 
    return TkTkBank->getMeasurement_code(tk);
    break;
  case 4: 
    return TkTkBank->getTrackType(tk);
    break;
  case 5: 
    return TkTkBank->getNumTEs(tk);
    break;
  case 6: 
    return TkTkBank->getCharge(tk);
    break;
  case 7: 
    return TkTkBank->getUnused(tk);
    break;
  case 8: 
    return TkTkBank->getNdf(tk);
    break;
  case 9: 
    return TkTkBank->getChi2(tk);
    break;
  case 10: 
    return TkTkBank->getLength(tk);
    break;
  case 11:
    return TkTkBank->getStartX(tk); 
    break;
  case 12:
    return TkTkBank->getStartY(tk); 
    break;
  case 13:
    return TkTkBank->getStartZ(tk); 
    break;
  case 14:
    return TkTkBank->getEndX(tk); 
    break;
  case 15:
    return TkTkBank->getEndY(tk); 
    break;
  case 16:
    return TkTkBank->getEndZ(tk); 
    break;
  case 17: 
    return TkTkBank->getCoord1_of_ref_point(tk);
    break;
  case 18: 
    return TkTkBank->getCoord2_of_ref_point(tk);
    break;
  case 19: 
    return TkTkBank->getCoord3_of_ref_point(tk);
    break;
  case 20: 
    return TkTkBank->getTheta(tk);
    break;
  case 21: 
    return TkTkBank->getPhi(tk);
    break;
  case 22: 
    return TkTkBank->getInvp(tk);
    break;
  case 23: 
    return TkTkBank->getCovmatrix1(tk);
    break;
  case 24: 
    return TkTkBank->getCovmatrix2(tk);
    break;
  case 25: 
    return TkTkBank->getCovmatrix3(tk);
    break;
  case 26: 
    return TkTkBank->getCovmatrix4(tk);
    break;
  case 27: 
    return TkTkBank->getCovmatrix5(tk);
    break;
  case 28: 
    return TkTkBank->getCovmatrix6(tk);
    break;
  case 29: 
    return TkTkBank->getCovmatrix7(tk);
    break;
  case 30: 
    return TkTkBank->getCovmatrix8(tk);
    break;
  case 31: 
    return TkTkBank->getCovmatrix9(tk);
    break;
  case 32: 
    return TkTkBank->getCovmatrix10(tk);
    break;
  case 33: 
    return TkTkBank->getCovmatrix11(tk);
    break;
  case 34: 
    return TkTkBank->getCovmatrix12(tk);
    break;
  case 35: 
    return TkTkBank->getCovmatrix13(tk);
    break;
  case 36: 
    return TkTkBank->getCovmatrix14(tk);
    break;
  case 37: 
    return TkTkBank->getCovmatrix15(tk);
    break;
  default: 
    throw runtime_error("tk attribute not valid");
  } 
}

int ireadtktkcpp(int attribute, int tk)
{

  tk = tk - 1;

  if(tk>TkTkBank->size()) return 0;

  switch (attribute) {
  case 1: 
    return (int)TkTkBank->getMod_ID(tk);
    break;
  case 2: 
    return (int)TkTkBank->getSubdetbits(tk);
     break;
  case 3: 
    return (int)TkTkBank->getMeasurement_code(tk);
    break;
  case 4: 
    return (int)TkTkBank->getTrackType(tk);
    break;
  case 5: 
    return (int)TkTkBank->getNumTEs(tk);
    break;
  case 6:
    return (int)TkTkBank->getCharge(tk);
    break;
  case 7: 
    return (int)TkTkBank->getUnused(tk);
    break;
  case 8: 
    return (int)TkTkBank->getNdf(tk);
    break;
  case 9: 
    return (int)TkTkBank->getChi2(tk);
    break;
  case 10:
    return (int)TkTkBank->getLength(tk);
    break;
  case 11: 
    return (int)TkTkBank->getStartX(tk);
    break;
  case 12: 
    return (int)TkTkBank->getStartY(tk);
    break;
  case 13: 
    return (int)TkTkBank->getStartZ(tk);
    break;
  case 14: 
    return (int)TkTkBank->getEndX(tk);
    break;
  case 15: 
    return (int)TkTkBank->getEndY(tk);
    break;
  case 16: 
    return (int)TkTkBank->getEndZ(tk);
    break;
  case 17: 
    return (int)TkTkBank->getCoord1_of_ref_point(tk);
    break;
  case 18: 
    return (int)TkTkBank->getCoord2_of_ref_point(tk);
    break;
  case 19: 
    return (int)TkTkBank->getCoord3_of_ref_point(tk);
    break;
  case 20: 
    return (int)TkTkBank->getTheta(tk);
    break;
  case 21: 
    return (int)TkTkBank->getPhi(tk);
    break;
  case 22: 
    return (int)TkTkBank->getInvp(tk);
    break;
  case 23: 
    return (int)TkTkBank->getCovmatrix1(tk);
    break;
  case 24: 
    return (int)TkTkBank->getCovmatrix2(tk);
    break;
  case 25: 
    return (int)TkTkBank->getCovmatrix3(tk);
    break;
  case 26: 
    return (int)TkTkBank->getCovmatrix4(tk);
    break;
  case 27: 
    return (int)TkTkBank->getCovmatrix5(tk);
    break;
  case 28: 
    return (int)TkTkBank->getCovmatrix6(tk);
    break;
  case 29: 
    return (int)TkTkBank->getCovmatrix7(tk);
    break;
  case 30: 
    return (int)TkTkBank->getCovmatrix8(tk);
    break;
  case 31: 
    return (int)TkTkBank->getCovmatrix9(tk);
    break;
  case 32: 
    return (int)TkTkBank->getCovmatrix10(tk);
    break;
  case 33: 
    return (int)TkTkBank->getCovmatrix11(tk);
    break;
  case 34: 
    return (int)TkTkBank->getCovmatrix12(tk);
    break;
  case 35: 
    return (int)TkTkBank->getCovmatrix13(tk);
    break;
  case 36: 
    return (int)TkTkBank->getCovmatrix14(tk);
    break;
  case 37: 
    return (int)TkTkBank->getCovmatrix15(tk);
    break;
  default: 
    throw runtime_error("tk attribute not valid");
  } 
}


int writetktkcpp(float value, int attribute, int tk){

  if(tk>TkTkBank->size()) return 1;

  tk = tk - 1;

  switch (attribute) {
  case 1: 
    TkTkBank->setMod_ID((int)value,tk);
    return 0;
    break;
  case 2: 
    TkTkBank->setSubdetbits((int)value,tk);
    return 0;
    break;
  case 3: 
    TkTkBank->setMeasurement_code((int)value,tk);
    return 0;
    break;
  case 4: 
    TkTkBank->setTrackType((int)value,tk);
    return 0;
    break;
  case 5: 
    TkTkBank->setNumTEs((int)value,tk);
    return 0;
    break;
  case 6: 
    TkTkBank->setCharge((int)value,tk);
    return 0;
    break;
  case 7: 
    TkTkBank->setUnused((int)value,tk);
    return 0;
    break;
  case 8: 
    TkTkBank->setNdf((int)value,tk);
    return 0;
    break;
  case 9: 
    TkTkBank->setChi2(value,tk);
    return 0;
    break;
  case 10: 
    TkTkBank->setLength(value,tk);
    return 0;
    break;
  case 11:
    TkTkBank->setStartX(value,tk);
    return 0;
    break;
  case 12:
    TkTkBank->setStartY(value,tk);
    return 0;
    break;
  case 13:
    TkTkBank->setStartZ(value,tk);
    return 0;
    break;
  case 14:
    TkTkBank->setEndX(value,tk);
    return 0;
    break;
  case 15:
    TkTkBank->setEndY(value,tk);
    return 0;
    break;
  case 16:
    TkTkBank->setEndZ(value,tk);
    return 0;
    break;
  case 17: 
    TkTkBank->setCoord1_of_ref_point(value,tk);
    return 0;
    break;
  case 18: 
    TkTkBank->setCoord2_of_ref_point(value,tk);
    return 0;
    break;
  case 19: 
    TkTkBank->setCoord3_of_ref_point(value,tk);
    return 0;
    break;
  case 20: 
    TkTkBank->setTheta(value,tk);
    return 0;
    break;
  case 21: 
    TkTkBank->setPhi(value,tk);
    return 0;
    break;
  case 22: 
    TkTkBank->setInvp(value,tk);
    return 0;
    break;
  case 23: 
    TkTkBank->setCovmatrix1(value,tk);
    return 0;
    break;
  case 24: 
    TkTkBank->setCovmatrix2(value,tk);
    return 0;
    break;
  case 25: 
    TkTkBank->setCovmatrix3(value,tk);
    return 0;
    break; 
  case 26: 
    TkTkBank->setCovmatrix4(value,tk);
    return 0;
    break;
  case 27: 
    TkTkBank->setCovmatrix5(value,tk);
    return 0;
    break;
  case 28: 
    TkTkBank->setCovmatrix6(value,tk);
    return 0;
    break;
  case 29: 
    TkTkBank->setCovmatrix7(value,tk);
    return 0;
    break;
  case 30: 
    TkTkBank->setCovmatrix8(value,tk);
    return 0;
    break;
  case 31: 
    TkTkBank->setCovmatrix9(value,tk);
    return 0;
    break;
  case 32: 
    TkTkBank->setCovmatrix10(value,tk);
    return 0;
    break;
  case 33: 
    TkTkBank->setCovmatrix11(value,tk);
    return 0;
    break;
  case 34: 
    TkTkBank->setCovmatrix12(value,tk);
    return 0;
    break;
  case 35: 
    TkTkBank->setCovmatrix13(value,tk);
    return 0;
    break;
  case 36: 
    TkTkBank->setCovmatrix14(value,tk);
    return 0;
    break;
  case 37: 
    TkTkBank->setCovmatrix15(value,tk);
    return 0;
    break;
  default: 
    cout << "attribute = " << attribute << endl;  
    throw runtime_error("tk attribute not valid");
  } 
  
}

int addtetktkcpp(int te , int tk) 
{
  tk = tk - 1;
  te = te - 1;

  TkTkBank->addTE(te ,tk);
  return 0;
}

int writetkitkdatcpp(int value, int attribute, int tk){
  tk = tk - 1;
  
  switch (attribute) {
  case 1: 
    TkTkBank->setPosOfFirstTEInTEList(value,tk);
    return 0;
    break;
  case 2: 
    TkTkBank->setNumOfTEs(value,tk);
    return 0;
    break;
  case 3: 
    TkTkBank->setTrackNo(value,tk);
    return 0;
    break;
  default: 
    std::cout << "attribute = " << attribute <<  std::endl;
    throw runtime_error("writetktkitkdatcpp: tk attribute not valid");
  } 

}

int readtkitkdatcpp(int attribute, int tk)
{

  tk = tk - 1;

  if(tk>TkTkBank->size()) return 0;

  switch (attribute) {
  case 1: 
    return TkTkBank->getPosOfFirstTEInTEList(tk);
    break;
  case 2: 
    return TkTkBank->getNumOfTEs(tk);
    break;
  case 3: 
    return TkTkBank->getTrackNo(tk);
    break;
  default:
    std::cout << "attribute = " << attribute << std::endl ;
    std::cout << "tk = " << tk << std::endl ;
    throw runtime_error("tk attribute not valid");
  } 

}
