// BcEnergyDensity class by A.Sapronov (sapronov@cern.ch)
// 18/07/2008

#ifndef BcEnergyDensity_h
#define BcEnergyDensity_h 1

#include "TROOT.h"

#include <vector>

using namespace std;

class BcEnergyDensity {
public:
  BcEnergyDensity(const char* inputfilename) { Init(inputfilename); }
  ~BcEnergyDensity() { Destroy(); };

  // returns energy density in cell at given coordinates
  // layers start from 1, backward direction - negative layer number
  // radius is in mm, phi in radians,
  // pEnDens - where the energy density will be written to (in GeV/mm2)
  // pEnDensError - its statistical error
  Bool_t GetEnergyDensity(const Int_t& rLayer, const Double_t& rRadius, const Double_t& rPhi, Double_t* pEnDens,
                          Double_t* pEnDensError) const;

private:
  static const int msSize = 3;
  typedef struct {
    Int_t Id[msSize];
    Double_t R;
    Double_t Phi;
    Double_t EnDens;
    Double_t EnDensErr;
  } CellType;

  typedef struct {
    Double_t Val; // value
    Double_t Lb;  // lower bound
    Double_t Ub;  // upper bound
  } ValErrType;

  vector<CellType*> mCellStorage;

  Int_t mNumberOfRings;
  Double_t mRingDeltaR;
  Double_t mRMin;
  Double_t mRMax;
  Double_t mPhiMin;
  vector<Int_t> mNumbersOfSegments;
  vector<Double_t> mSegmentDeltaPhi;

private:
  void Init(const char* inputfilename);
  void Destroy();
  Int_t GetCellNumber(const Int_t& rLayer, const Double_t& rRadius, const Double_t& rPhi) const;

  Bool_t MakeInterpolation(const Double_t& rX, const vector<Double_t>& rXData, const vector<ValErrType>& rYData,
                           ValErrType* pResult) const;
};

#endif
