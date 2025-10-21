#ifndef _G2CD_hh_
#define _G2CD_hh_

#include <EVENT/CalorimeterHit.h>
#include <TNtuple.h>
#include <TObject.h>
#include <fstream>
#include <iostream>
#include <marlin/Processor.h>
#include <string>

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
class TTree;

// namespace CALICE {

class G2CD : public marlin::Processor {
public:
  G2CD(const G2CD&) = delete;
  G2CD& operator=(const G2CD&) = delete;

  Processor* newProcessor() { return new G2CD; }

  G2CD();

  ~G2CD() {};

  void init();

  void processEvent(LCEvent* evtP);

  void end();

protected:
  std::string _treeFileName{};
  std::string _treeName{};
  std::string _colName{};
  std::vector<std::string> _hcalCollections{};
  std::vector<std::string> _outputHcalCollections{};
  std::vector<std::string> _ecalCollections{};
  std::vector<std::string> _outputEcalCollections{};
  std::vector<float> _ChargeSpatialDistri{};
  std::vector<float> _thresholdHcal{};
  std::vector<float> _calibCoeffEcal{};
  std::vector<int> _ShowerPositionShiftID{}; // should be of the form deltaI, J, K
  float _thresholdEcal{};
  int _NEcalThinLayer{};
  int _overwrite{};
  int _DigiCellSize{};
  int _UsingDefaultDetector{};
  float _PolyaParaA{}, _PolyaParaB{}, _PolyaParaC{};
  float _ChanceOfKink{}, _KinkHitChargeBoost{};
  TTree* _outputTree{};
  TH1F *_NH1stLayer{}, *_NH8thLayer{};
  TF1* _QPolya{};

  int _Num{};
  int _eventNr{};

  int _M{}, _S{}, _I{}, _J{}, _K{}, _Seg{};
  float _PosX{}, _PosY{}, _PosZ{};
  float _EDepo{}, _Charge{};
  int _NHit1mm{}, _NHit1mmCenter{}, _NHit1mmCorner{}, _NHit1mmSide{};
  int _TotalNHit1mm{}, _TotalNHit{}, _TotalMultiHit{};
  int _N1{}, _N2{}, _N3{};

  std::string _fileName{};
  std::ostream* _output{};
};

#endif
