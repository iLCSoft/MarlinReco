#include "DigiHitExtended.h"

DigiHitExtended::DigiHitExtended() { _simHitVec.clear(); }

DigiHitExtended::~DigiHitExtended() {}

void DigiHitExtended::addSimHit(SimCalorimeterHit* simhit) { _simHitVec.push_back(simhit); }

SimCalorimeterHitVec& DigiHitExtended::getSimHitVector() { return _simHitVec; }

void DigiHitExtended::setAmpl(float ampl) { _ampl = ampl; }

void DigiHitExtended::setCellID(int cellid) { _cellid = cellid; }

void DigiHitExtended::setPosition(float* pos) {
  _pos[0] = pos[0];
  _pos[1] = pos[1];
  _pos[2] = pos[2];
}

float DigiHitExtended::getAmpl() { return _ampl; }

int DigiHitExtended::getCellID() { return _cellid; }

float* DigiHitExtended::getPosition() { return _pos; }
