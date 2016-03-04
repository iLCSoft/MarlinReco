/*
 * PIDParticles.cc
 *
 *  Created on: Dec 23, 2015
 *      Author: strahinja
 */

#include "PIDParticles.hh"


PIDParticles::ParticleMap* PIDParticles::CreateParticleMap() {

  ParticleMap *parameterMap = new  ParticleMap;

  parameterMap->insert(std::pair<particleType, PIDParticle_base>
      (electron, electronProperties));

  parameterMap->insert(std::pair<particleType, PIDParticle_base>
      (muon, muonProperties));

  parameterMap->insert(std::pair<particleType, PIDParticle_base>
      (pion, pionProperties));

  parameterMap->insert(std::pair<particleType, PIDParticle_base>
      (kaon, kaonProperties));

  parameterMap->insert(std::pair<particleType, PIDParticle_base>
      (proton, protonProperties));

  return parameterMap;
}


PIDParticles::LLHypothesesMap* PIDParticles::CreateLLPIDMap(std::vector<float> priors) {

  LLHypothesesMap *parameterMap = new  LLHypothesesMap;

  parameterMap->insert(std::pair<particleType, LLPIDHypothesis>
      (electron, LLPIDHypothesis(electronProperties, priors.at(electron))));

  parameterMap->insert(std::pair<particleType, LLPIDHypothesis>
      (muon, LLPIDHypothesis(muonProperties, priors.at(muon))));

  parameterMap->insert(std::pair<particleType, LLPIDHypothesis>
      (pion, LLPIDHypothesis(pionProperties, priors.at(pion))));

  parameterMap->insert(std::pair<particleType, LLPIDHypothesis>
      (kaon, LLPIDHypothesis(kaonProperties, priors.at(kaon))));

  parameterMap->insert(std::pair<particleType, LLPIDHypothesis>
      (proton, LLPIDHypothesis(protonProperties, priors.at(proton))));

  return parameterMap;
}


PIDParticles::MVAHypothesesMap* PIDParticles::CreateMVAPIDMap() {

  MVAHypothesesMap *parameterMap = new  MVAHypothesesMap;

  parameterMap->insert(std::pair<particleType, MVAPIDHypothesis>
      (electron, MVAPIDHypothesis(electronProperties)));

  parameterMap->insert(std::pair<particleType, MVAPIDHypothesis>
      (muon, MVAPIDHypothesis(muonProperties)));

  parameterMap->insert(std::pair<particleType, MVAPIDHypothesis>
      (pion, MVAPIDHypothesis(pionProperties)));

  parameterMap->insert(std::pair<particleType, MVAPIDHypothesis>
      (kaon, MVAPIDHypothesis(kaonProperties)));

  parameterMap->insert(std::pair<particleType, MVAPIDHypothesis>
      (proton, MVAPIDHypothesis(protonProperties)));

  return parameterMap;
}


void PIDParticles::MVAPIDHypothesis::Evaluate(const TString &method) {
  _mva = _reader->EvaluateMVA(method);
  int mvaBin = _histoQ->FindFixBin(_mva);
  if(_mva-_histoQ->GetBinLowEdge(mvaBin) < _histoQ->GetBinWidth(mvaBin)) mvaBin--;
  _q = _histoQ->GetBinContent(mvaBin);

  _sigAbove = _histoSig->Integral(mvaBin, _histoSig->GetNbinsX());
}


/* test */
