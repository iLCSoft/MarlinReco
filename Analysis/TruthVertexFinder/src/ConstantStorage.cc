#include "ConstantStorage.hh"
using std::vector;
namespace TTbarAnalysis {
vector<int> ConstantStorage::myCharmedMesonsPDGs;
vector<int> ConstantStorage::myBottomMesonsPDGs;
vector<int> ConstantStorage::myStrangeMesonsPDGs;
vector<int> ConstantStorage::myBottomBaryonsPDGs;
vector<int> ConstantStorage::myBottomHadronsPDGs;
vector<int> ConstantStorage::myCharmedBaryonsPDGs;
vector<int> ConstantStorage::myCharmedHadronsPDGs;
vector<int> ConstantStorage::myTrackableParticlesPDGs;
vector<int> ConstantStorage::myTauLeptonPDGs;
vector<int> ConstantStorage::myNonTrackableParticlesPDGs;
vector<int> ConstantStorage::myEmptyPDGs;

ConstantStorage::_init ConstantStorage::_initializer;

const std::vector<int>& ConstantStorage::TAU_LEPTON_PDG() { return myTauLeptonPDGs; }
const std::vector<int>& ConstantStorage::CHARMED_MESONS_PDG() { return myCharmedMesonsPDGs; }
const std::vector<int>& ConstantStorage::STRANGE_MESONS_PDG() { return myStrangeMesonsPDGs; }
const std::vector<int>& ConstantStorage::BOTTOM_MESONS_PDG() { return myBottomMesonsPDGs; }
const std::vector<int>& ConstantStorage::BOTTOM_BARYONS_PDG() { return myBottomBaryonsPDGs; }
const std::vector<int>& ConstantStorage::BOTTOM_HADRONS_PDG() { return myBottomHadronsPDGs; }
const std::vector<int>& ConstantStorage::CHARMED_BARYONS_PDG() { return myCharmedBaryonsPDGs; }
const std::vector<int>& ConstantStorage::CHARMED_HADRONS_PDG() { return myCharmedHadronsPDGs; }
const std::vector<int>& ConstantStorage::TRACKABLE_PARTICLES_PDG() { return myTrackableParticlesPDGs; }
const std::vector<int>& ConstantStorage::NONTRACKABLE_PARTICLES_PDG() { return myNonTrackableParticlesPDGs; }
const std::vector<int>& ConstantStorage::GET_PDG(PDGTYPE type) {
  switch (type) {
  case BOTTOM_MESONS:
    return BOTTOM_MESONS_PDG();
  case CHARMED_MESONS:
    return CHARMED_MESONS_PDG();
  case STRANGE_MESONS:
    return STRANGE_MESONS_PDG();
  case BOTTOM_BARYONS:
    return BOTTOM_BARYONS_PDG();
  case BOTTOM_HADRONS:
    return BOTTOM_HADRONS_PDG();
  case CHARMED_HADRONS:
    return CHARMED_HADRONS_PDG();
  case CHARMED_BARYONS:
    return CHARMED_BARYONS_PDG();
  case TAU_LEPTON:
    return TAU_LEPTON_PDG();
  case TRACKABLE_PARTICLES:
    return TRACKABLE_PARTICLES_PDG();
  case NONTRACKABLE_PARTICLES:
    return NONTRACKABLE_PARTICLES_PDG();
  default:
    return ConstantStorage::myEmptyPDGs;
  }
}

} // namespace TTbarAnalysis
