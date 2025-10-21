#include "ConstantStorage.hh"
#include <EVENT/MCParticle.h>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#ifndef _DecayChain_hh
#define _DecayChain_hh
namespace TTbarAnalysis {
class DecayChain {
public:
  //
  //	Constants
  //

  //
  //	Constructors
  //
  DecayChain(const std::vector<EVENT::MCParticle*>* particles, std::string name, int pdg);
  DecayChain(std::string name, int pdg);
  virtual ~DecayChain() {};
  //
  //	Methods
  //
  void Add(EVENT::MCParticle* particle);
  EVENT::MCParticle* Get(int i);
  int GetSize() const;
  int GetParentPDG() const;
  std::string GetName() const;
  const std::vector<EVENT::MCParticle*>& GetAll() const;
  void Merge(DecayChain& other);
  EVENT::MCParticle* Find(PDGTYPE type) const;
  // oid Print();
private:
  //
  //	Data
  //
  std::vector<EVENT::MCParticle*> myParticles{};
  std::string myName{};
  int myPDG{};
  //
  //	Private methods
  //
};
} // namespace TTbarAnalysis
#endif
