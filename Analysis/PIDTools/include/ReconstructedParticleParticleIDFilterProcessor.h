#ifndef ReconstructedParticleParticleIDFilterProcessor_H
#define ReconstructedParticleParticleIDFilterProcessor_H 1

#include "marlin/Processor.h"

#include <string>
#include <vector>

class ReconstructedParticleParticleIDFilterProcessor : public marlin::Processor {
public:
  Processor* newProcessor() override { return new ReconstructedParticleParticleIDFilterProcessor(); }

  ReconstructedParticleParticleIDFilterProcessor();
  ReconstructedParticleParticleIDFilterProcessor(const ReconstructedParticleParticleIDFilterProcessor&)            = delete;
  ReconstructedParticleParticleIDFilterProcessor& operator=(const ReconstructedParticleParticleIDFilterProcessor&) = delete;
  ~ReconstructedParticleParticleIDFilterProcessor()                                                                = default;

  void processEvent(LCEvent* event) override;

private:
  std::string              m_inputCollName{};
  std::vector<std::string> m_filterPidAlgos{};
};

#endif  // ReconstructedParticleParticleIDFilterProcessor_H
