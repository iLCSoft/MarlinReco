#include "ReconstructedParticleImpl_CopyProcessor.h"
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "UTIL/LCRelationNavigator.h"
#include "marlin/Global.h"

ReconstructedParticleImpl_CopyProcessor aReconstructedParticleImpl_CopyProcessor;

ReconstructedParticleImpl_CopyProcessor::ReconstructedParticleImpl_CopyProcessor()
    : Processor("ReconstructedParticleImpl_CopyProcessor") {

  _description = "ReconstructedParticleImpl_CopyProcessor: Copies each entry of a LCIO collection of "
                 "ReconstructedParticleImpl. Select in the optional parameters which members should not be copied.";

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "InputCollection", "Pandora PFOs", _InputColName,
                          std::string("PandoraPFOs"));

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "OutputCollection", "New Pandora PFOs", _OutputColName,
                           std::string("NewPandoraPFOs"));

  registerOutputCollection(LCIO::LCRELATION, "RelationCollection", "Old to New Pandora PFOs Link", _RelationColName,
                           std::string("Old2NewPandoraPFOsLink"));

  registerProcessorParameter("copyType", "Copy type.", _copyType, bool(true));

  registerProcessorParameter("copyMomentum", "Copy momentum.", _copyMomentum, bool(true));

  registerProcessorParameter("copyEnergy", "Copy energy.", _copyEnergy, bool(true));

  registerProcessorParameter("copyCovMatrix", "Copy covariance matrix.", _copyCovMatrix, bool(true));

  registerProcessorParameter("copyMass", "Copy mass.", _copyMass, bool(true));

  registerProcessorParameter("copyCharge", "Copy charge.", _copyCharge, bool(true));

  registerProcessorParameter("copyReferencePoint", "Copy reference point.", _copyReferencePoint, bool(true));

  registerProcessorParameter("copyParticleIDs", "Copy all particle ID entries.", _copyParticleIDs, bool(true));

  registerProcessorParameter("copyParticleIDUsed", "Copy which particle ID is used.", _copyParticleIDUsed, bool(true));

  registerProcessorParameter("copyGoodnessOfPID", "Copy the goodness of the PID.", _copyGoodnessOfPID, bool(true));

  registerProcessorParameter("copyParticles", "Copy all particles.", _copyParticles, bool(true));

  registerProcessorParameter("copyClusters", "Copy all clusters.", _copyClusters, bool(true));

  registerProcessorParameter("copyTracks", "Copy all tracks.", _copyTracks, bool(true));

  registerProcessorParameter("copyStartVertex", "Copy start vertex.", _copyStartVertex, bool(true));
}

void ReconstructedParticleImpl_CopyProcessor::init() {
  // usually a good idea to
  printParameters();

  _nEvt = 0;
}

void ReconstructedParticleImpl_CopyProcessor::processEvent(LCEvent* evt) {
  LCCollection *incol{}, *relcol{};
  int n_incol = 0;

  try {
    incol = evt->getCollection(_InputColName);
    n_incol = incol->getNumberOfElements();
  } catch (DataNotAvailableException& e) {
    streamlog_out(MESSAGE) << "Input collection not found in event " << _nEvt << std::endl;
  }

  LCCollectionVec* outcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCRelationNavigator relnav(LCIO::RECONSTRUCTEDPARTICLE, LCIO::RECONSTRUCTEDPARTICLE);

  for (int i = 0; i < n_incol; ++i) {
    ReconstructedParticleImpl* inpart = dynamic_cast<ReconstructedParticleImpl*>(incol->getElementAt(i));
    ReconstructedParticleImpl* outpart = new ReconstructedParticleImpl;

    if (_copyType)
      outpart->setType(inpart->getType());
    if (_copyMomentum)
      outpart->setMomentum(inpart->getMomentum());
    if (_copyEnergy)
      outpart->setEnergy(inpart->getEnergy());
    if (_copyCovMatrix)
      outpart->setCovMatrix(inpart->getCovMatrix());
    if (_copyMass)
      outpart->setMass(inpart->getMass());
    if (_copyCharge)
      outpart->setCharge(inpart->getCharge());
    if (_copyReferencePoint)
      outpart->setReferencePoint(inpart->getReferencePoint());
    if (_copyParticleIDs)
      for (unsigned int j = 0; j < inpart->getParticleIDs().size(); ++j)
        outpart->addParticleID(inpart->getParticleIDs()[j]);
    if (_copyParticleIDUsed)
      outpart->setParticleIDUsed(inpart->getParticleIDUsed());
    if (_copyGoodnessOfPID)
      outpart->setGoodnessOfPID(inpart->getGoodnessOfPID());
    if (_copyParticles)
      for (unsigned int j = 0; j < inpart->getParticles().size(); ++j)
        outpart->addParticle(inpart->getParticles()[j]);
    if (_copyClusters)
      for (unsigned int j = 0; j < inpart->getClusters().size(); ++j)
        outpart->addCluster(inpart->getClusters()[j]);
    if (_copyTracks)
      for (unsigned int j = 0; j < inpart->getTracks().size(); ++j)
        outpart->addTrack(inpart->getTracks()[j]);
    if (_copyStartVertex)
      outpart->setStartVertex(inpart->getStartVertex());

    outcol->addElement(outpart);
    relnav.addRelation(inpart, outpart, 1);
  }

  relcol = relnav.createLCCollection();
  evt->addCollection(outcol, _OutputColName);
  evt->addCollection(relcol, _RelationColName);

  _nEvt++;
}

void ReconstructedParticleImpl_CopyProcessor::end() {}
