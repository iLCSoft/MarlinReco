#include "VertexRestorer.hh"

using namespace lcio;
using namespace marlin;
using EVENT::ParticleID;
using EVENT::Track;
using IMPL::ParticleIDImpl;
using IMPL::ReconstructedParticleImpl;
using IMPL::VertexImpl;
using std::abs;
using std::string;
using std::vector;
namespace TTbarAnalysis {
VertexRestorer aVertexRestorer;

VertexRestorer::VertexRestorer() : Processor("VertexRestorer") {

  // modify processor description
  _description = "VertexRestorer";

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter("ROOTFileName", "Output ROOT File Name", _hfilename, string("VertexRestorer.root"));

  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "PFOCollection", "Name of the Calorimeter hit collection",
                          _colName, string("PandoraPFOs"));

  registerInputCollection(LCIO::VERTEX, "PrimaryCollectionName", "Name of the PrimaryVertex collection", _colPriName,
                          std::string("PrimaryVertex"));
  registerOutputCollection(LCIO::VERTEX, "OutputCollectionName", "Name of the Vertex collection", _outputcolName,
                           std::string("RecoveredJets_vtx"));
  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "OutputJetCollectionName", "Name of the Jet collection",
                           _outputjetcolName, std::string("RecoveredJets"));
  registerOutputCollection(LCIO::LCRELATION, "OutputRelCollectionName", "Name of the jet vertex relation collection",
                           _outputjetrelcolName, std::string("RecoveredJets_rel"));
  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "OutputRelRPCollectionName",
                           "Name of the jet vertex RP collection", _outputjetrelRPcolName,
                           std::string("RecoveredJets_vtx_RP"));
  registerInputCollection(LCIO::VERTEX, "SecondaryCollectionName", "Name of the BuildUpVertex collection", _colSecName,
                          std::string("FinalJets_vtx"));
  registerInputCollection(LCIO::VERTEX, "SecondaryRPCollectionName", "Name of the BuildUpVertex collection",
                          _colSecRPName, std::string("FinalJets_vtx_RP"));
  registerInputCollection(LCIO::VERTEX, "V0CollectionName", "Name of the BuildUpVertex_V0 collection", _colV0Name,
                          std::string("BuildUpVertex_V0"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "JetCollectionName", "Name of the Jet collection", _colJetName,
                          std::string("FinalJets"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "MissedCollectionName", "Name of the Jet collection",
                          _colMissedName, std::string("MissedParticles"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "MCMissedCollectionName", "Name of the Jet collection",
                          _colMCMissedName, std::string("MCMissedParticles"));
  registerInputCollection(LCIO::TRACK, "NotUsedTracksCollectionName", "Name of tracks not passed to PFA collection",
                          _colTrashCanName, std::string("TracksFailBothCanFormPfoFlags"));
  registerInputCollection(LCIO::LCRELATION, "RecoMcTruthCollectionName", "Name of the RecoMcTruth collection",
                          _colRelName, std::string("RecoMCTruthLink"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "MissedVtxCollectionName", "Name of the Jet collection",
                          _colMissedVtxName, std::string("MissedParticlesVtx"));
  registerInputCollection(LCIO::MCPARTICLE, "BStarCollectionName", "Name of the BStar collection", _colBStarName,
                          std::string("BStar"));
  registerInputCollection(LCIO::MCPARTICLE, "EGProngsCollectionName", "Name of the EGProngs collection", _colEGPName,

                          // std::string("MissedParticles")
                          std::string("EGProngs"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "JetRelCollectionName", "Name of the Jet relation collection",
                          _colJetRelName, std::string("FinalJets_rel"));
  _angleCut = 0.05;
  registerProcessorParameter("angleCut", "angle cut for association", _angleCut, _angleCut);
  _offsetCut = 0.05;
  registerProcessorParameter("offsetCut", "Offset cut for association", _offsetCut, _offsetCut);
  _aParameter = 0.005;
  registerProcessorParameter("a", "a parameter of accuracy in mm", _aParameter, _aParameter);
  _bParameter = 0.01;
  registerProcessorParameter("b", "b parameter of accuracy in mm", _bParameter, _bParameter);

  ip[0] = 0.0;
  ip[1] = 0.0;
  ip[2] = 0.0;
}

void VertexRestorer::init() {
  // usually a good idea to
  printParameters();

  _nRun = 0;
  _nEvt = 0;

  hfile = new TFile(_hfilename.c_str(), "RECREATE", _hfilename.c_str());

  _Tree = new TTree("Stats", "tree");

  _Tree->Branch("detectedTotal", &_detectedTotal, "detectedTotal/I");
  _Tree->Branch("missedTotal", &_missedTotal, "missedTotal/I");
  _Tree->Branch("allmissedMomentum", _allmissedMomentum, "allmissedMomentum[missedTotal]/F");
  _Tree->Branch("bstarTotal", &_bstarTotal, "bstarTotal/I");

  _Tree->Branch("missedDetected", &_missedDetected, "missedDetected/I");
  _Tree->Branch("missedMomentum", _missedMomentum, "missedMomentum[missedDetected]/F");
  _Tree->Branch("missedError", _missedError, "missedError[missedDetected]/F");
  _Tree->Branch("missedOffset", _missedOffset, "missedOffset[missedDetected]/F");
  _Tree->Branch("missedAngle", _missedAngle, "missedAngle[missedDetected]/F");
  _Tree->Branch("missedObs", _missedObs, "missedObs[missedDetected]/F");
  _Tree->Branch("missedVtxHits", _missedVtxHits, "missedVtxHits[missedDetected]/I");
  _Tree->Branch("missedFtdHits", _missedFtdHits, "missedFtdHits[missedDetected]/I");
  _Tree->Branch("missedDeviation", _missedDeviation, "missedDeviation[missedDetected]/F");
  _Tree->Branch("missedSecDeviation", _missedSecDeviation, "missedSecDeviation[missedDetected]/F");
  _Tree->Branch("missedCostheta", _missedCostheta, "missedCostheta[missedDetected]/F");

  _Tree->Branch("bstarDetected", &_bstarDetected, "bstarDetected/I");
  _Tree->Branch("bstarMomentum", _bstarMomentum, "bstarMomentum[bstarDetected]/F");
  _Tree->Branch("bstarError", _bstarError, "bstarError[bstarDetected]/F");
  _Tree->Branch("bstarAngle", _bstarAngle, "bstarAngle[bstarDetected]/F");
  _Tree->Branch("bstarDeviation", _bstarDeviation, "bstarDeviation[bstarDetected]/F");

  /*_Tree->Branch("v0Detected", &_v0Detected, "v0Detected/I");
  _Tree->Branch("v0Momentum", _v0Momentum, "v0Momentum[v0Detected]/F");
  _Tree->Branch("v0Error", _v0Error, "v0Error[v0Detected]/F");
  _Tree->Branch("v0Angle", _v0Angle, "v0Angle[v0Detected]/F");
  _Tree->Branch("v0Deviation", _v0Deviation, "v0Deviation[v0Detected]/F");*/

  _Tree->Branch("fakeDetected", &_fakeDetected, "fakeDetected/I");
  _Tree->Branch("fakeMomentum", _fakeMomentum, "fakeMomentum[fakeDetected]/F");
  _Tree->Branch("fakeError", _fakeError, "fakeError[fakeDetected]/F");
  _Tree->Branch("fakeAngle", _fakeAngle, "fakeAngle[fakeDetected]/F");
  _Tree->Branch("fakeObs", _fakeObs, "fakeObs[fakeDetected]/F");
  _Tree->Branch("fakeVtxHits", _fakeVtxHits, "fakeVtxHits[fakeDetected]/I");
  _Tree->Branch("fakeFtdHits", _fakeFtdHits, "fakeFtdHits[fakeDetected]/I");
  _Tree->Branch("fakeDeviation", _fakeDeviation, "fakeDeviation[fakeDetected]/F");
  _Tree->Branch("fakeSecDeviation", _fakeSecDeviation, "fakeSecDeviation[fakeDetected]/F");
  _Tree->Branch("fakeCostheta", _fakeCostheta, "fakeCostheta[fakeDetected]/F");

  _Tree->Branch("purDetected", &_purDetected, "purDetected/I");
  _Tree->Branch("purOffset", _purOffset, "purOffset[purDetected]/F");
  _Tree->Branch("purAngle", _purAngle, "purAngle[purDetected]/F");
  _Tree->Branch("purCostheta", _purCostheta, "purCostheta[purDetected]/F");
  _Tree->Branch("purMomentum", _purMomentum, "purMomentum[purDetected]/F");
  _Tree->Branch("missedRecoTotal", &_missedRecoTotal, "missedRecoTotal/I");
  _Tree->Branch("missedRecoMomentum", _missedRecoMomentum, "missedRecoMomentum[missedRecoTotal]/F");
  _Tree->Branch("missedRecoTheta", _missedRecoTheta, "missedRecoTheta[missedRecoTotal]/F");

  /*
             _primaryTree =  new TTree( "Primaries", "tree" );
             _primaryTree->Branch("primeZ", &_primeZ, "primeZ/F");
             _primaryTree->Branch("primeT", &_primeT, "primeT/F");
             _primaryTree->Branch("primariesTotal", &_primariesTotal, "primariesTotal/I");
             _primaryTree->Branch("primeType", &_primeType, "primeType[primariesTotal]/I");
             _primaryTree->Branch("primeOffset", &_primeOffset, "primeOffset[primariesTotal]/F");
             _primaryTree->Branch("primeMomentum", &_primeMomentum, "primeMomentum[primariesTotal]/F");
             _primaryTree->Branch("primeDeviation", &_primeDeviation, "primeDeviation[primariesTotal]/F");
             _primaryTree->Branch("primeError", &_primeError, "primeError[primariesTotal]/F");
             _primaryTree->Branch("primeTeorError", &_primeTeorError, "primeTeorError[primariesTotal]/F");
             _primaryTree->Branch("primeChi2", &_primeChi2, "primeChi2[primariesTotal]/F");
             _primaryTree->Branch("primeCostheta", &_primeCostheta, "primeCostheta[primariesTotal]/F");
             _primaryTree->Branch("primeVtxHits", &_primeVtxHits, "primeVtxHits[primariesTotal]/I");
             _primaryTree->Branch("primeFtdHits", &_primeFtdHits, "primeFtdHits[primariesTotal]/I");
             _primaryTree->Branch("primeEtdHits", &_primeEtdHits, "primeEtdHits[primariesTotal]/I");
  */
  _purgatoryTree = new TTree("Purgatory", "tree");
  _purgatoryTree->Branch("purgatoryTotal", &_purgatoryTotal, "purgatoryTotal/I");
  _purgatoryTree->Branch("purgatoryOffset", _purgatoryOffset, "purgatoryOffset[purgatoryTotal]/F");
  _purgatoryTree->Branch("purgatoryMomentum", _purgatoryMomentum, "purgatoryMomentum[purgatoryTotal]/F");
  _purgatoryTree->Branch("purgatoryDeviation", _purgatoryDeviation, "purgatoryDeviation[purgatoryTotal]/F");
  _purgatoryTree->Branch("purgatoryError", _purgatoryError, "purgatoryError[purgatoryTotal]/F");
  _purgatoryTree->Branch("purgatoryTeorError", _purgatoryTeorError, "purgatoryTeorError[purgatoryTotal]/F");
  _purgatoryTree->Branch("purgatoryCostheta", _purgatoryCostheta, "purgatoryCostheta[purgatoryTotal]/F");
  _purgatoryTree->Branch("purgatoryVtxHits", _purgatoryVtxHits, "purgatoryVtxHits[purgatoryTotal]/I");
  _purgatoryTree->Branch("purgatoryParent", _purgatoryParent, "purgatoryParent[purgatoryTotal]/I");

  _trashTree = new TTree("Zombies", "tree");
  _trashTree->Branch("trashTotal", &_trashTotal, "trashTotal/I");
  _trashTree->Branch("trashVtxHits", _trashVtxHits, "trashVtxHits[trashTotal]/I");
  _trashTree->Branch("trashTpcHits", _trashTpcHits, "trashTpcHits[trashTotal]/I");
  _trashTree->Branch("trashHasMC", _trashHasMCParticle, "trashHasMC[trashTotal]/I");
  _trashTree->Branch("trashIsReco", _trashIsReco, "trashIsReco[trashTotal]/I");
  _trashTree->Branch("trashIsDublicate", _trashIsDublicate, "trashIsDublicate[trashTotal]/I");
  _trashTree->Branch("trashCostheta", _trashCostheta, "trashCostheta[trashTotal]/F");

  _TestTree = new TTree("Test", "tree");
  /*_TestTree->Branch("ntotaltest", &_ntotaltest, "ntotaltest/I");
  _TestTree->Branch("ngoodtest", &_ngoodtest, "ngoodtest/I");
  _TestTree->Branch("nlosttest", &_nlosttest, "nlosttest/I");
  _TestTree->Branch("nmultitest", &_nmultitest, "nmultitest/I");
  _TestTree->Branch("costhetaTestedTotal", _costhetaTestedTotal, "costhetaTestedTotal[ntotaltest]/F");
  _TestTree->Branch("costhetaTestedGood", _costhetaTestedGood, "costhetaTestedGood[ngoodtest]/F");
  _TestTree->Branch("costhetaTestedLost", _costhetaTestedLost, "costhetaTestedLost[nlosttest]/F");
  _TestTree->Branch("costhetaTestedMultiple", _costhetaTestedMultiple, "costhetaTestedMultiple[nmultitest]/F");*/
  _TestTree->Branch("ntest1", &_ntest, "ntest/I");
  _TestTree->Branch("parameter1", _testparameter, "parameter1[ntest]/F");
  _TestTree->Branch("test1", _testDistance, "test1[ntest]/F");
}

void VertexRestorer::processRunHeader(LCRunHeader* run) { _nRun++; }
float VertexRestorer::GetDeviation(const EVENT::ReconstructedParticle* particle, double* pos) {
  vector<float> direction = MathOperator::getDirection(particle->getMomentum());
  float accuracy = GetError(
      particle); // sqrt(_aParameter*_aParameter + _bParameter*_bParameter /( p * p * pow(sin(angles[1]), 4.0/3.0)) );
  float offset = MathOperator::getDistanceTo(ip, direction, pos);
  return offset / accuracy;
}
void VertexRestorer::PrintTrack(Track* track) {
  if (!track) {
    return;
  }
  std::cout << "|\t" << track->getD0() << "|\t" << track->getZ0() << "|\t" << track->getPhi() << "|\t"
            << track->getOmega() << "|\t" << track->getTanLambda() << '\n';
}

bool VertexRestorer::IsParticleFromIP(const EVENT::ReconstructedParticle* particle) {
  if (!particle) {
    return false;
  }
  const Vertex* vertex = particle->getStartVertex();
  if (vertex) {
    return vertex->isPrimary();
  }
  std::cout << "Particle has no vertex!\n";
  return false;
}
void VertexRestorer::CopyParameters(EVENT::LCParameters& from, EVENT::LCParameters& to) {
  vector<string> intkeys;
  from.getIntKeys(intkeys);
  vector<string> floatkeys;
  from.getFloatKeys(floatkeys);
  vector<string> stringkeys;
  from.getStringKeys(stringkeys);
  // std::cout << "Int keys size: " << intkeys.size() << "\n";
  for (unsigned int i = 0; i < intkeys.size(); i++) {
    vector<int> intvals;
    from.getIntVals(intkeys[i], intvals);
    // std::cout << "Int key "<< intkeys[i] << " has " << intvals.size() << " values\n";
    to.setValues(intkeys[i], intvals);
  }
  // std::cout << "Float keys size: " << floatkeys.size() << "\n";
  for (unsigned int i = 0; i < floatkeys.size(); i++) {
    vector<float> floatvals;
    from.getFloatVals(floatkeys[i], floatvals);
    // std::cout << "Float key "<< floatkeys[i] << " has " << floatvals.size() << " values\n";
    to.setValues(floatkeys[i], floatvals);
  }
  // std::cout << "String keys size: " << stringkeys.size() << "\n";
  for (unsigned int i = 0; i < stringkeys.size(); i++) {
    vector<string> stringvals;
    from.getStringVals(stringkeys[i], stringvals);
    // std::cout << "String key "<< stringkeys[i] << " has " << stringvals.size() << " values\n";
    to.setValues(stringkeys[i], stringvals);
  }
}
void VertexRestorer::PrintParticle(ReconstructedParticle* particle) {
  if (!particle) {
    return;
  }
  int vtxhits = 0;
  int tpchits = 0;
  std::cout << std::fixed << std::setw(6) << std::setprecision(3) << std::setfill(' ');
  // int id =  particle->getTracks()[0]->getSubdetectorHitNumbers()[4];
  if (abs(particle->getCharge()) > 0.5) {
    vtxhits = particle->getTracks()[0]->getSubdetectorHitNumbers()[0];
    tpchits = particle->getTracks()[0]->getSubdetectorHitNumbers()[6];
  }
  // float chi2 = (particle->getTracks().size() > 0)? particle->getTracks()[0]->getChi2() /
  // (float)particle->getTracks()[0]->getNdf() : -1.0;
  std::cout << "|" << vtxhits << ":" << tpchits << "\t\t|" << particle->getMass() << "\t\t|" << particle->getCharge()
            << "\t\t|" << particle->getEnergy() << "\t\t|\n";
}
void VertexRestorer::processEvent(LCEvent* evt) {
  _nEvt++;

  try {
    std::cout << "********VertexRecovery*" << _nEvt << "**********\n";
    // TrackOperator opera;
    // opera.test();
    LCCollection* pricol = evt->getCollection(_colPriName);
    Vertex* primary = dynamic_cast<Vertex*>(pricol->getElementAt(0));
    myPrimary = primary;
    // LCCollection* v0col = evt->getCollection( _colV0Name );
    LCCollection* maincol = evt->getCollection(_colName);
    LCCollection* seccol = evt->getCollection(_colSecName);
    // LCCollection* relation = evt->getCollection(_colRelName);
    LCCollection* secRPcol = evt->getCollection(_colSecRPName);

    RecoveryOperator recovery(primary, maincol);
    /*LCCollection* trashcol = NULL; //evt->getCollection( _colTrashCanName );
    vector< Vertex * > newvertices = recovery.RecoverBuildVertices(seccol,trashcol);
    IMPL::LCCollectionVec * newseccol = new IMPL::LCCollectionVec ( EVENT::LCIO::VERTEX ) ;
    IMPL::LCCollectionVec * newvtxRPcol = new IMPL::LCCollectionVec ( EVENT::LCIO::RECONSTRUCTEDPARTICLE ) ;
    for (int i = 0; i < newvertices.size(); i++)
    {
            newseccol->addElement(newvertices[i]);
            newvtxRPcol->addElement(newvertices[i]->getAssociatedParticle());
    }
    CopyParameters(seccol->parameters (), newseccol->parameters ());
    evt->addCollection( newseccol, _outputcolName ) ;
    CopyParameters(secRPcol->parameters (), newvtxRPcol->parameters ());
    evt->addCollection( newvtxRPcol, _outputjetrelRPcolName ) ;
    */
    LCCollection* jetcol = evt->getCollection(_colJetName);
    LCCollection* relcol = evt->getCollection(_colJetRelName);
    // vector < ReconstructedParticle * > * result = RestoreVerticesPFO(maincol, seccol);

    LCCollection* trashcol = NULL; // evt->getCollection( _colTrashCanName );
    // IMPL::LCCollectionVec * newjetcol = new IMPL::LCCollectionVec ( (const IMPL::LCCollectionVec&)(*jetcol) ) ;
    IMPL::LCCollectionVec* newjetcol = new IMPL::LCCollectionVec(EVENT::LCIO::RECONSTRUCTEDPARTICLE);
    newjetcol->setSubset();
    IMPL::LCCollectionVec* newjetrelcol = new IMPL::LCCollectionVec(EVENT::LCIO::LCRELATION);
    IMPL::LCCollectionVec* newseccol = new IMPL::LCCollectionVec(EVENT::LCIO::VERTEX);
    IMPL::LCCollectionVec* newvtxRPcol = new IMPL::LCCollectionVec(EVENT::LCIO::RECONSTRUCTEDPARTICLE);
    vector<Vertex*> newvertices = recovery.RecoverJetVertices(jetcol, relcol, seccol, trashcol, newjetrelcol);
    for (unsigned int i = 0; i < newvertices.size(); i++) {
      newseccol->addElement(newvertices[i]);
      newvtxRPcol->addElement(newvertices[i]->getAssociatedParticle());
    }
    for (int i = 0; i < jetcol->getNumberOfElements(); i++) {
      newjetcol->addElement(jetcol->getElementAt(i));
    }
    // EVENT::LCParameters * copiedparameters = new LCParametersImpl((const IMPL::LCParametersImpl&)jetcol->parameters
    // ()); EVENT::LCParameters * newparameters = newjetcol->parameters (); newparameters = copiedparameters;
    CopyParameters(jetcol->parameters(), newjetcol->parameters());
    newjetcol->setSubset(true);
    evt->addCollection(newjetcol, _outputjetcolName);
    CopyParameters(seccol->parameters(), newseccol->parameters());
    evt->addCollection(newseccol, _outputcolName);
    CopyParameters(relcol->parameters(), newjetrelcol->parameters());
    evt->addCollection(newjetrelcol, _outputjetrelcolName);
    CopyParameters(secRPcol->parameters(), newvtxRPcol->parameters());
    evt->addCollection(newvtxRPcol, _outputjetrelRPcolName);
    _totaltracks += recovery.GetStatistics();
    std::cout << "Total added tracks: " << _totaltracks << "\n";
    //*/
    // vector < MyReconstructedParticle * > * result = RefineVertices(maincol, seccol, v0col);
    // vector < ReconstructedParticle * > * result = RestoreOneTrackVertices(jetcol, relcol, seccol);

    // vector < ReconstructedParticle * > * result = RestoreVertices(primary, seccol);
    // vector < MyReconstructedParticle * > * result = RestoreVertices(jetcol, relcol, seccol);
    // LCCollection* misscol = evt->getCollection( _colMissedName );
    // LCCollection* missvtxcol =evt->getCollection( _colMissedVtxName );
    /*LCCollection* bstarcol = evt->getCollection( _colBStarName );
    //AnalyseSecondaries(maincol, misscol);
    LCCollection* trackrelcol = evt->getCollection( "MarlinTrkTracksMCTruthLink" );
    //CompareCollectionsRel(result, evt->getCollection(_colEGPName), bstarcol , relation); // relation
    //CompareCollections(result, misscol, bstarcol, missvtxcol);//*/
    //_Tree->Fill();
    // TestMethod2(evt->getCollection("MissedParticles"), seccol, evt->getCollection(_colRelName));
    //_TestTree->Fill();

    // WritePrimaries(opera, pricol, evt->getCollection(_colEGPName), relation);
    // WriteZombies(trashcol, maincol, evt->getCollection("MCParticlesSkimmed"), relation, trackrelcol);
    // WritePurgatory(pricol, evt->getCollection( "FinalVertex" ), maincol, evt->getCollection(_colEGPName), relation,
    // evt->getCollection("MCParticlesSkimmed")); _primaryTree->Fill(); _purgatoryTree->Fill();
    _trashTree->Fill();
    ClearVariables();
  } catch (DataNotAvailableException& e) {
    std::cout << "Whoops!....\n";
  }
}
void VertexRestorer::ClearVariables() {
  _primariesTotal = 0;
  _purgatoryTotal = 0;
  for (int i = 0; i < MAXP; i++) {
    _primeOffset[i] = -1.0;
    _primeType[i] = 0;
    _primeDeviation[i] = -1.0;
    _primeMomentum[i] = -1.0;
    _primeError[i] = -1.0;
    _primeTeorError[i] = -1.0;
    _primeChi2[i] = -1.0;
    _primeCostheta[i] = -2.0;
    _primeVtxHits[i] = -1.0;
    _primeFtdHits[i] = -1.0;
  }
  /*for (int i = 0; i < MAXN; i++)
  {
          _missedAngle[i] = -1.0;
          _fakeAngle[i] = -1.0;
  }*/
} /*
 void VertexRestorer::WriteZombies(LCCollection * trashcol, LCCollection * pfos, LCCollection * mccol, LCCollection *
 rel, LCCollection * trackrel)
 {
         LCRelationNavigator navigator(rel);
         _trashTotal = trashcol->getNumberOfElements();
         vector<MCParticle *> taken;
         for (int i = 0; i < _trashTotal; i++)
         {
                 Track * trashtrack = dynamic_cast< Track * >(trashcol->getElementAt(i));
                 _trashVtxHits[i] = trashtrack->getSubdetectorHitNumbers()[0];
                 _trashTpcHits[i] = trashtrack->getSubdetectorHitNumbers()[6];
                 ReconstructedParticle * trashreco = myTrackOperator.ReconstructParticle(trashtrack);
                 PrintParticle(trashreco);
                 MCParticle * trashmc = CompareToCollection(trashreco, mccol, trackrel, true);
                 if (!trashmc)
                 {
                         std::cout << "No mc particle!!\n";
                         _trashHasMCParticle[i] = 0;
                         continue;
                 }
                 _trashIsDublicate[i] = 0;
                 for (int j = 0; j < taken.size(); j++)
                 {
                         if (taken[j] == trashmc)
                         {
                                 std::cout << "Particle is dublicate\n";
                                 _trashIsDublicate[i] = 1;
                         }
                 }
                 taken.push_back(trashmc);
                  _trashHasMCParticle[i] = 1;
                 std::cout << "MC p: " << MathOperator::getModule(trashmc->getMomentum())
                           << " q: " << trashmc->getCharge()
                           << "\n";
                 vector<float> direction = MathOperator::getDirection(trashmc->getMomentum());
                 _trashCostheta[i] = std::cos( MathOperator::getAngles(direction)[1]);
                 vector< LCObject * > obj = navigator.getRelatedFromObjects(trashmc);
                 vector< float > weights = navigator.getRelatedFromWeights(trashmc);
                 float maxweight = 0.5;
                 ReconstructedParticle * trashreco2 = NULL;
                 for (int j = 0; j < obj.size(); j++)
                 {
                         if (weights[j] > maxweight)
                         {
                                 maxweight = weights[j];
                                 trashreco2 = dynamic_cast< ReconstructedParticle * >(obj[j]);
                         }
                 }
                 if (trashreco2)
                 {
                         std::cout << "Reco particle found:\n";
                         _trashIsReco[i] = 1;
                         PrintParticle(trashreco2);
                 }
                 else
                 {
                         _trashIsReco[i] = 0;
                         std::cout << "Reco particle NOT found!\n";
                 }
         }
         std::cout << "Total: " << _trashTotal << "\n";
 }

 void VertexRestorer::WritePrimaries(TrackOperator & opera, EVENT::LCCollection * pricol, EVENT::LCCollection * sec,
 EVENT::LCCollection * rel)
 {
         Vertex * primary = dynamic_cast< Vertex * >(pricol->getElementAt(0));
         double * pos = MathOperator::toDoubleArray(primary->getPosition(), 3);
         vector < ReconstructedParticle * > primaries = primary->getAssociatedParticle()->getParticles();
         _primariesTotal = primaries.size();
         _primeT = sqrt(primary->getPosition()[0] * primary->getPosition()[0] + primary->getPosition()[1] *
 primary->getPosition()[1]); _primeZ= primary->getPosition()[2]; for (unsigned int i = 0; i < primaries.size(); i++)
         {
                 //opera.test(primaries[i]);
                 if (sec && rel && CompareToCollection(primaries[i], sec,rel))
                 {
                         std::cout << "EXCLUDED!\n";
                         continue;
                 }
                 vector<float> direction = MathOperator::getDirection(primaries[i]->getMomentum());
                 double * start = opera.GetStartPoint(primaries[i]);
                 //PrintParticle(primaries[i]);
                 _primeType[i] = primaries[i]->getType();
                 _primeOffset[i] = MathOperator::getDistanceTo(pos, direction, start);//opera.GetOffset(primaries[i]);
                 //_primeOffset[i] = opera.GetOffset(primaries[i]);
                 _primeMomentum[i] = MathOperator::getModule(primaries[i]->getMomentum());

                 _primeError[i] = std::sqrt(opera.GetOffsetError(primaries[i], opera.GetStartPoint(primaries[i]),
 primary, _primeOffset[i])); _primeTeorError[i] = GetError(primaries[i]); _primeChi2[i] =
 primaries[i]->getTracks()[0]->getChi2()/(float)primaries[i]->getTracks()[0]->getNdf();
                 //_primeChi2[i] = primaries[i]->getTracks()[0]->getChi2();
                 _primeCostheta[i] = std::cos(MathOperator::getAngles(direction)[1]);
                 _primeVtxHits[i] = primaries[i]->getTracks()[0]->getSubdetectorHitNumbers()[0];
                 _primeFtdHits[i] = primaries[i]->getTracks()[0]->getSubdetectorHitNumbers()[5];
                 _primeEtdHits[i] = primaries[i]->getTracks()[0]->getSubdetectorHitNumbers()[10];
                 //std::cout << "Offset: " << _primeOffset[i] << "; error^2: " << _primeError[i] << '\n';
                 _primeDeviation[i] = _primeOffset[i] / _primeError[i];
                 //std::cout << "Offset: " << _primeOffset[i] << " deviation: " << _primeDeviation[i] << '\n';

         }

 }
 MCParticle * VertexRestorer::GetInterestingParent(MCParticle * daughter)
 {
         if (!daughter || daughter->getParents().size() < 1)
         {
                 //std::cout << "\tno parents\n";
                 return NULL;
         }
         MCParticle * parent = daughter->getParents()[0];
         int parentpdg = abs(parent->getPDG());
         if ((parentpdg > 200 && parentpdg < 300) ||
             (parentpdg > 10100 && parentpdg < 10225) ||
             (parentpdg > 110 && parentpdg < 120 )||
             parentpdg == 2212 ||
             parentpdg == 2112 ||
             parentpdg == 11 ||
             parentpdg == 13)
         {
                 //std::cout << "\tcontinue to search with " << parentpdg <<'\n';
                 return GetInterestingParent(parent);
         }
         //std::cout << "Satisfied with " << parentpdg <<'\n';
         return parent;
 }
 void VertexRestorer::WritePurgatory(LCCollection * pricol, LCCollection * sec, LCCollection * pfos, LCCollection
 *egprongs, LCCollection *rel, LCCollection *mccol)
 {
         Vertex * primary = dynamic_cast< Vertex * >(pricol->getElementAt(0));
         double * pos = MathOperator::toDoubleArray(primary->getPosition(), 3);
         vector < ReconstructedParticle * > primaries = primary->getAssociatedParticle()->getParticles();
         vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
         particles->reserve(particles->size() + primaries.size());
         particles->insert(particles->end(),primaries.begin(),primaries.end());
         int snumber = sec->getNumberOfElements();
         //int v0number = v0col->getNumberOfElements();
         for (int i = 0; i < snumber; i++)
         {
                 Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
                 vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
                 particles->reserve(particles->size() + secondaries.size());
                 particles->insert(particles->end(),secondaries.begin(),secondaries.end());
         }
         int pfocharge = 0;
         for (int i = 0; i < pfos->getNumberOfElements(); i++)
         {
                 ReconstructedParticle * purparticle = dynamic_cast< ReconstructedParticle * > (pfos->getElementAt(i));
                 if (std::abs(purparticle->getCharge()) < 0.9)
                 {
                         continue;
                 }
                 pfocharge++;
                 if (IsDublicate(purparticle, *particles))
                 {
                         continue;
                 }
                 if (sec && rel && CompareToCollection(purparticle, egprongs,rel))
                 {
                         std::cout << "EXCLUDED2!\n";
                         continue;
                 }
                 MCParticle * mcp = CompareToCollection(purparticle, mccol,rel);
                 if (mcp)
                 {
                         //std::cout << "Start looking for parent of " << abs(mcp->getPDG()) << "\n";
                         MCParticle * parent = GetInterestingParent(mcp);
                         _purgatoryParent[_purgatoryTotal] = (parent)? abs(parent->getPDG()): 0;
                         if (abs(mcp->getPDG()) == 13 && mcp->getParents().size() > 0 &&
 _purgatoryParent[_purgatoryTotal] == 92 )
                         {
                                 MCParticle * parent = mcp->getParents()[0];
                                 if (abs(parent->getPDG()) == 211)
                                 {
                                         _purgatoryParent[_purgatoryTotal] = 50000;
                                 }
                         }
                 }
                 else
                 {
                         _purgatoryParent[_purgatoryTotal] = 0;
                 }

                 vector<float> direction = MathOperator::getDirection(purparticle->getMomentum());
                 double * start = myTrackOperator.GetStartPoint(purparticle);
                 _purgatoryOffset[_purgatoryTotal] = MathOperator::getDistanceTo(pos, direction,
 start);//opera.GetOffset(primaries[i]); _purgatoryMomentum[_purgatoryTotal] =
 MathOperator::getModule(purparticle->getMomentum()); _purgatoryTeorError[_purgatoryTotal] = GetError(purparticle);
                 _purgatoryError[_purgatoryTotal] = std::sqrt(myTrackOperator.GetOffsetError(purparticle,
 myTrackOperator.GetStartPoint(purparticle), primary, _purgatoryOffset[_purgatoryTotal]));
                 _purgatoryCostheta[_purgatoryTotal] = std::cos(MathOperator::getAngles(direction)[1]);
                 _purgatoryVtxHits[_purgatoryTotal] = purparticle->getTracks()[0]->getSubdetectorHitNumbers()[0];
                 _purgatoryDeviation[_purgatoryTotal] = _purgatoryOffset[_purgatoryTotal] /
 _purgatoryError[_purgatoryTotal]; if ( purparticle->getStartVertex())
                 {
                         //std::cout << "Particle has start vertex  "<< purparticle->getStartVertex()->isPrimary() <<
 "\n";
                 }
                 else
                 {
                         //std::cout << "PRIMARY PARTICLE HAS NO VERTEX!!!\n";
                 }
                 _purgatoryTotal++;

         }
         std::cout << "Primaries: " << primaries.size()
                 << " primaries+secondaries: " << particles->size()
                 << " purgatories: " << _purgatoryTotal
                 << " PFOs: " << pfocharge
                 << "\n";

 }
 MCParticle * VertexRestorer::CompareToCollection(EVENT::ReconstructedParticle * particle, EVENT::LCCollection * prongs,
 EVENT::LCCollection * rel, bool byTrack)
 {
         LCRelationNavigator navigator(rel);
         vector< LCObject * > obj;// = navigator.getRelatedToObjects((byTrack)? particle->getTracks()[0]: particle);
         vector< float > weights;// = navigator.getRelatedToWeights((byTrack)? particle->getTracks()[0]: particle);
         if (byTrack)
         {
                 obj = navigator.getRelatedToObjects(particle->getTracks()[0]);
                 weights = navigator.getRelatedToWeights(particle->getTracks()[0]);
         }
         else
         {
                 obj = navigator.getRelatedToObjects(particle);
                 weights = navigator.getRelatedToWeights(particle);
         }
         if (obj.size() < 1)
         {
                 std::cout << "Particle was not generated!\n";
                 return false;
         }
         MCParticle * winner = NULL;
         float maxweight = 0.50;
         for (unsigned int i = 0; i < obj.size(); i++)
         {
                 if (weights[i] > maxweight)
                 {
                         winner = dynamic_cast< MCParticle * >(obj[i]);
                         maxweight = weights[i];
                 }
         }
         int numberOfProngs = prongs->getNumberOfElements();
         MCParticle * mcprong = NULL;
         for (int i = 0; i < numberOfProngs; i++)
         {
                 mcprong = dynamic_cast< MCParticle * >(prongs->getElementAt(i));
                 if (mcprong == winner)
                 {
                         return mcprong;
                 }
         }
         return NULL;
 }
 void VertexRestorer::CompareCollectionsRel(std::vector< MyReconstructedParticle * > * detected, EVENT::LCCollection *
 prongs, EVENT::LCCollection * bstar, EVENT::LCCollection *  rel)
 {
         _detectedTotal = detected->size();
         _missedTotal = prongs->getNumberOfElements();
         for (int i = 0; i < _missedTotal; i++)
         {
                 MCParticle * particle = dynamic_cast< MCParticle * >(prongs->getElementAt(i));
                 _allmissedMomentum[i] = MathOperator::getModule(particle->getMomentum());
         }
         _bstarTotal = bstar->getNumberOfElements();
         _missedDetected = 0;
         _bstarDetected = 0;
         _fakeDetected = 0;
         _purDetected = 0;
         bool byTracks = false;
         for (int i = 0; i < _detectedTotal; i++)
         {
                 float taken = false;
                 ReconstructedParticle * particle = detected->at(i)->Get();
                 if (CompareToCollection(particle, prongs, rel, byTracks))
                 {
                         taken = true;
                         std::cout << "Found a prong!\n";
                         _missedError[_missedDetected] = GetError(particle);
                         _missedMomentum[_missedDetected] = MathOperator::getModule(particle->getMomentum());
                         _missedDeviation[_missedDetected] =  detected->at(i)->GetOffset() /
 _missedError[_missedDetected]; _missedOffset[_missedDetected] =  detected->at(i)->GetOffset() ;
                         _missedAngle[_missedDetected] = detected->at(i)->GetAngle();
                         _missedObs[_missedDetected] = detected->at(i)->GetObservable();
                         _missedVtxHits[_missedDetected] = particle->getTracks()[0]->getSubdetectorHitNumbers()[0];
                         _missedFtdHits[_missedDetected] = particle->getTracks()[0]->getSubdetectorHitNumbers()[4];
                         _missedSecDeviation[_missedDetected] =  detected->at(i)->GetSecOffset();
                         _missedCostheta[_missedDetected] = detected->at(i)->GetCostheta();
                         _missedDetected++;
                 }
                 if (CompareToCollection(particle, bstar, rel, byTracks))
                 {
                         taken = true;
                         _bstarError[_bstarDetected] = GetError(particle);
                         _bstarMomentum[_bstarDetected] = MathOperator::getModule(particle->getMomentum());
                         _bstarDetected++;
                 }
                 if (!taken)
                 {
                         _fakeError[_fakeDetected] = GetError(particle);
                         _fakeMomentum[_fakeDetected] = MathOperator::getModule(particle->getMomentum());
                         _fakeDeviation[_fakeDetected] =  detected->at(i)->GetOffset() / _fakeError[_fakeDetected];
                         _fakeSecDeviation[_fakeDetected] =  detected->at(i)->GetSecOffset();
                         _fakeAngle[_fakeDetected] = detected->at(i)->GetAngle();
                         _fakeObs[_fakeDetected] = detected->at(i)->GetObservable();
                         _fakeCostheta[_fakeDetected] = detected->at(i)->GetCostheta();
                         _fakeVtxHits[_fakeDetected] = particle->getTracks()[0]->getSubdetectorHitNumbers()[0];
                         _fakeFtdHits[_fakeDetected] = particle->getTracks()[0]->getSubdetectorHitNumbers()[4];
                         _fakeDetected++;
                 }
                 if (!taken && !particle->getStartVertex())//&& !IsDublicate(particle,
 myPrimary->getAssociatedParticle()->getParticles()))
                 {
                         _purMomentum[_purDetected] = MathOperator::getModule(particle->getMomentum());
                         _purOffset[_purDetected] =  detected->at(i)->GetOffset();
                         _purAngle[_purDetected] =  detected->at(i)->GetAngle();
                         _purCostheta[_purDetected] =  detected->at(i)->GetCostheta();
                         std::cout << "Purgatory detected!\n";
                         _purDetected++;
                 }
         }
         std::cout << "Total detected: " << _detectedTotal
                 << "; missed total: " << _missedTotal
                 << "; missed detected: " << _missedDetected
                 << "; B* detected: " << _bstarDetected
                 << '\n';
 }

 void VertexRestorer::CompareCollections(vector< MyReconstructedParticle * > * detected, LCCollection * missed,
 LCCollection * bstar, LCCollection * missedvtx)
 {
         _missedTotal = missed->getNumberOfElements();
         for (int i = 0; i < _missedTotal; i++)
         {
                 ReconstructedParticle * particle = dynamic_cast< ReconstructedParticle * >(missed->getElementAt(i));
                 _allmissedMomentum[i] = MathOperator::getModule(particle->getMomentum());
         }
         _bstarTotal = bstar->getNumberOfElements();
         _detectedTotal = detected->size();

         _missedDetected = 0;
         _bstarDetected = 0;
         _fakeDetected = 0;

         for (int i = 0; i < _detectedTotal; i++)
         {
                 float taken = false;
                 ReconstructedParticle * detectedParticle = detected->at(i)->Get();
                 for (int j = 0; j < _missedTotal; j++)
                 {
                         ReconstructedParticle * particle = dynamic_cast< ReconstructedParticle * >
 (missed->getElementAt(j)); if (CompareParticles(particle, detectedParticle))
                         {
                                 taken = true;
                                 _missedError[_missedDetected] = GetError(detectedParticle);
                                 _missedMomentum[_missedDetected] = MathOperator::getModule(particle->getMomentum());
                                 _missedDeviation[_missedDetected] =  detected->at(i)->GetOffset() /
 _missedError[_missedDetected]; _missedSecDeviation[_missedDetected] =  detected->at(i)->GetSecOffset() /
 _missedError[_missedDetected]; _missedAngle[_missedDetected] = detected->at(i)->GetAngle(); _missedObs[_missedDetected]
 = detected->at(i)->GetObservable(); _missedCostheta[_missedDetected] = detected->at(i)->GetCostheta();
                                 _missedDetected++;
                         }
                 }
                 if (missedvtx && !taken)
                 {
                         int missedvtxTotal = missedvtx->getNumberOfElements();
                         for (int j = 0; j < missedvtxTotal; j++)
                         {
                                 ReconstructedParticle * particle = dynamic_cast< ReconstructedParticle * >
 (missedvtx->getElementAt(j)); if (CompareParticles(particle, detectedParticle))
                                 {
                                         taken = true;
                                         //_missedDetected++;
                                         _missedError[_missedDetected] = GetError(detectedParticle);
                                         _missedDeviation[_missedDetected] =  detected->at(i)->GetOffset() /
 _missedError[_missedDetected]; _missedAngle[_missedDetected] = detected->at(i)->GetAngle(); _missedObs[_missedDetected]
 = detected->at(i)->GetObservable(); _missedCostheta[_missedDetected] = detected->at(i)->GetCostheta();
                                         _missedSecDeviation[_missedDetected] =  detected->at(i)->GetSecOffset() /
 _missedError[_missedDetected]; _missedMomentum[_missedDetected++] = MathOperator::getModule(particle->getMomentum());
                                 }
                         }
                 }
                 for (int j = 0; j < _bstarTotal; j++)
                 {
                         MCParticle * particle = dynamic_cast< MCParticle * > (bstar->getElementAt(j));
                         if (CompareParticles(particle, detectedParticle))
                         {
                                 _bstarError[_bstarDetected] = GetError(detectedParticle);
                                 _bstarMomentum[_bstarDetected] = MathOperator::getModule(particle->getMomentum());
                                 //_bstarAngle[_bstarDetected] = detectedParticle->getParticleIDs()[0]->getLikelihood();
                                 _bstarDetected++;
                                 //_bstarDeviation[_bstarDetected] = GetDeviation(detectedParticle,
 MathOperator::toDoubleArray(detectedParticle->getStartVertex()->getPosition(), 3)); taken = true;
                         }
                 }
                 if (!taken)
                 {
                         _fakeError[_fakeDetected] = GetError(detectedParticle);
                         _fakeDeviation[_fakeDetected] =  detected->at(i)->GetOffset() / _fakeError[_fakeDetected];
                         _fakeSecDeviation[_fakeDetected] =  detected->at(i)->GetSecOffset() /
 _fakeError[_fakeDetected]; _fakeAngle[_fakeDetected] = detected->at(i)->GetAngle(); _fakeObs[_fakeDetected] =
 detected->at(i)->GetObservable(); _fakeCostheta[_fakeDetected] = detected->at(i)->GetCostheta();
                         _fakeMomentum[_fakeDetected++] = MathOperator::getModule(detectedParticle->getMomentum());
                 }
         }
         //_fakeDetected = _detectedTotal - _missedDetected - _bstarDetected;
         std::cout << "Total detected: " << _detectedTotal
                 << "; missed total: " << _missedTotal
                 << "; missed detected: " << _missedDetected
                 << "; B* detected: " << _bstarDetected
                 << '\n';

 }

 vector< MyReconstructedParticle * > * VertexRestorer::RefineVertices(LCCollection * main, LCCollection * sec,
 LCCollection * v0col)
 {
         vector < MyReconstructedParticle * > * result = new vector< MyReconstructedParticle * >();
         vector< ReconstructedParticle * > primaries;
         int mnumber = main->getNumberOfElements();
         for (int i = 0; i < mnumber; i++)
         {
                 primaries.push_back(dynamic_cast< ReconstructedParticle * > (main->getElementAt(i)));
         }
         vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
         int snumber = sec->getNumberOfElements();
         int v0number = v0col->getNumberOfElements();
         for (int i = 0; i < snumber; i++)
         {
                 Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
                 vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
                 particles->reserve(particles->size() + secondaries.size());
                 particles->insert(particles->end(),secondaries.begin(),secondaries.end());
         }
         for (int i = 0; i < v0number; i++)
         {
                 Vertex * secondary = dynamic_cast< Vertex * >(v0col->getElementAt(i));
                 vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
                 particles->reserve(particles->size() + secondaries.size());
                 //particles->insert(particles->end(),secondaries.begin(),secondaries.end());
         }
         for (int i = 0; i < snumber; i++)
         {
                 Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
                 //ParametrizeVertex(secondary);
                 vector< MyReconstructedParticle * > additional = AddParticles(primaries, secondary, particles);
                 for (unsigned int j = 0; j < additional.size(); j++)
                 {
                         if (!IsDublicate(additional[j], *result))
                         {
                                 result->push_back(additional[j]);
                         }
                 }
         }
         return result;
 }
 vector< float > VertexRestorer::ParametrizeVertex(const Vertex * sec)
 {
         const vector< ReconstructedParticle * > secondaries = sec->getAssociatedParticle()->getParticles();
         vector< float > result;
         if (secondaries.size() < 2)
         {
                 return result;
         }
         float mindprime = 1000;
         float maxdprime = 0;
         double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
         //std::cout << "Vertex position: " << MathOperator::getModule(sec->getPosition()) << " Chi2: " <<
 sec->getChi2() << "\n"; for (int i = 0; i < secondaries.size(); i++)
         {
                 double * secPos = myTrackOperator.GetStartPoint(secondaries[i]);
                 vector< float > secDir = MathOperator::getDirection(secondaries[i]->getMomentum());
                 double * vtxPos = MathOperator::toDoubleArray(sec->getPosition(),3);
                 vector<float> diff = MathOperator::getDirection(vtxPos, secPos);
                 double diif[3];
                 for (int m = 0; m < 3; m++)
                 {
                         diif[m] = diff[m];
                 }
                 float distance = MathOperator::getAngle(diif, secondaries[i]->getMomentum());
                 if (distance < mindprime)
                 {
                         mindprime = distance;
                 }
                 if (distance > maxdprime && distance < 2.0)
                 {
                         maxdprime = distance;
                 }

         }
         //std::cout << "\tdmin'= " << mindprime << "\tdmax'= " << maxdprime << "\n";
         result.push_back(mindprime);
         result.push_back(maxdprime);
         return result;

 }
 float VertexRestorer::GetMinDiffDistance(const ReconstructedParticle * candidate, const Vertex * sec, float & obs)
 {
         const vector< ReconstructedParticle * > secondaries = sec->getAssociatedParticle()->getParticles();
         double * candidatePos = myTrackOperator.GetStartPoint(candidate);
         vector< float > candidateDir = MathOperator::getDirection(candidate->getMomentum());
         float alpha = 0.0;
         float maxdistance = 0.0;
         float mindistance = 100.0;
         double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
         for (unsigned int i = 0; i < secondaries.size(); i++)
         {
                 double * secPos = myTrackOperator.GetStartPoint(secondaries[i]);
                 vector< float > secDir = MathOperator::getDirection(secondaries[i]->getMomentum());
                 float distance = MathOperator::getDistanceBtw(candidatePos, candidateDir, secPos, secDir);

                 //float angle = MathOperator::getAngle(candidate->getMomentum(), secondaries[i]->getMomentum());
                 if (distance < mindistance)
                 {
                         mindistance = distance;
                         obs = myTrackOperator.GetDprime(candidate, secondaries[i], primaryPosition);
                 }
         }

         return mindistance;
 }
 bool VertexRestorer::TakeParticle(EVENT::ReconstructedParticle * primary, const EVENT::Vertex * sec, vector<double> *
 parameters)
 {
                 if (std::abs(primary->getCharge() ) < 0.9)
                 {
                         return false;
                 }
                 const vector< ReconstructedParticle * > secondaries = sec->getAssociatedParticle()->getParticles();
                 double * trackPosition = myTrackOperator.GetStartPoint(primary);
                 float trackDistance = MathOperator::getModule(trackPosition);
                 if (trackDistance > 20.0)
                 {
                         return false;
                 }
                 //return true;
                 vector< float > direction = MathOperator::getDirection(primary->getMomentum());
                 double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
                 double * secondaryPosition = MathOperator::toDoubleArray(sec->getPosition(),3);
                 float primaryOffset = MathOperator::getDistanceTo(primaryPosition, direction, trackPosition);
                 float accuracy = GetError(primary);
                 //float accuracy = std::sqrt(myTrackOperator.GetOffsetError(primary, trackPosition, myPrimary,
 primaryOffset)); float secondaryOffset = MathOperator::getDistanceTo(secondaryPosition, direction, trackPosition);
                 vector<float> diff = MathOperator::getDirection(secondaryPosition, trackPosition);
                 double diif[3];
                 for (int m = 0; m < 3; m++)
                 {
                         diif[m] = diff[m];
                 }
                 float angle = MathOperator::getAngle(diif, primary->getMomentum());
                 float secOffset =  MathOperator::getDistanceTo(secondaryPosition, direction, trackPosition);
                 float costheta = std::cos(MathOperator::getAngles(direction)[1]); //std::abs(
 std::cos(MathOperator::getAngles(secDiraction)[1] ) ); float dprime = myTrackOperator.GetDprime(primary,
 secondaries[0], primaryPosition);
                 //float l =  GetMinDiffDistance(primary, sec, dprime);// std::abs( distance / secondaryOffset -0.5 ) *
 cos; vector< float > limits = ParametrizeVertex(sec); int vtxhits =
 primary->getTracks()[0]->getSubdetectorHitNumbers()[0]; int ftdhits =
 primary->getTracks()[0]->getSubdetectorHitNumbers()[5]; float anglecut = 0.1;//limits[1]*5.0+.05;
                 //bool result = //(primaryOffset /accuracy  > 1. + 40.*angle || angle < 0.005)
                 bool result = (primaryOffset /accuracy  > 0.5 + 15.0 * sqrt(angle) || angle < 0.005)
                         && (vtxhits> 3  || angle < 0.001) // ||( abs(costheta) < 0.9 && vtxhits+ftdhits > 5 ) || angle
 < 0.001)
                         && angle < anglecut; // angle < anglecut;// + 0.03 * sine;
                 if (result)
                 {
                         std::cout << "Found a track with offset " << primaryOffset
                                 << ", error " << accuracy//GetError(primary)
                                 << ", Angle " << angle //GetError(primary)
                                 << ", DISTANCE " << dprime - MathOperator::getModule(secondaryPosition)
 //distance//GetError(primary)
                                 << ", costheta " << costheta
                                 << ", vtxhits " << vtxhits
                                 //<< ", ftdhits " << ftdhits
                                 <<" :\n";
                         if (parameters)
                         {
                                 parameters->push_back(angle);//dprime - MathOperator::getModule(secondaryPosition));
                                 parameters->push_back(limits[1]);
                                 parameters->push_back(primaryOffset);
                                 parameters->push_back(secOffset);
                                 parameters->push_back(costheta);
                         }
                 }
                 else
                 {
                         if (parameters)
                         {
                                 delete parameters;
                         }
                 }
                 return result;
 }

 vector< MyReconstructedParticle * > VertexRestorer::AddParticles(const std::vector< EVENT::ReconstructedParticle * > &
 pri, EVENT::Vertex * sec, const std::vector< EVENT::ReconstructedParticle * > * toCompare)
 {
         vector< MyReconstructedParticle * > result;
         if (sec->getAssociatedParticle()->getParticles().size() < 2)
         {
                 return result;
         }
         const vector< ReconstructedParticle * > * particles = (toCompare)? toCompare :
 &(sec->getAssociatedParticle()->getParticles()); for (unsigned int i = 0; i < pri.size(); i++)
         {
                 vector<double> * parameters = new vector<double> ();
                 ReconstructedParticle * primary = pri[i];
                 if (TakeParticle(primary, sec, parameters))
                 {
                         PrintParticle(primary);
                         bool found = false;
                         for (unsigned int j = 0; j < particles->size(); j++)
                         {
                                 ReconstructedParticle * secondary = particles->at(j);
                                 if (CompareParticles(secondary, primary))
                                 {
                                         found = true;
                                         break;
                                 }
                         }

                         if (!IsDublicate(primary, *particles))
                         {
                                 std::cout << "The particle is unique!\n";
                                 //float errorvtx = (sec->getChi2() < 0.00001)? 0.00001: sec->getChi2();
                                 MyReconstructedParticle * newparticle = new MyReconstructedParticle(primary);
                                 newparticle->SetAngle(parameters->at(0));
                                 newparticle->SetObservable(parameters->at(1));
                                 newparticle->SetOffset(parameters->at(2));
                                 newparticle->SetSecOffset(parameters->at(3));
                                 newparticle->SetCostheta(parameters->at(4));
                                 result.push_back(newparticle);
                         }
                 }
         }
         return result;

 }
 void VertexRestorer::AddParticleID(EVENT::ReconstructedParticle * particle, float theta)
 {
         ParticleIDImpl * pid = new ParticleIDImpl();
         pid->setAlgorithmType(42);
         pid->setLikelihood(theta);
         particle->addParticleID(pid);
 }
 vector< ReconstructedParticle * > * VertexRestorer::RestoreVerticesPFO(LCCollection * main, EVENT::LCCollection * sec)
 {
         vector < ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
         vector< ReconstructedParticle * > primaries;
         int mnumber = main->getNumberOfElements();
         for (int i = 0; i < mnumber; i++)
         {
                 primaries.push_back(dynamic_cast< ReconstructedParticle * > (main->getElementAt(i)));
         }
         int snumber = sec->getNumberOfElements();
         vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
         for (int i = 0; i < snumber; i++)
         {
                 Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
                 vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
                 particles->reserve(particles->size() + secondaries.size());
                 particles->insert(particles->end(),secondaries.begin(),secondaries.end());
         }
         for (int i = 0; i < snumber; i++)
         {
                 Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
                 vector< ReconstructedParticle * > additional = GetAdditionalParticles(primaries, secondary, particles);
                 for (unsigned int j = 0; j < additional.size(); j++)
                 {
                         if (!IsDublicate(additional[j], *result))
                         {
                                 result->push_back(additional[j]);
                         }
                 }
         }
         return result;

 }

 void VertexRestorer::AnalyseSecondaries(EVENT::LCCollection * main, EVENT::LCCollection * missed)
 {
         int mnumber = main->getNumberOfElements();
         _missedTotal = missed->getNumberOfElements();
         vector< ReconstructedParticle * > lost;
         for ( int i = 0; i < _missedTotal; i++)
         {
                 ReconstructedParticle * secondary = dynamic_cast< ReconstructedParticle * > (missed->getElementAt(i));
                 bool found = false;
                 for (int j = 0; j < mnumber; j++)
                 {
                         ReconstructedParticle * primary = dynamic_cast< ReconstructedParticle * >
 (main->getElementAt(j)); if (CompareParticles(secondary, primary))
                         {
                                 found = true;
                                 break;
                         }
                 }
                 if (!found)
                 {
                         lost.push_back(secondary);
                 }
         }
         _missedRecoTotal = lost.size();
         for (int i = 0; i < _missedRecoTotal; i++)
         {
                 _missedRecoMomentum[i] = MathOperator::getModule(lost[i]->getMomentum());
                 vector<float> direction = MathOperator::getDirection(lost[i]->getMomentum());
                 _missedRecoTheta[i] = MathOperator::getAngles(direction)[1];
         }
 }
 bool VertexRestorer::CompareParticles(const ReconstructedParticle * particle1, const ReconstructedParticle * particle2,
 bool verbose)
 {
         if (particle1 == particle2)
         {
                 //std::cout << "Equal pointers!\n";
                 return true;
         }
         //return false;
         if (particle1->getCharge() * particle2->getCharge() < 0.0)
         {
                 return false;
         }
         float angle = MathOperator::getAngle(particle1->getMomentum(), particle2->getMomentum());
         if (angle > 0.02)
         {
                 return false;
         }
         float recomodule = MathOperator::getModule(particle2->getMomentum());
         float mcmodule = MathOperator::getModule(particle1->getMomentum());
         float ratio = (1 - recomodule/mcmodule > 0.0) ? 1 - recomodule/mcmodule : recomodule/mcmodule - 1.0;
         if (ratio > 0.1)
         {
                 return false;
         }
         return true;
 }
 bool VertexRestorer::CompareParticles(const MCParticle * particle1, const ReconstructedParticle * particle2)
 {
         if (particle1->getCharge() * particle2->getCharge() < 0.0)
         {
                 return false;
         }
         float angle = MathOperator::getAngle(particle1->getMomentum(), particle2->getMomentum());
         if (angle > 0.02)
         {
                 return false;
         }
         float recomodule = MathOperator::getModule(particle2->getMomentum());
         float mcmodule = MathOperator::getModule(particle1->getMomentum());
         float ratio = (1 - recomodule/mcmodule > 0.0) ? 1 - recomodule/mcmodule : recomodule/mcmodule - 1.0;
         if (ratio > 0.02)
         {
                 return false;
         }
         return true;
 }
 vector < EVENT::ReconstructedParticle * > * VertexRestorer::RestoreVertices(Vertex * primary, LCCollection * sec)
 {
         //vector < EVENT::Vertex * > * result = new vector< EVENT::Vertex* >();
         vector < ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
         vector< ReconstructedParticle * > primaries = primary->getAssociatedParticle()->getParticles();
         int number = sec->getNumberOfElements();
         vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
         for (int i = 0; i < number; i++)
         {
                 Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
                 vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
                 particles->reserve(particles->size() + secondaries.size());
                 particles->insert(particles->end(),secondaries.begin(),secondaries.end());
         }

         for (int i = 0; i < number; i++)
         {
                 Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
                 vector< ReconstructedParticle * > additional = GetAdditionalParticles(primaries, secondary, particles);
                 for (unsigned int j = 0; j < additional.size(); j++)
                 {
                         if (!IsDublicate(additional[j], *result))
                         {
                                 result->push_back(additional[j]);
                         }
                 }
         }
         return result;
 }
 vector< ReconstructedParticle * > * VertexRestorer::RestoreOneTrackVertices(LCCollection * jetcol, LCCollection * rel,
 LCCollection * sec)
 {
         vector < ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
         int number = jetcol->getNumberOfElements();
         LCRelationNavigator navigator(rel);
         PIDHandler pidh(jetcol);
         int alid = pidh.getAlgorithmID("lcfiplus");
         vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
         for (int i = 0; i < sec->getNumberOfElements(); i++)
         {
                 Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
                 vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
                 particles->reserve(particles->size() + secondaries.size());
                 particles->insert(particles->end(),secondaries.begin(),secondaries.end());
         }

         for (int i = 0; i < number; i++)
         {
                 ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle *>(jetcol->getElementAt(i));
                 int nvtx = navigator.getRelatedToObjects(jet).size();
                 const ParticleID& pid = pidh.getParticleID(jet,alid);
                 vector<float> params = pid.getParameters();
                 float btag = params[pidh.getParameterIndex(alid,"BTag")];
                 float ctag = params[pidh.getParameterIndex(alid,"CTag")];
                 if (btag < 0.3 || nvtx == 2 || ctag > 0.3)
                 {
                         continue;
                 }
                 std::cout << "Jet energy: " << jet->getEnergy()
                           << " b-tag: " << btag
                           << " c-tag: " << ctag
                           << " # of vtx: " << nvtx
                           <<  '\n';
                 vector< ReconstructedParticle * > additional = GetAdditionalParticles(jet->getParticles(), particles);
                 result->reserve(result->size() + additional.size());
                 result->insert(result->end(), additional.begin(),additional.end());

         }
         return result;
 }
 vector< ReconstructedParticle * > VertexRestorer::GetAdditionalParticles(const vector< ReconstructedParticle * > & pri,
 const vector< ReconstructedParticle * > * particles )
 {
         vector< ReconstructedParticle * > result;
         for (unsigned int i = 0; i < pri.size(); i++)
         {
                 ReconstructedParticle * primaryTrack = pri[i];
                 if (abs(primaryTrack->getCharge())  < 0.9)
                 {
                         continue;
                 }
                 vector< float > direction = MathOperator::getDirection(primaryTrack->getMomentum());
                 double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
                 double * trackPosition = myTrackOperator.GetStartPoint(primaryTrack);
                 float primaryOffset = MathOperator::getDistanceTo(primaryPosition, direction, trackPosition);
                 float accuracy = std::sqrt(myTrackOperator.GetOffsetError(primaryTrack, trackPosition, myPrimary,
 primaryOffset)); if (primaryOffset / accuracy > 5.0
                     && primaryOffset < 3.0)
                 {
                         std::cout << "Found a track with offset " << primaryOffset
                                 << ", deviation " << primaryOffset / accuracy
                                 << ", distance " << MathOperator::getModule(trackPosition)
                                 << ", error*e-4 " << accuracy
                                 <<" :\n";
                         PrintParticle(primaryTrack);
                         bool found = false;
                         for (unsigned int j = 0; j < particles->size(); j++)
                         {
                                 ReconstructedParticle * secondary = particles->at(j);
                                 if (CompareParticles(secondary, primaryTrack))
                                 {
                                         found = true;
                                         break;
                                 }
                         }
                         if (!found)
                         {
                                 std::cout << "The particle is unique!\n";
                                 result.push_back(primaryTrack);
                         }
                 }
         }
         return result;	    //&& primaryOffset / accuracy < 3.0

 }
 vector < MyReconstructedParticle * > * VertexRestorer::RestoreVertices(LCCollection * jetcol, LCCollection * rel,
 LCCollection * secvtx,  LCCollection * trash)
 {
         vector < MyReconstructedParticle * > * result = new vector< MyReconstructedParticle * >();
         int number = jetcol->getNumberOfElements();
         LCRelationNavigator navigator(rel);
         int snumber = secvtx->getNumberOfElements();
         vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
         vector< Vertex * > tagged;

         for (int i = 0; i < snumber; i++)
         {
                 Vertex * secondary = dynamic_cast< Vertex * >(secvtx->getElementAt(i));
                 tagged.push_back(secondary);
                 vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
                 particles->reserve(particles->size() + secondaries.size());
                 particles->insert(particles->end(),secondaries.begin(),secondaries.end());
         }
         vector< ReconstructedParticle * >  resurrected;
         if (trash)
         {
                 //std::cout << "Resurrected particles: \n";
                 int tracknumber = trash->getNumberOfElements();
                 for (int i = 0; i < tracknumber; i++)
                 {
                         Track * track = dynamic_cast<Track *>(trash->getElementAt(i));
                         ReconstructedParticle * particle = myTrackOperator.ReconstructParticle(track);
                         //PrintParticle(particle);
                         resurrected.push_back(particle);
                 }
                 for (int i = 0; i < resurrected.size(); i++)
                 {
                         for (int j = 0; j < i; j++)
                         {
                                 if (j != i && CompareParticles(resurrected[i],resurrected[j]))
                                 {
                                          //std::cout <<"Found same zombies i: "<< i << " j: " << j <<"!\n";
                                          //PrintParticle(resurrected[i]);
                                          //PrintParticle(resurrected[j]);
                                 }
                         }
                 }
                 std::cout << "\n";
         }
         for (int i = 0; i < number; i++)
         {
                 ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle *>(jetcol->getElementAt(i));
                 int nvtx = navigator.getRelatedToObjects(jet).size();
                 std::cout << "Jet energy: " << jet->getEnergy()
                          // << " b-tag: " << btag
                          // << " c-tag: " << ctag
                           << " # of vtx: " << nvtx
                           <<  '\n';
                 vector< LCObject * > objs = navigator.getRelatedToObjects(jet);
                 for (int j = 0; j < nvtx; j++)
                 {
                         Vertex * vertex = dynamic_cast< Vertex * >(objs[j]);
                         bool tag = false;
                         for (int i = 0; i < tagged.size(); i++)
                         {
                                 if (tagged[i] == vertex)
                                 {
                                         tag = true;
                                 }
                         }
                         if (!tag)
                         {
                                 continue;
                         }
                         vector< ReconstructedParticle * > toInject;
                         toInject.reserve(jet->getParticles().size() + resurrected.size());
                         toInject.insert(toInject.end(), jet->getParticles().begin(), jet->getParticles().end());
                         //toInject.insert(toInject.end(), resurrected.begin(), resurrected.end());

                         vector< MyReconstructedParticle * > additional = AddParticles(toInject, vertex, particles);
                         for (unsigned int j = 0; j < additional.size(); j++)
                         {
                                 if (!IsDublicate(additional[j], *result))
                                 {
                                         result->push_back(additional[j]);
                                 }
                         }
                 }
         }
         return result;
 }

 vector< ReconstructedParticle * > VertexRestorer::GetAdditionalParticles(const vector< ReconstructedParticle * > & pri,
 Vertex * sec, const vector< ReconstructedParticle * > * toCompare)
 {
         vector< ReconstructedParticle * > result;
         if (sec->getAssociatedParticle()->getParticles().size() < 2)
         {
                 return result;
         }
         const vector< ReconstructedParticle * > * particles = (toCompare)? toCompare :
 &(sec->getAssociatedParticle()->getParticles()); for (unsigned int i = 0; i < pri.size(); i++)
         {
                 ReconstructedParticle * primaryTrack = pri[i];
                 if (abs(primaryTrack->getCharge())  < 0.9)
                 {
                         continue;
                 }
                 vector< float > direction = MathOperator::getDirection(primaryTrack->getMomentum());
                 double * pos = MathOperator::toDoubleArray(sec->getPosition(), 3);
                 double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
                 double * trackPosition = myTrackOperator.GetStartPoint(primaryTrack);
                 float primaryOffset = MathOperator::getDistanceTo(primaryPosition, direction, trackPosition);
                 float secondaryOffset = MathOperator::getDistanceTo(pos, direction, trackPosition);
                 float accuracy = std::sqrt(myTrackOperator.GetOffsetError(primaryTrack, trackPosition, myPrimary,
 primaryOffset)); vector<float> diff = MathOperator::getDirection(pos, trackPosition); double diif[3]; for (int m = 0; m
 < 3; m++)
                 {
                         diif[m] = diff[m];
                 }
                 float tan  = secondaryOffset / MathOperator::getModule(diif);
                 float angle = MathOperator::getAngle(primaryTrack->getMomentum(), diif);
                 if (primaryOffset / accuracy < 3.0
                     //&& primaryOffset / accuracy < 3.0

                     //&& secondaryOffset < 5.0
                     && angle < 0.03)
                 {
                         std::cout << "Found a track with offset " << primaryOffset
                                 << ", deviation " << primaryOffset / accuracy
                                 << ", angle " << angle
                                 << ", distance " << MathOperator::getModule(trackPosition)
                                 << ", error*e-4 " << accuracy
                                 <<" :\n";
                         PrintParticle(primaryTrack);
                         bool found = false;
                         for (unsigned int j = 0; j < particles->size(); j++)
                         {
                                 ReconstructedParticle * secondary = particles->at(j);
                                 if (CompareParticles(secondary, primaryTrack))
                                 {
                                         found = true;
                                         break;
                                 }
                         }
                         if (!found)
                         {
                                 std::cout << "The particle is unique!\n";
                         //	if (!IsDublicate(primaryTrack,result))
                                 {
                                         result.push_back(CopyParticle(primaryTrack,sec, tan));
                                 }
                         }
                 }
         }
         return result;
 }

 ReconstructedParticle * VertexRestorer::CopyParticle(const ReconstructedParticle * particle, const Vertex * vertex,
 float theta)
 {
         ReconstructedParticleImpl * nouveau = new ReconstructedParticleImpl();
         VertexImpl * vimp = new VertexImpl();
         nouveau->setMomentum(particle->getMomentum());
         nouveau->setMass(particle->getMass());
         nouveau->setCharge(particle->getCharge());
         nouveau->setEnergy(particle->getEnergy());
         vimp->setPosition(vertex->getPosition());
         vimp->setPrimary(vertex->isPrimary());
         nouveau->setStartVertex(vimp);
         nouveau->addTrack(particle->getTracks()[0]); //CRUNCH
         AddParticleID(nouveau, theta);
         return nouveau;
 }
 bool VertexRestorer::IsDublicate(const ReconstructedParticle * particle, const vector< ReconstructedParticle * > &
 data)
 {
         bool dublicate = false;
         for (unsigned int j = 0; j < data.size(); j++)
         {
                 if (CompareParticles(particle, data[j]))
                 {
                         //std::cout << "Dublicate found!!!!\n";
                         dublicate = true;
                         break;
                 }
         }
         return dublicate;
 }
 bool VertexRestorer::IsDublicate( MyReconstructedParticle * particle, vector< MyReconstructedParticle * > & data)
 {
         bool dublicate = false;
         for (unsigned int j = 0; j < data.size(); j++)
         {
                 if (CompareParticles(particle->Get(), data[j]->Get()))
                 {
                         //std::cout << "Dublicate found!!!!\n";
                         dublicate = true;
                         break;
                 }
         }
         return dublicate;
 }
 float VertexRestorer::GetError(const ReconstructedParticle * particle)
 {
         if (!particle || particle->getTracks().size() < 1)
         {
                 std::cout << "The particle is null or 0 tracks!\n";
                 return 0.0;
         }
         float p = MathOperator::getModule(particle->getMomentum());
         vector<float> direction = MathOperator::getDirection(particle->getMomentum());
         vector<float> angles = MathOperator::getAngles(direction);
         float accuracy = sqrt(_aParameter*_aParameter + _bParameter*_bParameter /( p * p *
 pow(sin(angles[1]), 4.0/3.0)) ); return accuracy;
 }
 void VertexRestorer::TestMethod(EVENT::LCCollection * prongs, EVENT::LCCollection * pfo, EVENT::LCCollection * rel)
 {
         _nlosttest = 0;
         _ngoodtest = 0;
         _nmultitest = 0;
         _ntotaltest = 0;
         vector< MCParticle * > verifiedProngs;
         for (int i = 0; i < pfo->getNumberOfElements(); i++)
         {
                 ReconstructedParticle * recoparticle = dynamic_cast<ReconstructedParticle *>(pfo->getElementAt(i));
                 if (std::abs(recoparticle->getCharge()) < 0.01)
                 {
                         continue;
                 }
                 int times = 0;
                 MCParticle * mcprong = CompareToCollection(recoparticle, prongs, rel);
                 if (mcprong)
                 {
                         verifiedProngs.push_back(mcprong);
                         for (unsigned int j = 0; j < verifiedProngs.size(); j++)
                         {
                                 if (verifiedProngs[j] == mcprong)
                                 {
                                         times++;
                                 }
                         }
                 }
                 vector<float> direction = MathOperator::getDirection(recoparticle->getMomentum());
                 float cos = std::cos(MathOperator::getAngles(direction)[1]);
                 if (times == 1)
                 {
                         std::cout << "A good one found\n";
                         _costhetaTestedGood[_ngoodtest] = cos;
                         _ngoodtest++;
                 }
                 if (times == 0)
                 {
                         std::cout << "A lost one found\n";

                         _costhetaTestedLost[_nlosttest] = cos;
                         _nlosttest++;
                 }
                 if (times>1)
                 {
                         std::cout << "Several tracks found\n";
                         _costhetaTestedMultiple[_nmultitest] = cos;
                         _nmultitest++;
                 }

         }
         for (int i = 0; i < prongs->getNumberOfElements(); i++)
         {
                 MCParticle * mcparticle = dynamic_cast<MCParticle *>(prongs->getElementAt(i));
                 vector<float> direction = MathOperator::getDirection(mcparticle->getMomentum());
                 float cos = std::cos(MathOperator::getAngles(direction)[1]);
                 _costhetaTestedTotal[_ntotaltest] = cos;
                 _ntotaltest++;
         }
         return;
 }
 void VertexRestorer::TestMethod2(EVENT::LCCollection * prongs, EVENT::LCCollection * pfo, EVENT::LCCollection * rel)
 {
         double momentum[3];
         momentum[0] = 30.0;
         momentum[1] = 30.0;
         momentum[2] = 0.0;
         double zero[3];
         zero[0] = 0.0;
         zero[1] = 0.0;
         zero[2] = 0.0;
         double point[3];
         point[0] = 6.0;
         point[1] = 1.0;
         point[2] = 0.0;
         vector< float > direction = MathOperator::getDirection(momentum);
         float offset = MathOperator::getDistanceTo(point, direction, zero);
         std::cout << "Our offset is " << offset << "mm\n";
         int vtxnumber = pfo->getNumberOfElements();

         _ntest = 0;
         for (int i = 0; i < vtxnumber; i++)
         {
                 Vertex * recovtx = dynamic_cast< Vertex * > (pfo->getElementAt(i));
                 const vector< ReconstructedParticle * > particles = recovtx->getAssociatedParticle()->getParticles();
                 float minDistance = 10.0;
                 for (int j = 0; j < particles.size(); j++)
                 {
                         double * secPos = myTrackOperator.GetStartPoint(particles[j]);
                         vector< float > secDir = MathOperator::getDirection(particles[j]->getMomentum());
                         double * vtxPos = MathOperator::toDoubleArray(recovtx->getPosition(),3);
                         for (int k = 0; k < j; k++)
                         {
                                 if (k == j)
                                 {
                                         continue;
                                 }
                                 double * secPos1 = myTrackOperator.GetStartPoint(particles[k]);
                                 vector< float > secDir1 = MathOperator::getDirection(particles[k]->getMomentum());
                                 //float distance = MathOperator::getAngle(particles[j]->getMomentum(),
 particles[i]->getMomentum()); float distance = MathOperator::getDistanceBtw(secPos, secDir, secPos1, secDir1);
                                 std::cout << "i: " << i << " j: " << j << " distance: " << distance << "\n";
                                 if (distance < minDistance)
                                 {
                                         minDistance = distance;
                                 }
                         }
                 }
                 _testparameter[_ntest] = recovtx->getChi2();
                 _testDistance[_ntest++] = minDistance;
         }*/
/*int prongnumber = pfo->getNumberOfElements();
vector< ReconstructedParticle * > particles;
for (int i = 0; i < prongnumber; i++)
{
        MCParticle * mcparticle = dynamic_cast<MCParticle *>(prongs->getElementAt(i));
        LCRelationNavigator navigator(rel);
        vector< LCObject * > obj = navigator.getRelatedFromObjects(mcparticle);
        vector< float > weights = navigator.getRelatedFromWeights(mcparticle);
        if (obj.size() < 1)
        {
                std::cout << "Particle was not generated!\n";
                return false;
        }
        MCParticle * winner = NULL;
        float maxweight = 0.0;
        for (unsigned int j = 0; j < obj.size(); j++)
        {
                if (std::abs(dynamic_cast< ReconstructedParticle * >(obj[j])->getCharge()) < 0.9)
                {
                        continue;
                }
                if (weights[j] > maxweight)
                {
                        winner = dynamic_cast< ReconstructedParticle * >(obj[j]);
                        maxweight = weights[j];
                }
        }
        if (winner)
        {

        }
//}
}*/
void VertexRestorer::TestCollections(EVENT::LCCollection* tagged, EVENT::LCCollection* pfo) {
  int pfonumber = pfo->getNumberOfElements();
  const vector<ReconstructedParticle*> primaries = myPrimary->getAssociatedParticle()->getParticles();
  for (int i = 0; i < pfonumber; i++) {
  }
}

void VertexRestorer::check(LCEvent* evt) {}

void VertexRestorer::end() {
  hfile->Write();
  hfile->Close();
}
} // namespace TTbarAnalysis
