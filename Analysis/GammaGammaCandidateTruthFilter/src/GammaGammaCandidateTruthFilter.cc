#include "GammaGammaCandidateTruthFilter.h"

#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>


GammaGammaCandidateTruthFilter aGammaGammaCandidateTruthFilter ;


GammaGammaCandidateTruthFilter::GammaGammaCandidateTruthFilter() : Processor("GammaGammaCandidateTruthFilter") {
  
  // modify processor description
  _description = "GammaGammaCandidateTruthFilter does whatever it does ..." ;

  
  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "Printing" , 
                              "Print certain messages"  ,
                              _printing,
                               (int)1 ) ;

  std::string ggResonanceName = "Pi0";
  registerProcessorParameter( "GammaGammaResonanceName" , 
                              "Particle decaying to Gamma Gamma"  ,
                              _ggResonanceName,
                              ggResonanceName) ;

  registerInputCollection( LCIO::LCRELATION,
			   "MCTruth2RecoLinkCollectionName" , 
			   "true - reco relation collection"  ,
			   _trueToReco,
			   std::string("MCTruthRecoLink") ) ;


  registerInputCollection( LCIO::LCRELATION,
			   "Reco2MCTruthLinkCollectionName" , 
			   "reco - true relation collection"  ,
			   _recoToTrue,
			   std::string("RecoMCTruthLink") ) ;


  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the MCParticle input collection"  ,
			   _mcParticleCollectionName ,
			   std::string("MCParticle") ) ;


  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "GammaGammaCandidateCollection" , 
			   "Name of the gamma-gamma candidate input collection"  ,
			   _gammaGammaParticleCollectionName ,
			   std::string("GammaGammaCandidates") ) ;


  std::string outputParticleCollectionName = "TrueGammaGammaMesons";
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                             "OutputParticleCollectionName" , 
                             "Output Particle Collection Name "  ,
                             _outputParticleCollectionName,
                             outputParticleCollectionName);

}


void GammaGammaCandidateTruthFilter::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;

  // usually a good idea to
  printParameters() ;
  nEvt = 0;

//  gROOT->ProcessLine("#include <vector>");
}

void GammaGammaCandidateTruthFilter::processRunHeader( LCRunHeader* run) { 
    
}
 

void GammaGammaCandidateTruthFilter::processEvent( LCEvent * evt ) {

  // Make a new vector of particles
  LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  recparcol->setSubset(true);      

  streamlog_out(MESSAGE) << " start processing event " << std::endl;

  int typeRequired;
  if(_ggResonanceName=="Pi0")typeRequired = 111;
  if(_ggResonanceName=="Eta")typeRequired = 221;
  if(_ggResonanceName=="EtaPrime")typeRequired = 331;

/*
  int nMCPi0 = 0;
  int nRecoPi0 = 0;  

  streamlog_out(DEBUG) << " iterator and navigator " << std::endl;
  LCIterator<MCParticle> mcpIt( evt, _mcParticleCollectionName ) ;
  streamlog_out(DEBUG) << " got mcpIt with length of " << mcpIt.size() << std::endl;
  LCRelationNavigator mc2recoNav(evt->getCollection( _trueToReco )); 
  streamlog_out(DEBUG) << " got mc2recoNav from " << mc2recoNav.getFromType() << " to " << mc2recoNav.getToType() << std::endl;
  
  //----------------------------------------------------------------------------------------------------------------------------
  // loop over MCParticles
  //----------------------------------------------------------------------------------------------------------------------------
  
  int index_mcp = -1;   // Should be in synch with dumpevent utility

  int nMCEta = 0;
  int nMCEtaPrime = 0;
  double ESum_MCP = 0.0;

  int nRecoEta = 0;
  int nRecoEtaPrime = 0;

  while( MCParticle* mcp = mcpIt.next()  ) {
    
    if (!mcp) continue; 
    index_mcp++;
    
    bool keep = false;
    const EVENT::MCParticleVec& daughters = mcp->getDaughters();
    // keep pi0/eta/etaprime which decays to gammagamma if production vertex within 10 cm of nominal interaction point
    if (mcp->getPDG() == 111 || mcp->getPDG() == 221 || mcp->getPDG() == 331) {
      if (daughters.size() == 2 && daughters[0]->getPDG() == 22 && daughters[1]->getPDG() == 22){
          double rVertexSquared = 0.0;
          for(int i=0; i<=2; i++){
              rVertexSquared += pow(mcp->getVertex()[i],2);
          }
          if(rVertexSquared < 10000.0)keep = true;
      }
    }
    
    if (keep) {   
    
      streamlog_out(DEBUG) << " Found GammaGamma decaying parent MCParticle " << "Index = " << index_mcp << " pdgID = " << mcp->getPDG() 
                           << ", genstat = " << mcp->getGeneratorStatus() 
                           << ", E = " << mcp->getEnergy() << std::endl;
      streamlog_out(DEBUG) << " Momentum = " << mcp->getMomentum()[0] << " " << mcp->getMomentum()[1] << " " << mcp->getMomentum()[2] << std::endl;
      streamlog_out(DEBUG) << " Vertex =   " << mcp->getVertex()[0] << " " << mcp->getVertex()[1] << " " << mcp->getVertex()[2] << std::endl;
      streamlog_out(DEBUG) << " Endpoint = " << mcp->getEndpoint()[0] << " " << mcp->getEndpoint()[1] << " " << mcp->getEndpoint()[2] << std::endl;  
    
      nMCPi0++;
      ESum_MCP += mcp->getEnergy();
      if(mcp->getPDG() == 221)nMCEta++;
      if(mcp->getPDG() == 331)nMCEtaPrime++;
    
      gear::Vector3D v( mcp->getVertex()[0], mcp->getVertex()[1], mcp->getVertex()[2] );
      gear::Vector3D e( mcp->getEndpoint()[0], mcp->getEndpoint()[1], mcp->getEndpoint()[2] );
      gear::Vector3D p( mcp->getMomentum()[0], mcp->getMomentum()[1], mcp->getMomentum()[2] );
    
      trueE.push_back(mcp->getEnergy());
      trueP.push_back(p.r());
      truePt.push_back(p.trans());
      trueTheta.push_back(p.theta());
      truePhi.push_back(p.phi());
      truePDG.push_back(mcp->getPDG());
      if (mcp->getParents()[0]) {
        trueMother.push_back(mcp->getParents()[0]->getPDG());
      }  
      else {
        trueMother.push_back(-1);
      }
      
      int nseen = 0;
      int nseenphoton = 0;
      double weightsum = 0;
   
      for (unsigned int idaughter = 0; idaughter < daughters.size(); idaughter++) {

        MCParticle* mcd = daughters[idaughter];

        streamlog_out(DEBUG) << " Daughter MCParticle " << idaughter << " pdg ID = " << mcd->getPDG() 
                             << ", genstat = " << mcd->getGeneratorStatus() 
                             << ", E = " << mcd->getEnergy() << std::endl;
        streamlog_out(DEBUG) << " Momentum = " << mcd->getMomentum()[0] << " " << mcd->getMomentum()[1] << " " << mcd->getMomentum()[2] << std::endl;
        streamlog_out(DEBUG) << " Vertex =   " << mcd->getVertex()[0] << " " << mcd->getVertex()[1] << " " << mcd->getVertex()[2] << std::endl;
        streamlog_out(DEBUG) << " Endpoint = " << mcd->getEndpoint()[0] << " " << mcd->getEndpoint()[1] << " " << mcd->getEndpoint()[2] << std::endl; 

        streamlog_out(DEBUG) << " get reco particle for daughter " << idaughter << std::endl;  
        const EVENT::LCObjectVec& recovec = mc2recoNav.getRelatedToObjects(mcd);
        streamlog_out(DEBUG) << " recovec has length " << recovec.size() << std::endl;  
        const EVENT::FloatVec& recoweightvec = mc2recoNav.getRelatedToWeights(mcd);
        //streamlog_out(DEBUG) << " recoweightvec has length " << recoweightvec.size() << std::endl;  
        double maxcaloweight = 0;
        int imaxcaloweight = -1;
        // reconstructed at all?
        if (recovec.size() > 0)  nseen++; 
        for (unsigned int irel = 0; irel < recovec.size(); irel++) {
          ReconstructedParticle* dummy =  (ReconstructedParticle*) recovec.at(irel); 
          streamlog_out(DEBUG) << " irel " << irel << ", recoweight = " << int(recoweightvec.at(irel)) 
                               << ", type = " << dummy->getType()
                               << ", recoweight%10000 (track) = " << int(recoweightvec.at(irel))%10000 
                               << ", recoweight/10000 (calo) = " << int(recoweightvec.at(irel))/10000 << std::endl;  
          double caloweight = double((int(recoweightvec.at(irel))/10000)/1000.);
          if (dummy->getType() == 22 && caloweight > maxcaloweight) {
            imaxcaloweight = irel;
            maxcaloweight = caloweight;
          }
        }
      
       
        // do we really have a photon?
        if (maxcaloweight > 0.1) {
          streamlog_out(MESSAGE) << " found reco photon for mcd at imaxcaloweight = " << imaxcaloweight << " with weight = " << maxcaloweight << std::endl ;
          nseenphoton++;
          weightsum += maxcaloweight;
        }  
        
      }  // loop over daughters
      
      isSeen.push_back(nseen);
      isSeenAsPhotons.push_back(nseenphoton);
      weightToPhotons.push_back(weightsum);
    
    }  // if keep
    
  }  // loop over MCPs

  streamlog_out(DEBUG) << " nMCGammaGammaParticles Found = " << nMCPi0 
                       << " ( " << nMCPi0-nMCEta-nMCEtaPrime << " " << nMCEta << " " << nMCEtaPrime << " )" << " Esum = " << ESum_MCP << std::endl;
*/
      
  streamlog_out(DEBUG) << " reco iterator and navigator " << std::endl;
  LCIterator<ReconstructedParticle> ggpIt( evt, _gammaGammaParticleCollectionName ) ;
  streamlog_out(DEBUG) << " got ggpIt with length of " << ggpIt.size() << std::endl;
  LCRelationNavigator rec2mcNav(evt->getCollection( _recoToTrue )); 
  streamlog_out(DEBUG) << " got rec2mcNav from " << rec2mcNav.getFromType() << " to " << rec2mcNav.getToType() << std::endl;
  
  //----------------------------------------------------------------------------------------------------------------------------
  // loop over GammaGammaParticles
  //----------------------------------------------------------------------------------------------------------------------------
  
  double ESum_Reco = 0.0;
  double ESum_Reco_Correct = 0.0;
  double ESum_Reco_MC_Correct = 0.0;
  double ESum_Reco_Wrong = 0.0;
  double ESum_Reco_MC_All = 0.0;
  double ESum_Reco_Meas = 0.0;

  int nCorrectPi0 = 0;
  int nCorrectEta = 0;
  int nCorrectEtaPrime = 0;

  int nRecoPi0 = 0;
  int nRecoEta = 0;
  int nRecoEtaPrime = 0;

  while( ReconstructedParticle* ggp = ggpIt.next()  ) {

    double Egg_Meas = 0.0;
    
    if (!ggp) continue;

    streamlog_out(DEBUG) << " GammaGammaParticle " << nRecoPi0 << " type = " << ggp->getType() << " E = " << ggp->getEnergy() 
                         << " GoodnessOfPid " << ggp->getGoodnessOfPID() << std::endl; 
    streamlog_out(DEBUG) << " Momentum = " << ggp->getMomentum()[0] << " " << ggp->getMomentum()[1] << " " << ggp->getMomentum()[2] << std::endl; 
    
    if(ggp->getType()==111)nRecoPi0++;
    if(ggp->getType()==221)nRecoEta++;
    if(ggp->getType()==331)nRecoEtaPrime++;

    ESum_Reco += ggp->getEnergy();   // Fitted energy sum (regardless of whether the fit is correct
            
    // get photons
    const EVENT::ReconstructedParticleVec& gammas = ggp->getParticles ();
    int istrue = 0;
    int pdg = 0;
    int ntruephoton = 0;
    int ntruemeson = 0;
    double sumTrueE = 0;
    double sumTrueP[3] = {0, 0, 0};
    double sumE = 0;
    double sumP[3] = {0, 0, 0};
    double sumWeight = 0;
    MCParticle* mcgps[2];

    streamlog_out(DEBUG) << " Consitituent gammas.size() " << gammas.size() << std::endl;
   
    for (unsigned int igamma = 0; igamma < gammas.size(); igamma++) {

         streamlog_out(DEBUG) << " Gamma RP constituent of GGP " << igamma << " type = " << gammas[igamma]->getType() 
                              << " Emeas = " << gammas[igamma]->getEnergy() << std::endl; 
         streamlog_out(DEBUG) << " Momentum = " << gammas[igamma]->getMomentum()[0] << " " 
                                                << gammas[igamma]->getMomentum()[1] << " " 
                                                << gammas[igamma]->getMomentum()[2] << std::endl; 

      sumE += gammas[igamma]->getEnergy();

      Egg_Meas += gammas[igamma]->getEnergy();     // Measured energy (this is zeroed for each ggp ...)

//      ESum_Reco_Meas += gammas[igamma]->getEnergy();    // Measured Energy Sum of the Photons that are fitted to GammaGammaParticles
      for (int i = 0; i < 3; i++) sumP[i] += gammas[igamma]->getMomentum()[i];        
      streamlog_out(DEBUG) << " get mc particle for gamma " << igamma << std::endl;  
      const EVENT::LCObjectVec& truevec = rec2mcNav.getRelatedToObjects(gammas[igamma]);
      streamlog_out(DEBUG) << " truthvec has length " << truevec.size() << std::endl;  
      const EVENT::FloatVec& truthweightvec = rec2mcNav.getRelatedToWeights(gammas[igamma]);
      //streamlog_out(DEBUG) << " truthweightvec has length " << truthweightvec.size() << std::endl;  
      double maxcaloweight = 0;
      int imaxcaloweight = -1;
      double maxtrckweight = 0;
      int imaxtrckweight = -1;
      double maxweight = 0;
      int imaxweight = -1;
      MCParticle* mcg = 0;
      
      bool isConversion = false;
      // catch conversions
      if (truevec.size() == 2) {
        MCParticle* el1 =  (MCParticle*) truevec.at(0); 
        MCParticle* el2 =  (MCParticle*) truevec.at(1); 
        if ((el1->getPDG() == 11 && el2->getPDG() == -11) || (el1->getPDG() == -11 && el2->getPDG() == 11)) {
          MCParticle* mo1 = el1->getParents()[0];
          MCParticle* mo2 = el2->getParents()[0];
          if (mo1->getPDG() == 22 && mo1 == mo2) {
            streamlog_out(MESSAGE) << " found conversion - using mother photon: "  
                                   << " true PDGs = " << el1->getPDG() << ", " << el2->getPDG()
                                   << ", true parent PDG = " << mo1->getPDG()  << std::endl;  
            isConversion = true;
            maxweight = double( (int(truthweightvec.at(0))/10000 + int(truthweightvec.at(1))/10000)/1000. );
            // conversion mother is true photon!
            mcg = mo1;  
          }
        }
      }
        
      if (!isConversion) {  
        for (unsigned int irel = 0; irel < truevec.size(); irel++) {
          MCParticle* dummy =  (MCParticle*) truevec.at(irel);

// calculate cosine of the angle between the Photon PFO and this related MCParticle's momentum 
          double cosg = 0.0;
          double momentumMCSquared = 0.0;
          for(int i=0;i<3;++i){
              cosg += dummy->getMomentum()[i]*gammas[igamma]->getMomentum()[i];
              momentumMCSquared += pow(dummy->getMomentum()[i],2);
          }
          cosg = cosg/(sqrt(momentumMCSquared)*gammas[igamma]->getEnergy()); 
 
          streamlog_out(DEBUG) << " irel " << irel << " E = " << dummy->getEnergy() << ", truew = " << int(truthweightvec.at(irel)) 
                               << ", true PDG = " << dummy->getPDG() << " cos(RP - MC) " << cosg
                               << ", true par. PDG = " << dummy->getParents()[0]->getPDG()
                               << ", recow%10000 (track) = " << int(truthweightvec.at(irel))%10000 
                               << ", recow/10000 (calo) = " << int(truthweightvec.at(irel))/10000 << std::endl;  
          double truthcaloweight = double( (int(truthweightvec.at(irel))/10000)/1000. );
          if (dummy->getPDG() == 22 && truthcaloweight > maxcaloweight) {
            imaxcaloweight = irel;
            maxcaloweight = truthcaloweight;
          }
          double truthtrckweight = double((int(truthweightvec.at(irel))%10000)/1000.);
          if (dummy->getPDG() == 22 && truthtrckweight > maxtrckweight) {
            imaxtrckweight = irel;
            maxtrckweight = truthtrckweight;
          }
        } // irel loop
        streamlog_out(DEBUG) << " found true photon at imaxcaloweight = " << imaxcaloweight << " with weight = " << maxcaloweight << std::endl ;
        streamlog_out(DEBUG) << " found true photon at imaxtrckweight = " << imaxtrckweight << " with weight = " << maxtrckweight << std::endl ;
        
        maxweight = maxcaloweight;
        imaxweight = imaxcaloweight;
        if (maxtrckweight > maxcaloweight) {
          maxweight = maxtrckweight;
          imaxweight = imaxtrckweight;
        }  
        
        // skip if photon has too small weight
        if (maxweight < 0.1) {
          streamlog_out(DEBUG) << " weight of photon = " << maxweight << " is too small - don't count as correct" << std::endl ;
          continue;
        }
        
        streamlog_out(DEBUG) << " found true photon at imaxweight = " << imaxweight << " with weight = " << maxweight << std::endl ;
        mcg =  (MCParticle*) truevec.at(imaxweight); 
      } // !isConversion if clause
      
      
      if (!mcg) continue;
      
      // found a correct photon or a conversion
      ntruephoton++; 

      sumWeight += maxweight;  
      sumTrueE += mcg->getEnergy();
      for (int i = 0; i < 3; i++) sumTrueP[i] += mcg->getMomentum()[i]; 

      streamlog_out(DEBUG) << " Gamma MCParticle of GGP " << igamma  
                           << " Etrue = " << mcg->getEnergy() << std::endl;

      ESum_Reco_MC_All += mcg->getEnergy();

      // get parent   
      mcgps[igamma] = mcg->getParents()[0];
      if (!mcgps[igamma]) continue;
      if (mcgps[igamma]->getPDG() == typeRequired ) {
        // found photon from desired meson
        ntruemeson++; 
        streamlog_out(MESSAGE) << "(correct) PDG of photon's parent is = " << mcgps[igamma]->getPDG() << std::endl ;
        pdg = mcgps[igamma]->getPDG();
      }  
      else {
        streamlog_out(MESSAGE) << "(wrong) PDG of photon's parent is = " << mcgps[igamma]->getPDG() << std::endl ;
      } 
    }  // Loop over constituent reconstructed particles (gammas) in the GammaGammaCandidate

    streamlog_out(DEBUG) << " ntruephoton = " << ntruephoton << std::endl;

    // check if really from same mother and has correct pdgid  (Note that mcgps may not be valid ... ) 

    if(ntruephoton==2){
       if (mcgps[0] && mcgps[1] && mcgps[0] == mcgps[1] && mcgps[0]->getPDG() == typeRequired ) {
         double rVertexSquared = 0.0;
         for(int i=0; i<=2; i++){
             rVertexSquared += pow(mcgps[0]->getVertex()[i],2);
         }
         if(rVertexSquared < 10000.0)istrue = 1;
         if(istrue){
            pdg = mcgps[0]->getPDG();
            streamlog_out(MESSAGE) << "(fully correct) PDG of both photons parent is = " << pdg << " Etrue = " << mcgps[0]->getEnergy() 
                                   << " VTX squared = " << rVertexSquared << std::endl ;
            ESum_Reco_Correct += ggp->getEnergy();
            ESum_Reco_MC_Correct += mcgps[0]->getEnergy();
            ESum_Reco_Meas += Egg_Meas; 
            if(pdg == 111)nCorrectPi0++;
            if(pdg == 221)nCorrectEta++;
            if(pdg == 331)nCorrectEtaPrime++;

// Add this RP to the output collection
//         ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
            recparcol->addElement(ggp);
         }
         else{
            streamlog_out(MESSAGE) << "(Non-prompt associated) PDG of both photons parent is = " << pdg << " Etrue = " << mcgps[0]->getEnergy() 
                                << " VTX squared = " << rVertexSquared << std::endl ;         
         }
      }
      else{
         ESum_Reco_Wrong += ggp->getEnergy();
      }
    }           
  }  // loop over GammaGammaParticles
  
  if(typeRequired == 111)streamlog_out(MESSAGE) << "                SUMMARY 111 " << nCorrectPi0 << " " 
                                                << ESum_Reco_Correct << " " << ESum_Reco_MC_Correct << " " << ESum_Reco_Meas << endl;
  if(typeRequired == 221)streamlog_out(MESSAGE) << "                SUMMARY 221 " << nCorrectEta << " " 
                                                << ESum_Reco_Correct << " " << ESum_Reco_MC_Correct << " " << ESum_Reco_Meas << endl;
  if(typeRequired == 331)streamlog_out(MESSAGE) << "           SUMMARY 331 " << nCorrectEtaPrime << " " 
                                                << ESum_Reco_Correct << " " << ESum_Reco_MC_Correct << " " << ESum_Reco_Meas << endl;

/*
  streamlog_out(DEBUG) << " nRecoGammaGammaParticles Found = " << nRecoPi0 << 
                       " ( " << nRecoPi0-nRecoEta-nRecoEtaPrime << " " << nRecoEta << " " << nRecoEtaPrime << " ) " << 
                       " ESum Reco = " << ESum_Reco << std::endl;
  streamlog_out(DEBUG) << " SUMMARY = " << nMCPi0 << " ESum MCP = " << ESum_MCP << " " << nRecoPi0 
                       << " ESum Reco = " << ESum_Reco << "  Esum_Reco_Correct = " << ESum_Reco_Correct 
                       << " ESum Reco Wrong = " << ESum_Reco_Wrong 
                       << " ESum_Reco_MC_Correct " << ESum_Reco_MC_Correct
                       << " ESum_Reco_Meas " << ESum_Reco_Meas << 
                       " ESum_Reco_MC_All " << ESum_Reco_MC_All << std::endl;

  streamlog_out(DEBUG) << " SUMMARY2 = " << std::setw(2) << nMCPi0 << " " << std::setw(10) << ESum_MCP << " " 
                                         << std::setw(2) << nMCPi0-nMCEta-nMCEtaPrime << " " 
                                         << std::setw(2) << nMCEta << " " << std::setw(2) << nMCEtaPrime << " "
                                         << std::setw(2) << nRecoPi0 << " " << std::setw(10) << ESum_Reco << " "
                                         << std::setw(2) << nRecoPi0-nRecoEta-nRecoEtaPrime << " " 
                                         << std::setw(2) << nRecoEta << " " << std::setw(2) << nRecoEtaPrime << " "
                                         << std::setw(2) << nCorrectPi0 << " " 
                                         << std::setw(2) << nCorrectEta << " " << std::setw(2) << nCorrectEtaPrime << " " 
                                         << std::setw(10) << ESum_Reco_Correct << " " << std::setw(10) << ESum_Reco_Wrong << " " 
                                         << std::setw(10) << ESum_Reco_MC_Correct << " "
                                         << std::setw(10) << ESum_Reco_Meas << " " 
                                         << std::setw(10) << ESum_Reco_MC_All << std::endl;
*/

  nEvt++;

  // Add new collection to event
  evt->addCollection( recparcol , _outputParticleCollectionName.c_str() );
  
  cout << "======================================== event " << nEvt << std::endl ;
  
}


void GammaGammaCandidateTruthFilter::check( LCEvent * evt ) { 

}


void GammaGammaCandidateTruthFilter::end(){ 

}



