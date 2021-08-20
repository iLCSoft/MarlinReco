#include <stdlib.h>
#include "TaJetClustering.h"
#include <iostream>
#include <math.h>
#include <algorithm>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "TVector3.h"


// 090128 large change: angle based -> invmass based

using namespace lcio ;
using namespace marlin ;
using namespace std;

TaJetClustering aTaJetClustering ;

TaJetClustering::TaJetClustering() : Processor("TaJetClustering")
{

  // modify processor description
  _description = "Tau Jet Finding" ;
  
  // collections
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "PFOCollection",
			   "PandoraPFA PFO collection",
			   _pfoCollectionName,
			   std::string("PandoraPFOs") );
  
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "OutputTauCollection",
			    "Tau output collection",
			    _tauCollectionName,
			    std::string("TaJets") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "RemainPFOCollection",
			    "Remained PFO collection",
			    _remainCollectionName,
			    std::string("RemainPFOs") );
  
  // steering parameters
  //registerProcessorParameter( "MaxTaus", "Number of maximum taus to be accepted", _maxTaus,100);
  //registerProcessorParameter( "ReturnCutTaus", "Return cut taus to remain PFO collection", _returnCutTaus,1);
  
  // clustering params
  registerProcessorParameter( "TauMass",
			      "Tau mass for tau clustering [GeV]",
			      _tauMass,
			      2.0 );
  registerProcessorParameter( "TauCosAngle",
			      "Allowed cosine angle to be clustered",
			      _tauCosAngle,
			      0.98 );
  
  // selection switch
  registerProcessorParameter( "NoSelection",
			      "Neither primary nor cone cut at all if true. Only clustering.",
			      _noSelection,
			      0 );
  // primary cuts
  registerProcessorParameter( "MinimumTrackEnergy",
			      "Minimum track energy to be accepted as taus",
			      _minimumTrackEnergy,
			      2.0 );
  registerProcessorParameter( "MinimumJetEnergy",
			      "Minimum jet energy to be accepted as taus",
			      _minimumJetEnergy,
			      3.0 );
  registerProcessorParameter( "MinimumTrackEnergyAssoc",
			      "Minimum track energy to be counted",
			      _minimumTrackEnergyAssoc,
			      2.0 );
  registerProcessorParameter( "MaximumNeutralEnergyIn3Prong",
			      "Reject 3-prong taus if including neutral energy more than this number in GeV",
			      _maximumNeutralEnergyInThreeProng,
			      1.0 );
  registerProcessorParameter( "AcceptFlexibleLowEnergyTrack",
			      "Low energy tracks can be accepted as either charged and neutral if true",
			      _acceptFlexibleLowEnergyTrack,
			      1 );
  
  // cone cuts
  registerProcessorParameter( "ConeMinCosAngle",
			      "Minimum cosine angle for cone",
			      _coneMinCosAngle,
			      0.90 );
  registerProcessorParameter( "ConeMaxCosAngle",
			      "Maximum cosine angle for cone",
			      _coneMaxCosAngle,
			      1.0 );
  registerProcessorParameter( "ConeMaxEnergyFrac",
			      "Energy fraction of cone compared to central",
			      _coneMaxEnergyFrac,
			      0.1 );

  // following cuts are found not useful now, so commented out for simplicity
  /*  
  // impact parameter cuts
  registerProcessorParameter( "IPCutsMinimumTrackEnergy",
			      "Minimum track energy for impact parameter acception",
			      _ipCutsMinimumTrackEnergy,
			      1000.0 );
  registerProcessorParameter( "IPCutsSigmaOneProng",
			      "Minimum number of sigma for one prong tau tracks",
			      _ipCutsSigmaOneProng,
			      3. );
  registerProcessorParameter( "IPCutsSigmaThreeProngNoNeutrals",
			      "Minimum number of sigma for three prong tau without neutrals",
			      _ipCutsSigmaThreeProngNoNeutrals,
			      5. );
  registerProcessorParameter( "IPCutsSigmaThreeProngWithNeutrals",
			      "Minimum number of sigma for three prong tau with neutrals, used with cone cuts (so looser)",
			      _ipCutsSigmaThreeProngWithNeutrals,
			      3. );
  registerProcessorParameter( "IPMaxAbsolute",
			      "Maximum absolute d0/z0 for tau",
			      _ipMaxAbsolute,
			      1.0 );
  
  // lepton ID cuts
  registerProcessorParameter( "LeptonCutsMinimumTrackEnergy",
			      "Minimum track energy for lepton ID acception",
			      _leptonCutsMinimumTrackEnergy,
			      1000.0 );
  registerProcessorParameter( "MuMaxFracEcal",
			      "Max ECAL fraction for muon ID",
			      _muMaxFracEcal,
			      0.5 );
  registerProcessorParameter( "MuMaxCalByTrack",
			      "Max CAL/track energy ratio for muon ID",
			      _muMaxCalByTrack,
			      0.5 );
  registerProcessorParameter( "EMinFracEcal",
			      "Min ECAL fraction for electron ID",
			      _eMinFracEcal,
			      0.97 );
  registerProcessorParameter( "EMinCalByTrack",
			      "Min CAL/track energy ratio for electron ID",
			      _eMinCalByTrack,
			      0.9 );
  registerProcessorParameter( "EMaxCalByTrack",
			      "Max CAL/track energy ratio for electron ID",
			      _eMaxCalByTrack,
			      1.1 );
  */
}


void TaJetClustering::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void TaJetClustering::processRunHeader( LCRunHeader*  /*run*/ ) { 

  _nRun++ ;
} 

bool compareEnergy( ReconstructedParticle * const &p1, ReconstructedParticle * const &p2 )
{
  return p1->getEnergy() > p2->getEnergy();
}

void TaJetClustering::processEvent( LCEvent * evt ) { 

  // obtain Reconstructed Particles
  LCCollection* col = evt->getCollection(_pfoCollectionName);
  int nPart = col->getNumberOfElements();
  
  // output collection
  LCCollectionVec* pTaJetsCol= new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  LCCollectionVec* pRemainCol= new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
  pRemainCol->setSubset(true);
  
  vector< ReconstructedParticle * > parts;
  vector< ReconstructedParticle * > tajets;
  
  for( int i = 0; i < nPart; i++ ) {
    parts.push_back( dynamic_cast< ReconstructedParticle* >( col->getElementAt(i) ) );
  }
  
  sort( parts.begin(), parts.end(), compareEnergy );
  /*  
  cout << "PFOs in event #" << _nEvt << endl;
  for(int i = 0; i < nPart; i++ ){
    cout << parts[i]->getEnergy() << (parts[i]->getCharge() ? "C" : "N") << endl;
  }
  */
  // Jet finding: around charged particle
  // because we do not care neutral jets for tau finding
  //cout << "Tau reconstruction start." << endl;
  for( int i = 0; i < nPart; i++ )
    {
      ReconstructedParticle *pCur = parts[i];
      if( pCur == NULL ) continue;       // already associated particle skipped.
      if( !pCur->getCharge() ) continue; // neutral particle skipped.
      if( pCur->getEnergy() < _minimumTrackEnergy ) break;
      
      //cout << "Targetting to " << i << "th charged particle: energy = " << pCur->getEnergy() << endl;
      
      // create output jet
      vector< ReconstructedParticle * > jetTracks, jetNeutrals;
      //pJet->addParticle(pCur);
      jetTracks.push_back( pCur );
      
      double jetEnergy = pCur->getEnergy();
      TVector3 jetMomentum(pCur->getMomentum());
      double jetCharge = pCur->getCharge();
      int ntracks = 1;
      double eneutral = 0.;
      /*
      cout << "Direction: ( " << pCur->getMomentum()[0] << ", " << pCur->getMomentum()[1] << ", "
	   << pCur->getMomentum()[2] << " )" << endl;
      */
      int nPartAssoc = 1, nPartAssocPriv = 0;
      
      while( nPartAssoc != nPartAssocPriv ){
	for( int j = 0; j < nPart; j++ )
	  {
	    ReconstructedParticle *pCur2 = parts[j];
	    
	    // remove particles which is already included in the current jet 
	    if( pCur2 == NULL ) continue;	 // already associated particle skipped.
	    if( pCur2->getCharge() && find( jetTracks.begin(), jetTracks.end(), pCur2 ) != jetTracks.end() ) continue;
	    if( !pCur2->getCharge() && find( jetNeutrals.begin(), jetNeutrals.end(), pCur2 ) != jetNeutrals.end() ) continue;
	    
	    TVector3 vec2( pCur2->getMomentum() );
	    
	    double cosangle = vec2.Unit().Dot( jetMomentum.Unit() );
	    double invmass2 = pow( pCur2->getEnergy() + jetEnergy, 2 ) - ( jetMomentum + vec2 ).Mag2();
	    /*
	    cout << "Candidate particle found: energy = " << pCur2->getEnergy() << ", " ;
	    cout << "Direction: ( " << pCur2->getMomentum()[0] << ", " << pCur2->getMomentum()[1] << ", "
		 << pCur2->getMomentum()[2] << " ), charge = " << pCur2->getCharge()
		 << ", cosangle = " << cosangle << ", invmass = " << sqrt(invmass2) << endl;
	    */
	    if( cosangle < _tauCosAngle ) continue;	 // acceptance angle cut; skipped.
	    if( invmass2 > _tauMass*_tauMass ) continue;	// invmass cut; skipped.
	    
	    //cout << "Candidate particle associated" << endl;
	    
	    jetEnergy += pCur2->getEnergy();
	    jetMomentum += vec2;
	    //pJet->addParticle(pCur2);
	    if( pCur2->getCharge() != 0. ){
	      // charge of lower track energy is ignored
	      if( pCur2->getEnergy() > _minimumTrackEnergyAssoc || _acceptFlexibleLowEnergyTrack ){
		ntracks ++;
		jetCharge += pCur2->getCharge();
		jetTracks.push_back( pCur2 );
	      }
	      else{
		// have to think about addition to eneutral
		jetNeutrals.push_back( pCur2 );
	      }
	    }
	    else{
	      eneutral += pCur2->getEnergy();
	      jetNeutrals.push_back( pCur2 );
	    }
	  }
	nPartAssocPriv = nPartAssoc;
	nPartAssoc = jetTracks.size() + jetNeutrals.size();
      }
      
      if(!_noSelection){
	//cout << "primary selection" << endl;
	// primary selection
	if( jetEnergy < _minimumJetEnergy ){ continue; }
      
	// #	 tracks, charge selection
	if( _acceptFlexibleLowEnergyTrack ){
	  while( ( ntracks != 1 && ntracks != 3 ) || abs( int( jetCharge ) ) != 1 ){
	    ReconstructedParticle *p = jetTracks[ jetTracks.size()-1 ];
	    if( p->getEnergy() > _minimumTrackEnergyAssoc ) break;
	    
	    //	 removal of the track
	    jetCharge -= p->getCharge();
	    ntracks --;
	    
	    jetTracks.pop_back();
	    
	    //	 add to neutral instead
	    jetNeutrals.push_back(p);
	    
	    //	cout << "Flexible low energy tracks: track treated as neutral: e = " << p->getEnergy() << ", c = " << p->getCharge();
	    //cout	 << ", ntracks now " << ntracks << ", jc = " << jetCharge << endl;
	    
	  }
	}
	if( ( ntracks != 1 && ntracks != 3 ) || abs( int( jetCharge ) ) != 1 ){ continue; }
      }
      // set RP
      ReconstructedParticleImpl *pJet = new ReconstructedParticleImpl;
	
      pJet->setEnergy( jetEnergy );
      double mom[3];
      mom[0] = jetMomentum.x();
      mom[1] = jetMomentum.y();
      mom[2] = jetMomentum.z();
      pJet->setMomentum( mom );
      pJet->setCharge( jetCharge );
      
      // put into reconstructed jet: awkward method...
      for( unsigned int n = 0; n < jetTracks.size(); n++ )
	pJet->addParticle( jetTracks[n] );
      for( unsigned int n = 0; n < jetNeutrals.size(); n++ )
	pJet->addParticle( jetNeutrals[n] );

      if(!_noSelection){
      
	// 3prong + neutral rejection
	if( ntracks == 3 && eneutral > _maximumNeutralEnergyInThreeProng ){ delete pJet; continue; }
 
	//cout << "cone selection" << endl;
	// calculate cone energy
	double econe = 0.;
	for( int j = 0; j < nPart; j++ ) {
	  ReconstructedParticle *p = dynamic_cast< ReconstructedParticle* >( col->getElementAt(j) );
	
	  // remove self
	  if( find( pJet->getParticles().begin(), pJet->getParticles().end(), p ) != pJet->getParticles().end() ) continue;
	
	  TVector3 v( p->getMomentum() );
	  double cosangle = cos( v.Angle( jetMomentum ) );
	
	  if( cosangle > _coneMinCosAngle && cosangle < _coneMaxCosAngle ){
	    econe += p->getEnergy();
	    //cout << "particle added to cone; costheta = " << cosangle << ", e = " << p->getEnergy() << ", charge = " << p->getCharge() << endl;
	  }
	}
	//cout << "econe = " << econe << ", jet energy = " << jetEnergy << endl;
      
	if( econe > _coneMaxEnergyFrac * jetEnergy ){ delete pJet; continue; }
      }

      /*      
      if( ntracks == 3 && eneutral > 1. && econe > _coneMaxEnergyFrac * jetEnergy ){ delete pJet; continue; }
      
      if( (ntracks == 3 && eneutral > 1. ) || ( econe > _coneMaxEnergyFrac * jetEnergy ) ){
	//if( ntracks == 3 && eneutral > 1. ) cout << "cone cut passed: move to IP cut." << endl;
	//else cout << "cone cut rejected: move to IP cut." << endl;
	
	// IP cut
	bool ipcutaccepted = false;
	for( unsigned int j = 0; j < pJet->getParticles().size(); j++){
	  ReconstructedParticle *p = pJet->getParticles()[j];
	  
	  if( p->getCharge() == 0 ) continue; // track only
	  if( p->getTracks().size() == 0 ) continue; // skip with no tracks
	  if( p->getEnergy() < _ipCutsMinimumTrackEnergy ) continue; // too small energy track
	  
	  double d0 = p->getTracks()[0]->getD0();
	  double z0 = p->getTracks()[0]->getZ0();
	  
	  double d0err = sqrt(p->getTracks()[0]->getCovMatrix()[0]);
	  double z0err = sqrt(p->getTracks()[0]->getCovMatrix()[9]);
	  
	  //cout << "Track: d0 " << d0 << " d0err " << d0err << " z0 " << z0 << " z0err " << z0err << endl;
	  
	  // 3-prong, neutral
	  if( ntracks == 3 && eneutral > 1. &&
	     ( ( fabs( d0/d0err ) > _ipCutsSigmaThreeProngWithNeutrals || fabs( z0/z0err ) > _ipCutsSigmaThreeProngWithNeutrals ) &&
	     ( d0 < _ipMaxAbsolute && z0 < _ipMaxAbsolute ) ) ) ipcutaccepted = true;
	  
	  // 3-prong, no neutral
	  if( ntracks == 3 && eneutral < 1. &&
	     ( ( fabs( d0/d0err ) > _ipCutsSigmaThreeProngNoNeutrals || fabs( z0/z0err ) > _ipCutsSigmaThreeProngNoNeutrals ) &&
	     ( d0 < _ipMaxAbsolute && z0 < _ipMaxAbsolute ) ) ) ipcutaccepted = true;
	  
	  // 1-prong
	  if( ntracks == 1 &&
	     ( ( fabs( d0/d0err ) > _ipCutsSigmaOneProng || fabs( z0/z0err ) > _ipCutsSigmaOneProng ) && 
	     ( d0 < _ipMaxAbsolute && z0 < _ipMaxAbsolute ) ) ) ipcutaccepted = true;
	}
	
	if( ntracks == 1 && eneutral < 1. && ipcutaccepted == false ){
	  //cout << "IP cut rejected: move to lepton ID cut." << endl;
	  
	  bool leptoncutaccepted = false;
	  for( unsigned int j = 0; j < pJet->getParticles().size(); j++ ){
	    ReconstructedParticle *p = pJet->getParticles()[j];
	    
	    if( p->getCharge() == 0 ) continue; // track only
	    if( p->getTracks().size() == 0 ) continue; // skip with no tracks
	    if( p->getClusters().size() == 0 ) continue; // skip with no neutrals
	    if( p->getEnergy() < _leptonCutsMinimumTrackEnergy ) continue; // too small energy track
	    
	    double ecal = p->getClusters()[0]->getSubdetectorEnergies()[0];
	    double hcal = p->getClusters()[0]->getSubdetectorEnergies()[1];
	    double tre = p->getEnergy();
	    
	    double fracecal = ecal / ( ecal + hcal );
	    double calbytrack = ( ecal + hcal ) / tre;
	    
	    //cout << "track for lepton ID: e = " << p->getEnergy() << " fracecal = " << fracecal << " calbytrack = " << calbytrack << endl;
	    
	    // muon ID
	    if( fracecal < _muMaxFracEcal && calbytrack < _muMaxCalByTrack ) leptoncutaccepted = true;
	    // electron ID
	    if( fracecal > _eMinFracEcal && calbytrack > _eMinCalByTrack && calbytrack < _eMaxCalByTrack ) leptoncutaccepted = true;
	  }
	  if( leptoncutaccepted == false ){
	    //cout << "lepton cut rejected: the jet is finally rejected." << endl;
	    delete pJet; continue;
	  }
	}
	else if( ipcutaccepted == false ){
	  //cout << "IP cut rejected: the jet is rejected." << endl;
	  delete pJet; continue;
	}
      }
      */
      
      //cout << "The jet is accepted for tau!";
      //cout << " e = " << jetEnergy << " : pz = " << jetMomentum.z() << " npart = " << pJet->getParticles().size() << endl;

      pTaJetsCol->addElement( pJet );
      // remove from parts[]
      for( int j = 0; j < nPart; j++ ) {
	ReconstructedParticle *p = parts[j];
	if(p == 0) continue;
	
	// remove self
	if( find( pJet->getParticles().begin(), pJet->getParticles().end(), p ) != pJet->getParticles().end() ) parts[j] = 0;
      }
      
    }
  //cout << "Tau reconstruction result: " << pTaJetsCol->getNumberOfElements() << " jets found." << endl;
  evt->addCollection( pTaJetsCol,_tauCollectionName );
  
  // make remaining collection
  for( int j = 0; j < nPart; j++ ) {
    ReconstructedParticle *p = parts[j];
    if(p == 0) continue;
    
    pRemainCol->addElement(p);
  }
  evt->addCollection( pRemainCol, _remainCollectionName );

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber()
		       << "   in run:  " << evt->getRunNumber()
		       << std::endl ;

  _nEvt ++ ;
}



void TaJetClustering::check( LCEvent *  /*evt*/ ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TaJetClustering::end(){ 
  
//   std::cout << "TaJetClustering::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

