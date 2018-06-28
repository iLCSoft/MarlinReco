#include "Adjusted_TrueJet_Parser.h"
#include <stdlib.h>
#include <math.h>
//#include <cmath>
#include <iostream>
#include <iomanip>


// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA


using namespace lcio ;
using namespace marlin ;


  Adjusted_TrueJet_Parser::Adjusted_TrueJet_Parser() {
   
    intvec=new IntVec()    ;
    mcpartvec=new MCParticleVec()    ;
    _COUNT_FSR=1;
}

Adjusted_TrueJet_Parser::~Adjusted_TrueJet_Parser() {
}



//*************************///
const double* Adjusted_TrueJet_Parser::p4seen(int ijet) { 
  p4[0]=Eseen(ijet);
  const double* mom = pseen(ijet);
   for (int kk=1 ; kk<=3 ; kk++ ) {
     p4[kk]=mom[kk-1];
   }
   return p4 ;
}


double Adjusted_TrueJet_Parser::Etrue(int ijet) {
  static FloatVec www ;
  LCObjectVec mcpvec = reltjmcp->getRelatedToObjects( jets->at(ijet) );
  www = reltjmcp->getRelatedToWeights( jets->at(ijet));
  double E=0.0;
  for ( unsigned kk=0 ; kk<mcpvec.size() ; kk++ ) {
    MCParticle* mcp  = dynamic_cast<MCParticle*>(mcpvec[kk]);
    if ( _COUNT_FSR ) {
      if (  abs(www[kk]) == 1.0) {
        E+=mcp->getEnergy();
      }
    } else {
      if (  www[kk] == 1.0) {
        E+=mcp->getEnergy();
      }
    }
  }
  return E;
}
double Adjusted_TrueJet_Parser::Mtrue(int ijet) {
  const double* p4 = p4true(ijet);
  double psqr=0;
  for (int kk=1 ; kk<=3 ; kk++ ) {
     psqr+=p4[kk]*p4[kk];
  }
  double M=sqrt(p4[0]*p4[0]-psqr);
  return M ;
}
const double* Adjusted_TrueJet_Parser::ptrue(int ijet) {
  static FloatVec www ;
  LCObjectVec mcpvec = reltjmcp->getRelatedToObjects( jets->at(ijet) );
  www =  reltjmcp->getRelatedToWeights( jets->at(ijet));
  //double* p=0; 
  p3[0]=0.; p3[1]=0.; p3[2]=0.; 
  for ( unsigned kk=0 ; kk<mcpvec.size() ; kk++ ) {
    MCParticle* mcp  = dynamic_cast<MCParticle*>(mcpvec[kk]);
    if ( _COUNT_FSR ) {
      if (  abs(www[kk]) == 1.0) {
        const double* mom = mcp->getMomentum();
        for (int kk=0 ; kk<3 ; kk++ ) {
          p3[kk]+=mom[kk];
        }
      }
    } else {
      if (  www[kk] == 1.0) {
        const double* mom = mcp->getMomentum();
        for (int kk=0 ; kk<3 ; kk++ ) {
          p3[kk]+=mom[kk];
        }
      }
    }
  }
  return p3;
}

const double* Adjusted_TrueJet_Parser::p4true(int ijet) {
  // double* p4=0; 
  const double* mom = ptrue(ijet);
  for (int kk=1 ; kk<=3 ; kk++ ) {
     p4[kk]=mom[kk-1];
  }
  p4[0]=Etrue(ijet);
  return p4 ;
}

const double* Adjusted_TrueJet_Parser::pquark(int ijet) {
  if (final_elementon(ijet) != NULL ) {
    return final_elementon(ijet)->getMomentum() ;
  } else {
    p3[0]=0. ;p3[1]=0. ;p3[2]=0. ;
    return p3 ;
  }
}

const double* Adjusted_TrueJet_Parser::p4quark(int ijet) {
  //  double* p4; 
  const double* mom = pquark(ijet);
  for (int kk=1 ; kk<=3 ; kk++ ) {
     p4[kk]=mom[kk-1];
  }
  p4[0]=Equark(ijet);
  return p4 ;
}


double Adjusted_TrueJet_Parser::Etrueseen(int ijet) {
  if (  reltrue_tj == 0 ) {
    LCCollection* rmclcol = NULL;
    try{
     rmclcol = evt->getCollection( get_recoMCTruthLink() );
    }
    catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << get_recoMCTruthLink()   << " collection not available" << std::endl;
        rmclcol = NULL;
    }
    reltrue_tj = new LCRelationNavigator( rmclcol );
  }  


  LCObjectVec mcpvec = reltjmcp->getRelatedToObjects( jets->at(ijet) );

  double E=0.0;
  for ( unsigned kk=0 ; kk<mcpvec.size() ; kk++ ) {
    MCParticle* mcp  = dynamic_cast<MCParticle*>(mcpvec[kk]);
    LCObjectVec recovec = reltrue_tj->getRelatedFromObjects( mcp);
   if ( recovec.size() > 0 ) { // if reconstructed
     E+=mcp->getEnergy();
    }
  }
  return E;
}
double Adjusted_TrueJet_Parser::Mtrueseen(int ijet) {
  const double* p4 = p4trueseen(ijet);
  double psqr=0;
  for (int kk=1 ; kk<=3 ; kk++ ) {
     psqr+=p4[kk]*p4[kk];
  }
  double M=sqrt(p4[0]*p4[0]-psqr);
  return M ;
}
const double* Adjusted_TrueJet_Parser::ptrueseen(int ijet) {
  if (  reltrue_tj == 0 ) {
    LCCollection* rmclcol = NULL;
     try{
      rmclcol = evt->getCollection( get_recoMCTruthLink() );
    }
    catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) << get_recoMCTruthLink()   << " collection not available" << std::endl;
        rmclcol = NULL;
    }
    reltrue_tj = new LCRelationNavigator( rmclcol );
  }  

  LCObjectVec mcpvec = reltjmcp->getRelatedToObjects( jets->at(ijet) );


  //double* p=0; 
  p3[0]=0 ; p3[1]=0 ; p3[2]=0 ; 
  for ( unsigned kk=0 ; kk<mcpvec.size() ; kk++ ) {
    MCParticle* mcp  = dynamic_cast<MCParticle*>(mcpvec[kk]);
    LCObjectVec recovec = reltrue_tj->getRelatedFromObjects( mcp);
    if ( recovec.size() > 0 ) { // if reconstructed

      const double* mom = mcp->getMomentum();
      for (int kk=0 ; kk<3 ; kk++ ) {
        p3[kk]+=mom[kk];
      }
    }
  }
  return p3;
}

const double* Adjusted_TrueJet_Parser::p4trueseen(int ijet) {
  //double* p4=0; 
  const double* mom = ptrueseen(ijet);
  for (int kk=1 ; kk<=3 ; kk++ ) {
     p4[kk]=mom[kk-1];
  }
  p4[0]=Etrueseen(ijet);
  return p4 ;
}


ReconstructedParticleVec* Adjusted_TrueJet_Parser::getJets(){
  //LCObjectVec* Adjusted_TrueJet_Parser::jets(){



    ReconstructedParticleVec* tjs = new ReconstructedParticleVec();
    //    LCObjectVec* tjs = new LCObjectVec();
    int ntj = tjcol->getNumberOfElements()  ;
    // std::cout << " n jets " << ntj << std::endl;
    
    for(int j=0; j< ntj ; j++){

      ReconstructedParticle* tj = dynamic_cast<ReconstructedParticle*>( tjcol->getElementAt( j ) ) ;
      tjs->push_back(tj);
      tj->ext<JetIndex>()=j;


    }
    return tjs;
}

 ReconstructedParticleVec* Adjusted_TrueJet_Parser::getFinalcn(){


    ReconstructedParticleVec* fcns = new  ReconstructedParticleVec();
    int nfcn = fcncol->getNumberOfElements()  ;
    //std::cout << " n fcn " << nfcn << std::endl;
    
    for(int j=0; j< nfcn ; j++){

      ReconstructedParticle* fcn = dynamic_cast<ReconstructedParticle*>( fcncol->getElementAt( j ) ) ;
      fcns->push_back(fcn);
      fcn->ext<FcnIndex>()=j;
    }
    return fcns;
}
 ReconstructedParticleVec* Adjusted_TrueJet_Parser::getInitialcn(){


    ReconstructedParticleVec* icns = new  ReconstructedParticleVec();
    int nicn = icncol->getNumberOfElements()  ;
    //std::cout << " n icn " << nicn << std::endl;
    
    for(int j=0; j< nicn ; j++){

      ReconstructedParticle* icn = dynamic_cast<ReconstructedParticle*>( icncol->getElementAt( j ) ) ;
      icns->push_back(icn);
      icn->ext<IcnIndex>()=j;
    }
    return icns;
}
const IntVec&  Adjusted_TrueJet_Parser::final_siblings( int ijet ) {
  intvec->clear();
  IntVec* sibl=intvec;
  LCObjectVec fcnvec = relfcn->getRelatedToObjects( jets->at(ijet) );
  int nsibl=0;
  for ( unsigned kk=0 ; kk<fcnvec.size() ; kk++ ) {
    ReconstructedParticleVec jetvec= dynamic_cast<ReconstructedParticle*>(fcnvec[kk])->getParticles();
    for ( unsigned jj=0 ; jj<jetvec.size() ; jj++ ) {
      if ( jetvec[jj] != jets->at(ijet) ) {
        // sibling
        int jjet=jetvec[jj]->ext<JetIndex>();
        sibl->push_back(jjet);
        nsibl++;
      }
    }
  }
  return *sibl;
  //

} 
const IntVec& Adjusted_TrueJet_Parser::initial_siblings( int ijet ){
  intvec->clear();
  IntVec* sibl=intvec;
  LCObjectVec icnvec = relicn->getRelatedToObjects( jets->at(ijet) );
  int nsibl=0;
  for ( unsigned kk=0 ; kk<icnvec.size() ; kk++ ) {
    ReconstructedParticleVec jetvec= dynamic_cast<ReconstructedParticle*>(icnvec[kk])->getParticles();
    for ( unsigned jj=0 ; jj<jetvec.size() ; jj++ ) {
      if ( jetvec[jj] != jets->at(ijet) ) {
        int jjet=jetvec[jj]->ext<JetIndex>();
        sibl->push_back(jjet);
        nsibl++;
      }
    }
  }
  return *sibl;
  //

} 
int Adjusted_TrueJet_Parser::final_cn( int ijet ) {
   LCObjectVec fcnvec = relfcn->getRelatedToObjects( jets->at(ijet) );
   int fcn ;
   if (fcnvec.size() > 0 ) {
     fcn=fcnvec[0]->ext<IcnIndex>();
   } else {
     fcn = -1 ;
   }
   return fcn;
} 
int Adjusted_TrueJet_Parser::initial_cn( int ijet ) {
   LCObjectVec icnvec = relicn->getRelatedToObjects( jets->at(ijet) );
   int icn ;
   if (icnvec.size() > 0 ) {
      icn=icnvec[0]->ext<IcnIndex>();
   } else {
     icn = -1 ;
   }

   return icn;
} 

const IntVec&  Adjusted_TrueJet_Parser::jets_of_final_cn( int ifcn ) {
  intvec->clear();
  IntVec* jets_of_fcn=intvec;
  // way to find: jet-to-icn link icn->reco : jets, find notthis.
  //  index<->jet : LCExtension
  LCObjectVec jetvec = relfcn->getRelatedFromObjects(  finalcns->at(ifcn) );
    // ReconstructedParticleVec jetvec= dynamic_cast<ReconstructedParticle*>(fcnvec[ifcn])->getParticles();
  for ( unsigned kk=0 ; kk<jetvec.size() ; kk++ ) {
    int jjet=jetvec[kk]->ext<JetIndex>();
    jets_of_fcn->push_back(jjet);
  }
  return *jets_of_fcn;
  //

} 
const IntVec&  Adjusted_TrueJet_Parser::jets_of_initial_cn( int iicn ) {
  intvec->clear();
  IntVec* jets_of_icn=intvec;
  // way to find: jet-to-icn link icn->reco : jets, find notthis.
  //  index<->jet : LCExtension
  //LCObjectVec jetvec = relfcn->getRelatedFromObjects(  finalcns->at(ifcn) );
  ReconstructedParticleVec jetvec= dynamic_cast<ReconstructedParticle*>(initialcns->at(iicn))->getParticles();
  for ( unsigned kk=0 ; kk<jetvec.size() ; kk++ ) {
    int jjet=jetvec[kk]->ext<JetIndex>();
    jets_of_icn->push_back(jjet);
  }
  return *jets_of_icn;
  //

} 

const IntVec&  Adjusted_TrueJet_Parser::pdg_icn_comps(int iicn) {
  intvec->clear();
  for ( unsigned ipid=1 ; ipid <  initialcns->at(iicn)->getParticleIDs().size() ; ipid++ ) {
    intvec->push_back( initialcns->at(iicn)->getParticleIDs()[ipid]->getPDG());
  }
  return *intvec;
  //

} 

const IntVec&  Adjusted_TrueJet_Parser::type_icn_comps(int iicn) {
  intvec->clear();
  for ( unsigned ipid=1 ; ipid <  initialcns->at(iicn)->getParticleIDs().size() ; ipid++ ) {
    intvec->push_back( initialcns->at(iicn)->getParticleIDs()[ipid]->getType());
  }
  return *intvec;
  //

}


const double* Adjusted_TrueJet_Parser::p4_icn(int iicn) { 
  //double* p4=0; 
  const double* mom = p_icn(iicn);
   for (int kk=1 ; kk<=3 ; kk++ ) {
     p4[kk]=mom[kk-1];
   }
   p4[0]=E_icn(iicn);
   return p4 ;
}

const IntVec&  Adjusted_TrueJet_Parser::pdg_fcn_comps(int ifcn) {
  intvec->clear();
  for ( unsigned ipid=1 ; ipid <  finalcns->at(ifcn)->getParticleIDs().size() ; ipid++ ) {
    intvec->push_back( finalcns->at(ifcn)->getParticleIDs()[ipid]->getPDG());
  }
  return *intvec;
  //

} 
const IntVec&  Adjusted_TrueJet_Parser::type_fcn_comps(int ifcn) {
  intvec->clear();
  for ( unsigned ipid=1 ; ipid <  finalcns->at(ifcn)->getParticleIDs().size() ; ipid++ ) {
    intvec->push_back( finalcns->at(ifcn)->getParticleIDs()[ipid]->getType());
  }
  return *intvec;
  //

} 



const double* Adjusted_TrueJet_Parser::p4_fcn(int ifcn) { 
  //double* p4=0; 
  const double* mom = p_fcn(ifcn);
   for (int kk=1 ; kk<=3 ; kk++ ) {
     p4[kk]=mom[kk-1];
   }
   p4[0]=E_fcn(ifcn);
   return p4 ;
}

MCParticle* Adjusted_TrueJet_Parser::initial_elementon(int ijet){
//  int icn=initial_cn(ijet);
//  LCObjectVec elementonvec = relip->getRelatedToObjects(initialcns->at(icn));
  LCObjectVec elementonvec = relip->getRelatedToObjects(jets->at(ijet));
  if (elementonvec.size() > 0 ) {
    MCParticle* mcp  = dynamic_cast<MCParticle*>(elementonvec[0]);
    return mcp;
  } else {
    return NULL;
  }

}
MCParticle* Adjusted_TrueJet_Parser::final_elementon(int ijet){
//  int fcn=final_cn(ijet);
//  LCObjectVec elementonvec = relfp->getRelatedToObjects(finalcns->at(fcn));
  LCObjectVec elementonvec = relfp->getRelatedToObjects(jets->at(ijet));
  if (elementonvec.size() > 0 ) {
    MCParticle* mcp  = dynamic_cast<MCParticle*>(elementonvec[0]);
    return mcp;
  } else {
    return NULL;
  }

}
const MCParticleVec& Adjusted_TrueJet_Parser::elementons_final_cn(int ifcn){
  mcpartvec->clear();
  
  MCParticleVec* elementons= mcpartvec;
  IntVec jetind=jets_of_final_cn(ifcn );
  for (unsigned kk=0 ; kk <jetind.size() ; kk++ ) {
    MCParticle* mcp=final_elementon(jetind[kk]);
    if ( mcp != NULL ) {
      elementons->push_back(mcp);
    }
  }
  return *elementons;
}

const MCParticleVec& Adjusted_TrueJet_Parser::elementons_initial_cn(int iicn){
  mcpartvec->clear();

  MCParticleVec* elementons= mcpartvec;
  IntVec jetind=jets_of_initial_cn(iicn );
  for (unsigned kk=0 ; kk <jetind.size() ; kk++ ) {
    MCParticle* mcp=initial_elementon(jetind[kk]);
    if ( mcp != NULL ) {
      elementons->push_back(mcp);
    }
  }
  return *elementons;
}

void Adjusted_TrueJet_Parser::getall( LCEvent * event ) { 


    evt=event;
     // get TrueJets
    try{
      tjcol = evt->getCollection( _trueJetCollectionName);
    }
    catch( lcio::DataNotAvailableException e )
    {
      streamlog_out(WARNING) <<    _trueJetCollectionName  << " collection not available 1" << std::endl;
      tjcol = NULL;
    }

    jets=getJets();

     // get  FinalColourNeutrals
    try{
      fcncol = evt->getCollection( _finalColourNeutralCollectionName);
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) <<    _finalColourNeutralCollectionName << " collection not available" << std::endl;
        fcncol = NULL;
    }

    finalcns=getFinalcn();

     // get  InitialColourNeutrals
    try{
      icncol = evt->getCollection( _initialColourNeutralCollectionName);
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) <<    _initialColourNeutralCollectionName << " collection not available" << std::endl;
        icncol = NULL;
    }

    initialcns=getInitialcn();

     // get  FinalColourNeutralLink
    LCCollection* fcnlcol = NULL;
    try{
      fcnlcol  = evt->getCollection(  _finalColourNeutralLink );
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) <<  _finalColourNeutralLink   << " collection not available" << std::endl;
        fcnlcol  = NULL;
    }
    relfcn = new LCRelationNavigator( fcnlcol );

     // get  InitialColourNeutralLink
    LCCollection* icnlcol = NULL;
    try{
      icnlcol  = evt->getCollection(  _initialColourNeutralLink );
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) <<  _initialColourNeutralLink  << " collection not available" << std::endl;
        fcnlcol  = NULL;
    }
    relicn = new LCRelationNavigator( icnlcol );

     // get  FinalElementonLink
    LCCollection* fplcol = NULL;
    try{
      fplcol  = evt->getCollection(  _finalElementonLink );
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) <<  _finalElementonLink   << " collection not available" << std::endl;
        fcnlcol  = NULL;
    }
    relfp = new LCRelationNavigator( fplcol );

     // get  InitialElementonLink
    LCCollection* iplcol = NULL;
    try{
      iplcol  = evt->getCollection(  _initialElementonLink );
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) <<  _initialElementonLink   << " collection not available" << std::endl;
        fcnlcol  = NULL;
    }
    relip = new LCRelationNavigator( iplcol );

     // get  TrueJetPFOLink
    LCCollection* tjrecolcol = NULL;
    try{
      tjrecolcol  = evt->getCollection(  _trueJetPFOLink );
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) <<  _trueJetPFOLink   << " collection not available" << std::endl;
        fcnlcol  = NULL;
    }
    reltjreco = new LCRelationNavigator( tjrecolcol );


     // get  TrueJetMCParticleLink
    LCCollection* tjmcplcol = NULL;
    try{
      tjmcplcol  = evt->getCollection(  _trueJetMCParticleLink );
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) <<  _trueJetMCParticleLink   << " collection not available" << std::endl;
        fcnlcol  = NULL;
    }
    reltjmcp = new LCRelationNavigator( tjmcplcol );
    reltrue_tj =NULL ;

}
void Adjusted_TrueJet_Parser::delall( ) { 
    if (  relfcn!= NULL ) delete relfcn;
    if (  relicn!= NULL ) delete relicn;
    if (  relfp!= NULL ) delete relfp;
    if (  relip!= NULL ) delete relip;
    if (  reltjreco != NULL) delete reltjreco;
    if (  reltjmcp != NULL) delete reltjmcp;
    if (  jets != NULL) delete  jets;
    if (  finalcns != NULL) delete  finalcns;
    if (  initialcns!= NULL ) delete   initialcns;
    if ( reltrue_tj != NULL ) delete reltrue_tj;
}
