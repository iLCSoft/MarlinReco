//Filename RootBCalTagEfficiency.cc

#include "BCalTagEfficiency.h"
#include <iostream>
#include <cmath>
#include <stdio.h>


#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud2D.h>
#include <AIDA/IHistogram2D.h>
#endif

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <EVENT/LCRelation.h>
#include <UTIL/LCRelationNavigator.h>

#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/BField.h>
#include "gear/CalorimeterParameters.h"

// include what BcEnergyDensity class needs
//#include "BcEnergyDensity.h"

// include root headers
#include "TRandom3.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio ;
using namespace marlin ;
using namespace std;
using namespace CLHEP ;

// Declaring the function to call FORTRAN routine


extern   "C"  {

  void bcalhit_(float pin_[3], float& q_,float vin_[3] ,
 float& B_, float& eBeam_, float& zbcal_, float pout_[3], float vout_[3]);
} 


BCalTagEfficiency aBCalTagEfficiency ;


BCalTagEfficiency::BCalTagEfficiency() : Processor("BCalTagEfficiency") {

  // modify processor description
  _description = "ForwardVeto propagates electrons or photons from IP to the BeamCal, and it gives the Energy and Momentum of the detected ones to the output ReconstructedParticles collection" ;
  

  // register steering parameters: name, description, class-variable,
  // default value

  registerInputCollection( LCIO::MCPARTICLE,
			   "CollectionName" , 
			   "Name of the MCParticle collection"  ,
			   _MCParticleName ,
			   std::string("MCParticlesSkimmed") ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "BCALInputTruthLinkName" , 
			   "Name of the truth link input collection (from SGV or Mokka)"  ,
			   _BCALInputTruthLinkName ,
			   std::string("BCALMCTruthLink") ) ;

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "BCALParticleName" ,
                           "Name of the tagged e/gamma collection"  ,
                           _BCALParticleName ,
                           std::string("BCALEffiParticles") ) ;

  registerOutputCollection( LCIO::CLUSTER,
                           "BCALClusterName" ,
                           "Name of the collection of clusters of tagged e's/gammas"  ,
                           _BCALClusterName ,
                           std::string("BCALEffiClusters") ) ;


  registerOutputCollection( LCIO::LCRELATION,
			    "BCALEffiMCTruthLinkName",
                           "Name of the tagged e/gamma  to mc-particle relation collection"  ,
			    _BCALEffiMCTruthLinkName, 
                           std::string("BCALEffiMCTruthLink") ) ;


  registerProcessorParameter("BackgroundFilename" ,
                             "Name of input file for background energy density"  ,
                             backgroundfilename ,
                             std::string("bg_aver_LDC_4T_14mrad_AntiDID.root") ) ;

  // eBeam = Energy of the beam in CMS (GeV)
  registerProcessorParameter("eBeam" ,
                             "" ,
                             eBeam,
                             float(250));


  // zbcal = z point from the IP (mm)
  registerProcessorParameter("zbcal" ,
                             "position of the BeamCal" ,
                             zbcal,
                             float(3594.9));

  registerProcessorParameter("thresholdmin",
                             "Minimum energy of detected particles",
                             thresholdMin,
                             float(40));


  registerProcessorParameter("thresholdmax",
                             "Maximum energy of detected particles" ,
                             thresholdMax,
                             float(260));


  registerProcessorParameter("smearEnergy",
                             "true: electrons/photons written to new collection with smeared energy (and momentum)",
                             smearEnergy,
                             bool(1));

  registerProcessorParameter("detectAll",
                             "true: all electrons/photons in BCAL acceptance will be written to new collection",
                             detectAll,
                             bool(1));

  registerProcessorParameter("densityScaling",
                             "Scaling of the energy density w.r.t default setup",
                             densityScaling,
                             float(1.));
                             
  registerProcessorParameter("DBDsample",
                             "true: assume DBD-style event sample with MCParticles with crossing angle, false: LoI-Style, MCParticles head-on",
                             DBDsample,
                             bool(1));

  registerProcessorParameter("newMap",
                             "true: new map with bcal layers 1 - 30, pair mon layer 31; false: old map with pair monitor in layer1",
                             newMap,
                             bool(1));

  registerProcessorParameter("writeTree",
                             "true: write root tree with debug info",
                             writeTree,
                             bool(1));

 registerProcessorParameter("writeSGVMap",
                            "true: write ascii map for SGV",
                             writeSGVMap,
                             bool(0));

 registerProcessorParameter("SGVMapFilename",
                            "name of file for SGV map",
                             SGVmapfilename,
                             std::string("SGVmap.root") );

 registerProcessorParameter("useBCALInputTruthLink",
                            "true: read BCALMCTruthLink from SGV / Mokka to find correct MCParticles",
                             useInputClusters,
                             bool(0));

}

void BCalTagEfficiency::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  // instanciate background handler object
  bc_en = new BcEnergyDensity(backgroundfilename.c_str());
  
  streamlog_out(DEBUG) << "   background initialized   " 
		       << std::endl ;
                       
  bField = Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();
  streamlog_out(MESSAGE) << "BCalTagEfficiency::init: B-FIELD = " << bField << std::endl;
  
  //read parameters from gear file
  const gear::CalorimeterParameters& bcparam = Global::GEAR->getBeamCalParameters(); //in bcparam are read all parameters from the gearfile; the gear file is already opened by default by Marlin (the name is known from the Marlin steering file) 
  xingangle = bcparam.getDoubleVal("beam_crossing_angle");
   
  // parameters of the Lorentz transformation / rotation matrix
  alpha = xingangle/2000;  // half x-ing angle in radians
  gamma = sqrt(1 + tan(alpha)*tan(alpha));
  betagamma = tan(alpha);
        
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
  PI = 4.*atan(1.);

  if (writeTree) { 
  
    // root stuff
    rootfile = new TFile("BCalEffi.root","recreate");
    tree = new TTree("Treename","Tree for electrons");

    tree->Branch("mcp",    &mcp,"mcp/I");
    tree->Branch("pdg",    pdg,"pdg[mcp]/I");
    tree->Branch("energy", energy,"energy[mcp]/F");
    tree->Branch("ePrime", ePrime,"ePrime[mcp]/F");
    tree->Branch("pxPrime",pxPrime,"pxPrime[mcp]/F");
    tree->Branch("pxIP",   pxIP,"pxIP[mcp]/F");
    tree->Branch("pyIP",   pyIP,"pyIP[mcp]/F");
    tree->Branch("pzIP",   pzIP,"pzIP[mcp]/F");
    tree->Branch("scaleP", scaleP,"scaleP[mcp]/F");
    tree->Branch("phiIP",  phiIP,"phiIP[mcp]/F");
    tree->Branch("theIP",  theIP,"theIP[mcp]/F");
    tree->Branch("lposx",  lposx,"lposx[mcp]/F");
    tree->Branch("lposy",  lposy,"lposy[mcp]/F");
    tree->Branch("lposz",  lposz,"lposz[mcp]/F");
    tree->Branch("gposx",  gposx,"gposx[mcp]/F");
    tree->Branch("gposy",  gposy,"gposy[mcp]/F");
    tree->Branch("gposz",  gposz,"gposz[mcp]/F");
    tree->Branch("radius", radius,"radius[mcp]/F");
    tree->Branch("phi",    phi,"phi[mcp]/F");
    tree->Branch("ebkg",   ebkg,"ebkg[mcp]/F");
    tree->Branch("ebkg_err",ebkg_err,"ebkg_err[mcp]/F");
    tree->Branch("efficiency", efficiency,"efficiency[mcp]/F");
    tree->Branch("rand" ,  rand,"rand[mcp]/F");
    tree->Branch("tag" ,   tag,"tag[mcp]/I");
  }

  if (writeSGVMap) { 
    nbinx = 50;
    xmin = -140.;
    xmax = 160.;
    nbiny = 50;
    ymin = -150.;
    ymax = 150.;
    TFile *SGVmaprootfile = new TFile(SGVmapfilename.c_str(),"recreate");
    SGVmapP = new TH2D("hmapP","mapP", nbinx, xmin, xmax, nbiny, ymin, ymax);
    SGVmapN = new TH2D("hmapN","mapN", nbinx, xmin, xmax, nbiny, ymin, ymax);
    SGVmapP->Sumw2();
    SGVmapN->Sumw2();
    double wx = (xmax-xmin)/nbinx;
    double wy = (ymax-ymin)/nbiny;
    double vloc[3]; // local coordinates
    for (int iz = -1; iz < 2; iz += 2) {
      float halfxingangle = iz * alpha;   
      for (int ix = 0; ix < nbinx; ++ix) {
      // bin center in x
        double x = xmin + (ix+0.5)*wx; 
        // convert to local: rotation around y axis of -7 mrad. 
        // rotation goes in different directions for +z / -z
        vloc[0] = cos(-halfxingangle)* x + sin(-halfxingangle)* (iz*zbcal);
        vloc[2] = -sin(-halfxingangle)* x + cos(-halfxingangle)* (iz*zbcal);
        for (int iy = 0; iy < nbiny; ++iy) {
          // bin center in y
          vloc[1] = ymin + (iy+0.5)*wy; 
          // Change from cartesian to cylindrical coordinates !LOCAL!
          Hep3Vector lvec (vloc[0], vloc[1], vloc[2]);
          double rloc = lvec.perp();
          double philoc = lvec.phi(); // gives [-PI, PI]
          if (philoc < 0) philoc += 2*PI;  // have [0, 2 PI] now
          
          double sumEDens = 0;
          double sumEDensErr = 0;
          double EDens = 0;
          double EDensErr = 0;
          bool ok = false;
          
          for (int layer=1; layer<31; layer++) {
            if (!newMap && layer ==1) continue;  // skip layer if using old map!
            ok = bc_en->GetEnergyDensity( layer * iz, rloc, philoc, &EDens, &EDensErr);
            if (ok) {
              sumEDens += EDens;
              sumEDensErr += EDensErr/55.*EDensErr/55.;
            } 
            // in keyhole region
            else {
              sumEDens = -1.;
              sumEDensErr = -1.;
              layer = 31;
            } 
          } 
          sumEDensErr = (sumEDensErr>=0) ? sqrt(sumEDensErr) : -1.;


          //Scaling of the Energy density with scaling factor
          if (sumEDens > 0) {
            sumEDens *= densityScaling;
            sumEDensErr *= densityScaling;
          }  
          if (iz > 0) SGVmapP->Fill(x,vloc[1],sumEDens);
          else SGVmapN->Fill(x,vloc[1],sumEDens);
	  streamlog_out(MESSAGE) << "SGVmap " << x << "  " << vloc[1] << "  " << iz*zbcal << "  " << sumEDens << "  " << sumEDensErr << std::endl;
        }
      }
    }
    
    SGVmapP->Write();
    SGVmapN->Write();
    SGVmaprootfile->Write();
    SGVmaprootfile->Close();
    //delete SGVmaprootfile;
    //delete SGVmapP;
    //delete SGVmapN;
    
  } // if write SGV map
  
}

void BCalTagEfficiency::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void BCalTagEfficiency::processEvent( LCEvent * evt ) { 

  streamlog_out(DEBUG) << "+++++++++++++++++ start processing event: " << evt->getEventNumber() 
                    << "   in run:  " << evt->getRunNumber() 
                    << std::endl ;
 
  // Create output collection
  LCCollectionVec* BCALCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  //LCCollectionVec* BCALEffiMCTruthLink = new LCCollectionVec(LCIO::LCRELATION);
  LCCollection* BCALEffiMCTruthLink = 0;
  LCCollectionVec* BCALClusters = new LCCollectionVec(LCIO::CLUSTER);
  LCRelationNavigator bcalTruthRelNav(LCIO::RECONSTRUCTEDPARTICLE , LCIO::MCPARTICLE  ) ;

  LCCollection* col = 0;
  // read MCParticles and track to BeamCal
  if (!useInputClusters) {
    streamlog_out(DEBUG) << "Using MCParticle collection!" << std::endl;
    col = evt->getCollection( _MCParticleName ) ;
    if (!col) streamlog_out(ERROR) << "MCParticle collection not found!!!" << std::endl;
  }
  // loop over BCalParticles, find corresponding MCParticle, track to BeamCal
  else {
    streamlog_out(DEBUG) << "Using BCALTruthLink collection!" << std::endl;
    col = evt->getCollection( _BCALInputTruthLinkName ) ;
    if (!col) streamlog_out(ERROR) << "TruthLinker not found!!!" << std::endl;
  }
 
  if( col != 0 ){
   
    int nMCP = col->getNumberOfElements()  ;
    
    // Arguments to be called from the FORTRAN routine
    float q =0;
    float pin[3]= {0, 0, 0};    // momenta at IP, !in coord sys of input! 
    float vin[3]= {0, 0, 0};    // coordinates of PCA, !in coord sys of input! 
    float gvout[3]= {0, 0, 0};	 // coordinates at BeamCal surface, global coord sys 
    float lvout[3]= {0, 0, 0};	 // coordinates at BeamCal surface, local coord sys 
    float pout[3]= {0, 0, 0};    // momenta at BeamCal surface !in coord sys of input! WILL NOT BE USED!
   
    mcp = 0; 

    for(int i=0; i < nMCP ; i++){
   
      MCParticle* p = 0;
      LCRelation* rp = 0;
   
      if (!useInputClusters) {
        // take all genStat = 1 electrons & photons
        p = dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;
        if ( p->getGeneratorStatus()!=1 ) continue;
        pdg[mcp] = p->getPDG();
        if ( abs(pdg[mcp]) != 11 && abs(pdg[mcp]) != 22 ) continue;
        streamlog_out(DEBUG) << "found MCP..." <<  std::endl;
      } 
      else {
        streamlog_out(DEBUG) << "looking for MCP..." <<  std::endl;
        // pick out correct MCP using SGV / Mokka information!
        rp =  dynamic_cast<LCRelation*>( col->getElementAt( i ) ) ;
        if (!rp) {
          streamlog_out(WARNING) << " entry in LCRelation not found! " << std::endl;
          continue;
        } 
        if (rp->getWeight()<0.5) {
          streamlog_out(WARNING) << " relation weight low:  " << rp->getWeight() << std::endl;
          continue;
        }  
        p =  dynamic_cast<MCParticle*>(rp->getTo());
        streamlog_out(DEBUG) << "found MCP with PDG = " << p->getPDG() ;
        streamlog_out(DEBUG) << " and energy = " << p->getEnergy() << std::endl;
      }
      
      // Get initial energy 
      energy[mcp] = p->getEnergy();

      // Only use the electron if it's energy > 40 GeV
      if((energy[mcp] > thresholdMin && energy[mcp] < thresholdMax)) {
      
        // get momentum
        pin[0] = p->getMomentum()[0];
        pin[1] = p->getMomentum()[1];
        pin[2] = p->getMomentum()[2];

        pxIP[mcp] = pin[0];
        pyIP[mcp] = pin[1];
        pzIP[mcp] = pin[2];
      
        // this is with / without x-ing angle for DBD / LoI samples
        HepLorentzVector p4v(pin[0], pin[1], pin[2], energy[mcp]);
        phiIP[mcp] = p4v.phi();
        streamlog_out(DEBUG) << "phi  of electron/photon before +=2PI " << phiIP[mcp] << std::endl;
        if (phiIP[mcp] < 0) phiIP[mcp] += 2*PI;
        theIP[mcp] = p4v.theta();
        
        // theta larger than 0.06 (r/z = 200mm / 3500m = 0.057) will go into main detector
        if (theIP[mcp] > 0.06 && theIP[mcp] < PI-0.06 || theIP[mcp] < 0.003 || theIP[mcp] > PI-0.003) continue;  
        streamlog_out(DEBUG) << "energy  of electron/photon = " << energy[mcp] << std::endl;
        streamlog_out(DEBUG) << "phi at IP of electron/photon = " << phiIP[mcp] << std::endl;
        streamlog_out(DEBUG) << "pt of electron/photon = " << p4v.perp() << std::endl;

        // gets de particle charge  
        q = p->getCharge();

        // call to the FORTRAN routine bcalhit        
        // gives pout and vout 

        if( abs(pdg[mcp]) == 11) {   // electron, full B field tracking
          bcalhit_(pin,q,vin,bField,eBeam,zbcal,pout,gvout);      
        }
        else {  // photon  - why do we have to treat this separately?
          gvout[0]=(pin[0]/abs(pin[2]))*zbcal;
          gvout[1]=(pin[1]/abs(pin[2]))*zbcal; 
          gvout[2]=(pin[2]/abs(pin[2]))*zbcal;
          pout[0]=pin[0];
          pout[1]=pin[1];
          pout[2]=pin[2];
        }
        
        // rotation goes in different directions for +z / -z
        float halfxingangle = pin[2]>0 ? alpha : -alpha;   
        
        // now sort out whether we have crossing angle already in MCParticles or not
        if (DBDsample) {
          // calculate local BeamCal coordinates in order to compare with background map
          // Rotation around y axis of -7 mrad. (back from global to local!)
          lvout[0] = cos(-halfxingangle)* gvout[0] + sin(-halfxingangle)* gvout[2];
          lvout[1] = gvout[1];
          lvout[2] = -sin(-halfxingangle)* gvout[0] + cos(-halfxingangle)* gvout[2];
 
          // don't need lpout? -> check! 
        }
        else {
          // move gvout to local variables 
          // (global coordinates without crossing-angle correspond directly to local 
          //  BeamCal coordinates!)
          // calculate global coordinates 
          
          lvout[0] = gvout[0];
          lvout[1] = gvout[1];
          lvout[2] = gvout[2];
          
          // Rotation around y axis of 7 mrad. (forward from local to global!)
          gvout[0] = cos(halfxingangle)* lvout[0] + sin(halfxingangle)* lvout[2];
          gvout[1] = lvout[1];
          gvout[2] = -sin(halfxingangle)* lvout[0] + cos(halfxingangle)* lvout[2];
                      
        } 
         
        streamlog_out(DEBUG) << "entry point on BeamCal surface, global ILD coordinates:" << std::endl;
        streamlog_out(DEBUG) << "gvout[0] = " << gvout[0] << ", gvout[1] = " << gvout[1] << ", gvout[2] = " << gvout[2] << std::endl;
        streamlog_out(DEBUG) << "entry point on BeamCal surface, local coordinates:" << std::endl;
        streamlog_out(DEBUG) << "lvout[0] = " << lvout[0] << ", lvout[1] = " << lvout[1] << ", lvout[2] = " << lvout[2] << std::endl;

        
        // Change from cartesian to cylindrical coordinates !LOCAL!
        Hep3Vector lvvec (lvout[0], lvout[1], lvout[2]);
        
        radius[mcp] = lvvec.perp();
        streamlog_out(DEBUG) << "radius on BeamCal surface in local coordinates = " << radius[mcp] << std::endl;
        
        phi[mcp] = lvvec.phi(); // gives [-PI, PI]
        streamlog_out(DEBUG) << "phi on BeamCal surface in local coordinates = " << phi[mcp] << std::endl;
        if (phi[mcp] < 0) phi[mcp] += 2*PI;  // have [0, 2 PI] now
        streamlog_out(DEBUG) << "phi after +2PI = " << phi[mcp] << std::endl;
        if (p4v.perp() < 1E-3) phi[mcp] = 0;
                          
        
        // just debugging!
        Hep3Vector gvvec (gvout[0], gvout[1], gvout[2]);
        float gphi = gvvec.phi(); // gives [-PI, PI]
        streamlog_out(DEBUG) << "gphi on BeamCal surface in global coordinates = " << gphi << std::endl;
        if (gphi < 0) gphi += 2*PI;  // have [0, 2 PI] now
        streamlog_out(DEBUG) << "gphi after +2PI = " << gphi << std::endl;
        float gtheta = gvvec.theta();
        streamlog_out(DEBUG) << "gtheta on BeamCal surface in global coordinates = " << gtheta << std::endl;
        streamlog_out(DEBUG) << "radius on BeamCal surface in global coordinates = " << gvvec.perp()*lvvec.mag()/gvvec.mag() << std::endl;
                          
        // Loop to sum the energy density over all layers.
        // ebkg = Energy of background
        // If res == 1 GetEnergyDensity succeeded
          
        ebkg[mcp] = 0;
        ebkg_err[mcp] = 0;
        bool res = 0;
        double en_dens = 0;
        double en_dens_err = 0;

        for (int layer=1; layer<31; layer++) {
        if (!newMap && layer ==1) continue;  // skip layer if using old map!

          int sign = 1; 
          if(lvout[2] < 0) sign = -1;
          res = bc_en->GetEnergyDensity( layer * sign, radius[mcp], phi[mcp], &en_dens, &en_dens_err);
          streamlog_out(DEBUG) << "energy density in layer " << layer * sign << "  = " << en_dens << std::endl;
 
          if (res) {
            ebkg[mcp] += en_dens;
            ebkg_err[mcp] += en_dens_err*en_dens_err;
          }  
          
        } 
        if (ebkg_err[mcp]>=0) {
          ebkg_err[mcp] = sqrt(ebkg_err[mcp]);
        }
        else {
          ebkg_err[mcp] = 0;
        }  
        streamlog_out(DEBUG) << "background energy density sum = "    << ebkg[mcp] << std::endl;


        //Scaling of the Energy density with scaling factor
        ebkg[mcp] *= densityScaling;
        
        if (densityScaling != 1.)  {  
           streamlog_out(DEBUG) << "scaling background energy density sum with factor "  << densityScaling 
                                << ", scaled density = "  << ebkg[mcp] <<std::endl;
        }
                  

        // Here starts the Eff_routine that parametrizes the efficiency
         
        //assuming that efficiency is parametrised in GeV/cm^2 instead of GeV/mm^2 ...
        //alternative: efficiency is in GeV / cell -> *= 55
        ebkg[mcp] *= 100.;

        efficiency[mcp] = 0;

        if ( ebkg[mcp] <= 0 && useInputClusters) {  // if == 0, then particle not in BCAL acceptance!
          
          streamlog_out(DEBUG) << "Particle made SGV cluster,  but ebkg[mcp] = " << ebkg[mcp] 
                               << ", so not in map: radius[mcp] = " << radius[mcp] 
                               << ", phi[mcp] = " <<  phi[mcp]*180/3.1415 << std::endl;
        }                        
         
        else if ( ebkg[mcp] > 0 ) {  // if == 0, then particle not in BCAL acceptance!
        
          float p0, p1, p2;
        
          //energy[mcp] = sqrt(pow(pout[0],2)+pow(pout[1],2)+pow(pout[2],2));

          if (ebkg[mcp]>=0 && ebkg[mcp]<35 && energy[mcp]>45 && energy[mcp]<345) {
            if (energy[mcp]>255) {
              streamlog_out(WARNING) << "not in original range of parametrisation , energy = " << energy[mcp] << std::endl;
            }

            p0 = -3.2754e-2 + 8.62e-4*energy[mcp] -2.4424e-6*energy[mcp]*energy[mcp];
            p1 = -8.936e-3 + 6.05e-4*energy[mcp] -1.719e-6*energy[mcp]*energy[mcp];
            p2 = 6.112e-2 + 2.84e+3/(energy[mcp]*energy[mcp]);

            efficiency[mcp] = p0 + 1/(1+p1*exp(p2*ebkg[mcp]));
            
            if (efficiency[mcp] > 1.) efficiency[mcp] = 1.;
            if (efficiency[mcp] < 0.) efficiency[mcp] = 0.;
        
          }
          else if (ebkg[mcp]>=0 && ebkg[mcp]<35 && energy[mcp]>45 && energy[mcp]<550) {
            efficiency[mcp] = 1.;        
          }
          else {
        
            streamlog_out(ERROR) << "not in valid range of parametrisation , ebkg = "  
                                 << ebkg[mcp] << ", energy = " << energy[mcp] << std::endl;
        
          }
     
          // End of the Eff_routine
   
          // Selection algorithm

          char energystring [50];
          char charseed[5];
          int seed;

          // seed setting 
          sprintf(energystring,"%f",energy[mcp]);
          
          charseed[0]=energystring[4];
          charseed[1]=energystring[5];
          charseed[2]=energystring[6];
          charseed[3]=energystring[7];
          charseed[4]=energystring[8];

          seed = atoi(charseed);

          TRandom3 *randGen = new TRandom3(seed);

          rand[mcp] = randGen->Rndm();

          tag[mcp] = 0;
          if (rand[mcp] < efficiency[mcp] || (detectAll == 1 && efficiency[mcp]>0.0) ) {
          
            ReconstructedParticleImpl* particle = new ReconstructedParticleImpl;
            ClusterImpl* cluster = new ClusterImpl;
            streamlog_out(DEBUG) << "adding relation between particle " << particle << "and MCParticle " << p << endl;
            bcalTruthRelNav.addRelation( particle , p , 1.0) ;
            //LCRelationImpl* MCrel  = new LCRelationImpl(particle,p,0.5);
            //MCrel->setFrom (particle);
            //MCrel->setTo (p);
            
            //FG: take generator mass ( not necessarily equal to PDG mass ) 
            const double m  = p->getMass() ;
            particle->setMass(m);
            particle->setCharge(q);
            
            if (DBDsample) {
              pxPrime[mcp] = pxIP[mcp];  // no boost needed
            }
            else {
              // DO BOOST! 
              // before the transformation: p4v 
              // (will store 4 vector at IP, not BCAL surface!)
              // the following is taken from Mokka:
              // PrimaryGeneratorAction::ApplyLorentzTransformation
              pxPrime[mcp] = betagamma * energy[mcp] + gamma * pxIP[mcp];
            }
            // py and pz remain the same, E changes implicitly with px

            scaleP[mcp] = 1.;
            
            // energy resolution:
            // dE = 25% * sqrt(E) (intrinsic resolution) + 15% * E_bgk_dens * 55mm^2 (background fluctuations)
            // alternative: dE = 130% sqrt(E) in total
            //double sigmaE = 1.3*sqrt(energy[mcp]);
            //double sigmaE = 0.25*sqrt(energy[mcp]) + 0.15 * ebkg[mcp] * 55; // ????
            double sigmaE = 0.25*sqrt(energy[mcp]) + ebkg_err[mcp] * 55;
            if (smearEnergy == true) { 
              double deltaE = randGen->Gaus() * sigmaE;
              scaleP[mcp] = (1+deltaE/energy[mcp]);
              if (scaleP[mcp] < 0) scaleP[mcp] = 0;
            }  

            double pnew[3] = {scaleP[mcp]*pxPrime[mcp], scaleP[mcp]*pyIP[mcp], scaleP[mcp]*pzIP[mcp]};
            particle->setMomentum(pnew);
            particle->setEnergy(scaleP[mcp]*sqrt(m*m + pxPrime[mcp]*pxPrime[mcp] +pyIP[mcp]*pyIP[mcp] +pzIP[mcp]*pzIP[mcp]));
            
            // store efficiency as "GoodnessOfPID"
            particle->setGoodnessOfPID(efficiency[mcp]);
            
            cluster->setEnergy(particle->getEnergy());
            cluster->setEnergyError(sigmaE);
            cluster->setPosition(gvout);  // global ILD coord
            cluster->setITheta(gtheta);
            cluster->setIPhi(gphi);
            cluster->subdetectorEnergies().resize(6); cluster->subdetectorEnergies()[5] = particle->getEnergy(); 

            particle->addCluster(cluster);

            BCALCol->addElement(particle);
            //BCALEffiMCTruthLink->addElement(MCrel);
            BCALClusters->addElement(cluster);

            // store positions as single values in tree
            lposx[mcp]=lvout[0];  // local BeamCal coord
            lposy[mcp]=lvout[1];
            lposz[mcp]=lvout[2];
            gposx[mcp]=gvout[0];  // global ILD coord
            gposy[mcp]=gvout[1];
            gposz[mcp]=gvout[2];
            tag[mcp] = 1;
            

            ePrime[mcp] = particle->getEnergy();
    

            streamlog_out(DEBUG) << "Detected particle with ID " << pdg[mcp] << endl;
            streamlog_out(DEBUG) << "Energy before smearing = " << energy[mcp] << endl;
            streamlog_out(DEBUG) << "Energy after smearing = " << ePrime[mcp] << endl;
            streamlog_out(DEBUG) << "Efficiency= " << efficiency[mcp] << endl;
            streamlog_out(DEBUG) << "E background= " << ebkg[mcp] << endl;
            streamlog_out(DEBUG) << "+++++++++++++++++++" << endl;
 
          }
     
          mcp++;
          
          if (writeTree) {
            streamlog_out(DEBUG) << "filling tree, mcp = " << mcp << endl;
            tree->Fill();
          }  

          delete randGen;
   
        } // beamcal acceptance
          
      }  // energy thresholds
              
    }  // end of MCparticle loop
    
  }
  
  BCALEffiMCTruthLink = bcalTruthRelNav.createLCCollection() ;
  streamlog_out(DEBUG) << "name of BCALEffiMCTruthLink: " << _BCALEffiMCTruthLinkName << ", number of elements = " << BCALEffiMCTruthLink->getNumberOfElements() << endl;
  streamlog_out(DEBUG) << "name of BCALCol: " << _BCALParticleName << ", number of elements = " << BCALCol->getNumberOfElements() << endl;

  //if ( BCALCol->getNumberOfElements() <= 0 ) {
     //delete  BCALCol; 
     //delete  BCALEffiMCTruthLink; 
     //delete  BCALClusters; 
  //} else {
    evt->addCollection(BCALCol ,_BCALParticleName) ;
    evt->addCollection(BCALEffiMCTruthLink ,_BCALEffiMCTruthLinkName) ;
    evt->addCollection(BCALClusters ,_BCALClusterName) ;
  //}


  streamlog_out(DEBUG) << "+++++++++++++++++ finished processing event: " << evt->getEventNumber() 
                    << "   in run:  " << evt->getRunNumber() 
                    << std::endl ;


  _nEvt ++ ;
}



void BCalTagEfficiency::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void BCalTagEfficiency::end(){ 
  
  delete bc_en;
  
  if (writeTree) {
    tree->Write();
    rootfile->Write();
  }  

}

