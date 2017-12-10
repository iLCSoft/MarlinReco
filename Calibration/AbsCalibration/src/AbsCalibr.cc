#include "AbsCalibr.h"
#include "Calibration.h"

#include <cstdio>
#include <cmath>
#include <cstring>
#include "lcio.h"
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/CalorimeterParameters.h>
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/MCParticle.h"
#include "IMPL/LCCollectionVec.h"


#include <cstdlib>

//  To create collection can be output by LCIO 
#include "UTIL/LCFixedObject.h"

// #define ASCII_OUTPUT
// #define USE_HBOOK

#ifdef USE_HBOOK
//+++++++++++++ HBOOK +++  HPLOT  ++++++++++++++++
#ifndef f2cFortran
#define f2cFortran 1
#endif
#include <cfortran.h>
#include <graflib.h>
//#include <packlib.h>
//#include <kernlib.h>
#include <hbook.h>
#include <higz.h>
#define PAWC_SIZE 5000000
typedef struct { float PAW[PAWC_SIZE]; } PAWC_DEF;
#define PAWC COMMON_BLOCK(PAWC,pawc)
COMMON_BLOCK_DEF(PAWC_DEF,PAWC);
//+++++++++++++ HBOOK +++  HPLOT  ++++++++++++++++
#endif


// using namespace std ;
using namespace EVENT ;
using namespace IMPL ;
using namespace UTIL;

#ifdef USE_HBOOK
//            +++++++++++++ HBOOK +++  HPLOT  ++++++++++++++++
//============================================================================
static void plot_init(){ // Open HIGZ windows and Initialize HBOOK
//============================================================================
  static bool is_init=0;
  if(is_init)  return;
  HLIMIT(PAWC_SIZE);
  HPLINT(4);
  is_init=1;
}  //   End of plot_init()

//============================================================================
static void hist_init(){ // Book Histograms
//============================================================================
  static bool is_init=0;
  if(is_init)  return;
  
  
  //             1000 GeV ECM
//  HBOOK1(1,"ECAL hit energy sum ",700,0.0,1200.,0.); 
//  HBOOK1(2,"HCAL hit energy sum ",700,0.0,1200.,0.); 
//  HBOOK1(3,"Calorimeter hit energy sum ",700,0.0,1200.,0.); 
//  HBOOK2(4,"E-ECAL vs E-HCAL",125,0.0,1000.0, 125,0.0,1000.0, 0.0);
//  HBOOK1(5,"Calorimeter hit energy sum ",350,0.0,1200.,0.); 
//  HBOOK1(6,"Calorimeter energy - Available",120,-120.0,120.,0.); 
//!!!
  //             500 GeV ECM
  HBOOK1(1,"ECAL hit energy sum ",700,0.0,700.,0.); 
  HBOOK1(2,"HCAL hit energy sum ",700,0.0,700.,0.); 
  HBOOK1(3,"Calorimeter hit energy sum ",700,0.0,700.,0.); 
  HBOOK2(4,"E-ECAL vs E-HCAL",125,0.0,500.0, 125,0.0,500.0, 0.0);
  HBOOK1(5,"Calorimeter hit energy sum ",350,0.0,700.,0.); 
  HBOOK1(6,"Calorimeter energy - Available",120,-60.0,60.,0.); 

   //   HPLOPT("GRID",1);
   //   HPLOPT("STAT",1);
  is_init=1;
}   //End of hist_init()
//         +++++++++++++ HBOOK +++  HPLOT  ++++++++++++++++
#endif

namespace marlin{

double  Balance(LCEvent * evt);

  AbsCalibr aAbsCalibr ;
  AbsCalibr::AbsCalibr() : Processor("AbsCalibr") {
    _description = " Calorimeter Abslute Energy Calibration" ;
    vector<int> nlayer;
    nlayer.push_back(30);
    nlayer.push_back(10);
    nlayer.push_back(40);
    registerProcessorParameter("NLayer","Number of layers in zone",_nlayer,nlayer);
    vector<float> coeff;
    coeff.push_back(33.02346);
    coeff.push_back(93.56822);
    coeff.push_back(21.196262);
    registerProcessorParameter("Coeff","Calorimeter coeffs",_coeff,coeff);
    vector<float> cuts;
    cuts.push_back(170.0e-6 * 0.5);
    cuts.push_back(170.0e-6 * 0.5);
    cuts.push_back(875.0e-6 * 0.5);
    registerProcessorParameter("Cuts","Calorimeter cuts",_cuts,cuts);
  }

//============================================================================
  void AbsCalibr::init() { 
//============================================================================
    std::cout << "AbsCalibr::init()  " << name() << std::endl 
	      << "  parameters: " << std::endl << *parameters()  ;



#ifdef USE_HBOOK
//         +++++++++++++ HBOOK +++  HPLOT  ++++++++++++++++
    plot_init();
    hist_init();

    std::cout << 
    "++++++++++++++ HBOOK and HPLOT were initiaslized +++++++++++++++++++ "
    << std::endl;
//         +++++++++++++ HBOOK +++  HPLOT  ++++++++++++++++
#endif

    _nRun = 0 ;
    _nEvt = 0 ;
}

//============================================================================
  void AbsCalibr::processRunHeader( LCRunHeader* run) { 
//============================================================================
    std::cout << "AbsCalibr::processRun()  " << name() 
	      << " in run " << run->getRunNumber() << "; description is " 
	      <<run->getDescription()<<std::endl ;
    _nRun++ ;
  } 
//============================================================================
  void AbsCalibr::processEvent( LCEvent * evt ) { 
//============================================================================

  LCCollectionVec* abc = new LCCollectionVec(LCIO::LCGENERICOBJECT);

  CalorimeterHit* hit; // pointer of class CalorimeterHit

  int nHits[3];
  double    en[3];

  nHits[ECAL1]=nHits[ECAL2]=nHits[HCAL]=0;
  en[ECAL1]=en[ECAL2]=en[HCAL]=0.;

//          Find collection names to fill hits arrayes 
  const std::vector<std::string>* strVec = evt->getCollectionNames() ;
  for( std::vector<std::string>::const_iterator name = strVec->begin() ; 
       name != strVec->end() ; name++) {
    const std::string & tname = evt->getCollection( *name )->getTypeName();
    //    std::cout<<"Collection type: "<<tname<<std::endl;

//++++++++++++++++++++++   ECAL   hit filling   +++++++++++++++++++++++++++++++++++++++
    if(tname == LCIO::CALORIMETERHIT ){
//-----------------------------------------------------------------------
      if(!std::strncmp((*name).data(),"ECAL",4)){ // Resolve collection name
//-----------------------------------------------------------------------
//	std::cout << "  ECAL:     " <<  *name << std::endl;
	std::basic_string<char>  coll = *name;
	EVENT::LCCollection* col_ECAL = evt->getCollection(coll) ;

	unsigned n =  col_ECAL->getNumberOfElements();       // Get number of hits in ECAL

	for( unsigned i = 0 ; i < n ; i++ ){             // then fill it
	  hit = dynamic_cast<CalorimeterHit*>( col_ECAL->getElementAt(i));
	  float e = hit->getEnergy();
	  unsigned keyGE = hit->getCellID0() ;
	  int l = (keyGE & 0x3F000000)>>24;
	  if(l<_nlayer[ECAL1]){
	    if(e>_cuts[ECAL1]){
	      en[ECAL1]+=e*_coeff[ECAL1];
	      nHits[ECAL1]++;
	    }
	  } else {
	    if(e>_cuts[ECAL2]){
	      en[ECAL2]+=e*_coeff[ECAL2];
	      nHits[ECAL2]++;
	    }
	  }  //  end ECAL hit loop
	}     // Resolve collection name

//+++++++++++++++++++++  HCAL hit filling +++++++++++++++++++++++++++++++++++++++++++++
      } else if(!std::strncmp((*name).data(),"HCAL",4)){// Resolve collection name
//-----------------------------------------------------------------------
//	std::cout << "  HCAL:     " <<  *name << std::endl;
	std::basic_string<char> coll = *name;
	EVENT::LCCollection* col_HCAL = evt->getCollection(coll) ;

	nHits[HCAL] =  col_HCAL->getNumberOfElements() ; //Get number of hits in HCAL

	for( int i=0 ; i< nHits[HCAL] ; i++ ){     // then fill it
	  hit = dynamic_cast<CalorimeterHit*>( col_HCAL->getElementAt(i));
	  float e = hit->getEnergy();
	  if(e>_cuts[HCAL])
	    en[HCAL]+=e*_coeff[HCAL];
	}  //  end HCAL hit loop
      }    // Resolve collection name
      else
	std::cout << "  !!! UNKNOWN !!!:     " <<  *name << std::endl;
    } // along calorimeter collections in LCIO
  }   // along all collections in LCIO

  double E_Real = Balance(evt);

//-----------------------------------------------------------------------
//      Create collection to store it in *.slcio
//-----------------------------------------------------------------------
  Calibration* b = new Calibration(_nEvt,nHits[0],nHits[1],nHits[2],en[0],en[1],en[2],E_Real);
  abc->addElement(b);
  evt->addCollection(abc, "AbsCalibr");
//-----------------------------------------------------------------------

#ifdef ASCII_OUTPUT
   static FILE *f=NULL;
   if(!f) f=fopen("abs_calibr.root","w");
   fprintf(f,"%10i %10u %10u %10u %15.3f %15.3f %15.3f %15.3f\n",
	   _nEvt,nHits[ECAL1],nHits[ECAL2],nHits[HCAL],
	   en[ECAL1],en[ECAL2],en[HCAL],E_Real);
#endif



#ifdef USE_HBOOK
/*******
   if (E_Real> 0.5) {
   HF1(1,E_Ecal,1.0);
   HF1(2,E_Hcal,1.0);
   HF1(3,E_Ecal+E_Hcal,1.0);
   HF1(5,E_Ecal+E_Hcal,1.0);
   HF2(4,E_Ecal,E_Hcal,1.0);
   HF1(6,E_Ecal+E_Hcal-E_Real,1.0);
   }
   
      if((_nEvt%20) == 0){
     //  Draw current Histograms 
     HPLOPT("STAT",1);
     // HPLOPT("LOGY",1);
     HPLZON(2,2,1," ");
     HPLOT(5,"BOX ","HIST",0);
     HPLOT(6,"BOX ","HIST",0);
     HPLOT(3,"BOX ","HIST",0);
     HPLOT(4,"BOX ","HIST",0);
     IUWK(0,1);
     HPLOPT("NSTA",1);
        }
******/
#endif
   //   std::cout<< "   ---> Calo energy =  " <<  E_Ecal <<" + "<< E_Hcal
   //	    <<" = "<<  E_Ecal+E_Hcal<<"; Available E = "<<E_Real<<std::endl;

#ifdef USE_HBOOK
//     Wait for next event 
//std::cout << "        ++++++++++ TYPE o TO SAVE HISTOGRAMS, IF YOU WISH, +++++++" << std::endl ;
  //   int c = getchar();
  //   if(c=='o')   HRPUT(0,"abs_calibr.his"," ");  // Save all existing histograms on disk
#endif

  _nEvt ++ ;
  return ;  // from void AbsCalibr::processEvent( LCEvent * evt ) 
}

//============================================================================
void AbsCalibr::check( LCEvent *  /*evt*/ ) { ;}

//============================================================================
void AbsCalibr::end(){ 
//============================================================================
  std::cout << "AbsCalibr::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    << std::endl ;
#ifdef USE_HBOOK
  HRPUT(0,"abs_calibr.his"," ");  // Save all existing histograms on disk
#endif
}

//============================================================================
double Balance(LCEvent * evt){
//============================================================================
  int idpdg;
  const double* mom;
  float enr;
  double mass;
  LCCollection* mcpCol = evt->getCollection("MCParticle" ) ;  
//-----------------------------------------------------------------------
// Calculate balance at IP taking into account everything
//-----------------------------------------------------------------------
  double px,py,pz,pt,ttet;

  double e_to_tube  = 0.;
  double e_to_tubex = 0.;
  double e_to_tubey = 0.;
  double e_to_tubez = 0.;
  int    n_to_tube = 0; 

  double e_neutr = 0.;   
  double e_neutrx= 0.;   
  double e_neutry= 0.;   
  double e_neutrz= 0.;   
  int    n_neutr= 0;    
   
  double e_muon = 0.;  
  double e_muonx= 0.;  
  double e_muony= 0.;  
  double e_muonz= 0.;  
  int    n_muon= 0;   
   
  double e_elect = 0.; 
  double e_electx= 0.; 
  double e_electy= 0.; 
  double e_electz= 0.; 
  int    n_elect= 0;  
   
  double e_photon = 0.;
  double e_photonx= 0.;
  double e_photony= 0.;
  double e_photonz= 0.;
  int    n_photon= 0; 
   
  double e_pi0 = 0.;
  double e_pi0x= 0.;
  double e_pi0y= 0.;
  double e_pi0z= 0.;
  int    n_pi0= 0; 
   
  double e_llhadr = 0.;
  double e_llhadrx= 0.;
  double e_llhadry= 0.;
  double e_llhadrz= 0.;
  int    n_llhadr= 0; 
   
  double e_slhadr = 0.;
  double e_slhadrx= 0.;
  double e_slhadry= 0.;
  double e_slhadrz= 0.;
  int    n_slhadr= 0; 
   
  double e_chadr = 0.; 
  double e_chadrx= 0.; 
  double e_chadry= 0.; 
  double e_chadrz= 0.; 
  int    n_chadr= 0;  

  for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom = imc-> getMomentum (); 
    enr = imc-> getEnergy (); 
    mass = imc-> getMass (); 
    if( imc-> getGeneratorStatus() == 1) { // stable particles only   
      px = mom[0]; 
      py = mom[1]; 
      pz = mom[2];
      pt = hypot(px,py);
      ttet = atan2(pt,pz);
      if ((fabs(ttet) < 0.1) || (fabs(M_PI-ttet) < 0.1)) {
	e_to_tube  += enr;
	e_to_tubex += px;
	e_to_tubey += py;
	e_to_tubez += pz;
	n_to_tube ++;
	continue;
      } 
      if((abs(idpdg)==12)||(abs(idpdg)==14)||(abs(idpdg)==16)) {
	e_neutr  += enr;
	e_neutrx += px;
	e_neutry += py;
	e_neutrz += pz;
	n_neutr ++;
	continue;
      } 
      if(abs(idpdg)==13) { // mu+ mu- 
	e_muon  += enr;
	e_muonx += px;
	e_muony += py;
	e_muonz += pz;
	n_muon ++;
	continue;
      } 
      if(abs(idpdg)==11) { //  e+ e-
	e_elect  += enr;
	e_electx += px;
	e_electy += py;
	e_electz += pz;
	n_elect ++;
	continue;
      } 
      if(idpdg == 111) { // Pi0 as stable 
	e_pi0  += enr;
	e_pi0x += px;
	e_pi0y += py;
	e_pi0z += pz;
	n_pi0 ++;
	continue;
      } 
      if(idpdg == 22) { // photon
	e_photon  += enr;
	e_photonx += px;
	e_photony += py;
	e_photonz += pz;
	n_photon ++;
	continue;
      } 
      if(    // long lived neutral hadrons
	 (abs(idpdg)==2112)|| // neutron
	 (abs(idpdg)== 130)   // KoL
	 ) {                     
	e_llhadr  += enr;
	e_llhadrx += px;
	e_llhadry += py;
	e_llhadrz += pz;
	n_llhadr ++;
	continue;
      }
      if(  // short lived neutral hadrons
	 (abs(idpdg)== 310)|| // KoS
	 (abs(idpdg)==3122)|| // Lambda0
	 (abs(idpdg)==3212)|| // Sigma0
	 (abs(idpdg)==3322)   // Xi0
	 ) {
	e_slhadr  += enr;
	e_slhadrx += px;
	e_slhadry += py;
	e_slhadrz += pz;
	n_slhadr ++;
	continue;
      }
      if(!(abs(idpdg)==12) && !(abs(idpdg)==14) && !(abs(idpdg)==16) &&  // neutrinos   
	 !(abs(idpdg)==13) && // mu+ mu- 
	 !(abs(idpdg)==11) && //  e+ e-
	 !(idpdg == 111) && // Pi0
	 !(idpdg == 22) &&  // photon
	 !(abs(idpdg)==2112) && !(abs(idpdg)== 311) && // neutral hadrons
	 !(abs(idpdg)== 130) && !(abs(idpdg)== 310) && // neutral hadrons
	 !(abs(idpdg)==3122) && !(abs(idpdg)==3212)) { // neutral hadrons
	e_chadr  += enr;
	e_chadrx += px;
	e_chadry += py;
	e_chadrz += pz;
	n_chadr ++;
	continue;
      }
      std::cout <<" Unknow for this program  ID is " <<idpdg<< std::endl;
    }    // Stable particles only
  }      // End for for MCParticles 
  double e_sum = e_elect+e_muon+e_chadr+e_pi0+e_photon+e_llhadr+e_slhadr+e_neutr;
  double evt_energy = e_sum;
  //   Minus muon loosed energy = 1.3 GeV in average
  double e_mu_lost = e_muon  - n_muon*1.3;
  double e_lost = e_neutr + e_mu_lost + e_to_tube;
  double e_real = evt_energy - e_lost;
  if(e_real< 0.0) e_real = 0.000001;

  /**  
  std::cout <<" =============================================================="<< std::endl;
  std::cout << " ========   Record Balance  ======="<< std::endl ;
  std::cout <<" =============================================================="<< std::endl;
  std::cout <<" ==============  Possible lost  ==================="<< std::endl;
  std::cout <<"  Neutrino energy      = "<<e_neutr<<",  in "<< n_neutr<<" neutrinos"<< std::endl;
  std::cout <<"  Energy to beam tube  = "<<e_to_tube<<",  in "<<n_to_tube<<" particles"<< std::endl;
  std::cout <<"  Muons energy lost    = "<<e_mu_lost<<"  in "<<n_muon<<" muons"<< std::endl;
  std::cout <<"  --------------------------------------------------"<< std::endl;
  std::cout <<"  Total Event energy at IP = "<<evt_energy<<" [GeV]"<< std::endl;
  std::cout <<"  --------------------------------------------------"<< std::endl;
  std::cout <<"  Whole lost Energy        = "<<e_lost<< std::endl;
  std::cout <<"  Available Energy in calo.= "<<e_real<< std::endl;
  std::cout <<" =============================================================="<< std::endl;
  std::cout <<"  Muon energy               = "<<e_muon  <<'\t'<<"  in "<< n_muon  <<" muons"<< std::endl;
  std::cout <<"  Electron energy           = "<<e_elect <<'\t'<<"  in "<< n_elect <<" electrons"<< std::endl;
  std::cout <<"  Charged hadron energy     = "<<e_chadr <<'\t'<<"  in "<< n_chadr <<" hadrons"<< std::endl;
  std::cout <<"  -------------------------------------------------------------"<< std::endl;
  std::cout <<"  Pi0 energy (if stable)    = "<<e_pi0   <<'\t'<<"  in "<< n_pi0   <<" Pi zeros"<< std::endl;
  std::cout <<"  Photon energy             = "<<e_photon<<'\t'<<"  in "<< n_photon<<" photons"<< std::endl;
  std::cout <<"  -------------------------------------------------------------"<< std::endl;
  std::cout <<"  Long lived hadron energy  = "<<e_llhadr<<'\t'<<"  in "<< n_llhadr<<" hadrons"<< std::endl;
  std::cout <<"  Short lived hadron energy = "<<e_slhadr<<'\t'<<"  in "<< n_slhadr<<" hadrons"<< std::endl;
  std::cout <<" =============================================================="<< std::endl;
  */

  return e_real;

 } // End  MC_Balance

}// namespace marlin
