#ifndef PFOID_h
#define PFOID_h 1

//#define root_out 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/ReconstructedParticle.h>

#include "PDF.h"

#include <vector>
#include <string>
#include <fstream>

#ifdef root_out
#include "TFile.h"
#include "TH1D.h"
#endif

using namespace lcio ;
using namespace marlin ;


/** === PFOID Processor === <br>
 *
 * What it does: For a set of reconstructed particle (by Wolf) it <br>
 *               continues to identify particles with the help of <br>
 *               the likelihood methode. The according histograms<br>
 *               (probability density functions PDF) must be available<br>
 *               in an ASCII file. <br>
 * What it needs: It needs a collection of reconstructed particles <br>
 *               created for example by Wolf and the PDF file. <br>
 * What' the Output: A new reconstracted particle collection. Each <br>
 *               particle has now an ParticleID object - i.e. a PDG. <br>
 *               The types, which are accessible via getType(), are now <br>
 *                      0  :  electron / positron (PDG: +-11) <br>
 *                      1  :  muon / antimuon (PDG: +-13) <br>
 *                      2  :  pion / antipion (PDG: +-211) <br>
 *                      3  :  gamma (PDG: 22) <br>
 *                      4  :  neutral Kaon (PDG: 130) <br>
 *              according the PDGs they have the corresponding masses. <br>
 * What can be adjusted: <br>
 *              "RecoParticleCollection" type="string" <br>
 *                                            default: RecoParticles <br>
 *                  from Wolf
 *              "NewRecoParticleCollection" type="string" <br>
 *                                            default: NewRecoParticles <br>
 *                  name of the new reconstructed particle collection <br>
 *              "FilePDFName" type="string" default: pdf.txt <br>
 *                  name of the ASCII file containing the PDFs <br>
 *                  for particles with track (charged). <br>
 *              "neutralFilePDFName" type="string" default: npdf.txt <br>
 *                  name of the ASCII file containing the PDFs <br>
 *                  for neutral particles. <br>
 *
 * written by Martin Ohlerich and Aliaksei Raspiareza (okt. 2006)
 */


class PDFs;

class PFOID : public Processor {

 public:

  virtual Processor* newProcessor() { return new PFOID ; }

  PFOID() ;

  virtual void init() ;

  virtual void processRunHeader( LCRunHeader * run ) ;

  virtual void processEvent( LCEvent * evt ) ;

  virtual void check( LCEvent * evt ) ;

  virtual void end() ;


 protected:

  int _nRun ;
  int _nEvt ;

  std::string _recoCol ;    // reconstructed particle Collection
  std::string _newrecoCol ;
  std::string _filename ;   // Name of file containing the pdfs (charged)
  std::string _filename1 ;  // Name of file containing the pdfs (neautral)
  float _bField;
  int noClusterParticle;

  struct info_t{
    double px, py, pz;   // momentum
    double Eecal, Ehcal; // subdetector energies
    double ex;           // excentricity
    double dmean;        // mean distance of hits to helix
    double Edmean;        // energy weighted mean distance of hits to helix
    double Necal, Nhcal;    // numbers of hits in subdetectors
    double L1,L2,L3;        // deepest Layers of HCal with hits
    double EtoN_ecal, EtoN_hcal, EecalToEtot; // ratios
    bool withTrack;   // flag wether it is 
    double EL1,EL2,EL3;        // first Layers of ECal with hits
  };
  info_t info;

  void init_info();
  void fill_info(int i, ReconstructedParticle *rp);

  PDF *pdf, *npdf;

#ifdef root_out
  TFile *file;
  TH1D *his_e;
  TH1D *his_mu;
  TH1D *his_pi;
  TH1D *his_g;
  TH1D *his_nh;
#endif
};

#endif


