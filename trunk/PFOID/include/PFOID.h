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

using namespace lcio ;
using namespace marlin ;


/** === PFOID Processor === <br>
 * Processor performs particle identification and discriminates between 
 * electrons, muons and charged hadrons in the case when calorimeter
 * cluster has an associated track and between photons and neutral hadrons 
 * in the case when no track is associated to cluster. 
 * Particle identification is based solely on the calorimeter information.
 * The likelihoods for different particle hypotheses are constructed from the 
 * variables which distinguish between various particle types. Highest likelihood
 * defines particle type. The following variables are used to construct likelihoods : <br>
 * - mean distance of hits to the associated track extrapolated into 
 * calorimeter volume (this variable is used only for charged particles); <br>
 * - ratio of energy deposited in ECAL to the total cluster energy; <br>
 * - cluster eccentricity (ratio of the cluster width to the cluster length); <br>
 * - average hit energy in ECAL; <br>
 * - average hit energy in HCAL; <br>
 * - first ECAL layer index in the calorimeter cluster 
 * (this variable is used only for neutral particles); <br>
 * - last HCAL layer index in the calorimeter cluster. <br>
 * For each particle reconstructed by particle flow algorithm, processor 
 * assigned ParticleID object, containing information about particle type.
 * The type of particle,  PDG code and likelihood can be accessed using
 * corresponding getter methods of the class ParticleID. Example below illustrates
 * how the information about particle ID can be accessed. <br>
 * ...
 * ReconstructedParticle * part = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
 * ParticleIDVec partIDVec = part->getParticleIDs();
 * ParticleID * partID = partIDVec[0];
 * float likelihood = partID->getLikelihood();
 * int type = partID->getType();
 * int pdg  = partID->getPDG();
 * ...
 * Types of the particles and their PDG codes are listed below:  <br>
 *   0  :  electron / positron (PDG: +-11) <br>
 *   1  :  muon- / muon+ (PDG: +-13) <br>
 *   2  :  pion+ / pion- (PDG: +-211) <br>
 *   3  :  gamma (PDG: 22) <br>
 *   4  :  neutral Kaon (PDG: 130) <br>
 * Note that every charged hadron is declared pion and every neutral   
 * hadron - long lived kaon.
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires collection of reconstructed particles. The name
 * of collection is specified via processor parameter "RecoParticleCollection".
 * Furthermore, a set of ASCII files, containing PDF's of variables entering likelihood
 * must be provided separately for charged and neutral particles for different 
 * energy bins. Ranges of energy bins are defined via processor parameter "EnergyBoundaries".
 * The names of ASCII files are provided by processor parameters
 * "FilePDFName" (for charged particles) and neutralFilePDFName (for neutral particles). <br>
 * <h4>Output</h4>
 * Processor assigns ParticleID objects to all reconstructed particles in event. <br>
 * @param RecoParticleCollection name of input collection of particles reconstructed by PFA (e.g. PandoraPFA) <br>
 * (default parameter value : "RecoParticles") <br>
 * @param FilePDFName names of ASCII files with pdf's of charged particles. Each file corresponds to certain energy bin. 
 * The order of files must be consistent with the energy bins. <br>
 * (default parameter value : "pdf_ILD00.txt") <br>
 * @param neutralFilePDFName names of ASCII files with pdf's of neutral particles. Each file corresponds to certain energy bin. 
 * The order of files must be consistent with the energy bins. <br>
 * (default parameter value : "npdf_ILD00.txt" . Only one energy bin, covering entire energy range, is assumed) <br> 
 * @param EnergyBoundaries the vector of energy boundaries, corresponding to different pdf files. The number of ranges 
 * should be consistent with the number of ASCII pdf files for charged and neutral particles (The number of boundaries equals
 * to the number of ASCII pdf files plus one). <br>
 * (default parameter values : 0.0 10000000.0 . Only one energy bin, covering the whole energy range, is assumed) <br>
 * @param Debug  debugging option <br>
 * (default parameter value : 0) <br>
 * <br>
 * @author M. Ohlerich (DESY Zeuthen) and A. Raspereza (MPI Munich)<br>
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
  int _debug;

  std::string _recoCol ;    // reconstructed particle Collection
  //  std::string _newrecoCol ;
  //  std::string _filename ;   // Name of file containing the pdfs (charged)
  //  std::string _filename1 ;  // Name of file containing the pdfs (neautral)
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
    bool withTrack;   // flag whether it is with track
    double EL1,EL2,EL3;        // first Layers of ECal with hits
  };
  info_t info;

  void init_info();
  void fill_info(int i, ReconstructedParticle *rp);

  std::vector<std::string> _filesCharged; // Names of files, containing pdfs (charged) 
  std::vector<std::string> _filesNeutral; // Names of files, containing pdfs (neutral)

  std::vector<PDF*> _pdfCharged; // set of PDFs for various energy bins (charged)
  std::vector<PDF*> _pdfNeutral; // set of PDFs for various energy bins (neutral)

  std::vector<float> _energyBoundaries; // energy boundaries defining binning in energy for PDFs

  int _nEnergyBins;

  PDF *pdf, *npdf;

};

#endif


