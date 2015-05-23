#ifndef TruthVertexFinder_h
#define TruthVertexFinder_h 1
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Vertex.h>
#include "marlin/VerbosityLevels.h"
#include "marlin/Processor.h"
#include "lcio.h"

// ----- include for verbosity dependend logging ---------
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <map>

#include "ConstantStorage.hh"
#include "MathOperator.hh"
#include "MCOperator.hh"
#include "VertexMCOperator.hh"
#include "MyVertex.hh"
using namespace lcio ;
using namespace marlin ;


namespace TTbarAnalysis 
{
	////////////////////////////////////////////////////////
	/// This processor is designed to extruct a secondary///
	/// vertices from collection of generated particles, ///
	/// by using a PDG type and daughter/parent relation ///
	/// relations. 					     ///
	/// TruthVertexFinder has two main outputs:	     ///
	/// EGProngs - MCParticle type with parameters	     ///
	/// MCVertex - Vertex type with particles	     ///
	/// For more info and usage please see doc/ folder   ///
	////////////////////////////////////////////////////////
	/// Author: BILOKIN Sviatoslav, PhD student	     ///
	///	    POESCHL Roman, Supervisor		     ///
	///	    RICHARD Francois, Supervisor	     ///
	/// 				designed: 2015-2017  ///
	////////////////////////////////////////////////////////
	class TruthVertexFinder : public Processor 
	{
	  
	 public:
	  
	  virtual Processor*  newProcessor() { return new TruthVertexFinder ; }
	  
	  
	  TruthVertexFinder() ;
	  
	  virtual void init() ;
	  virtual void processRunHeader( LCRunHeader* run ) ;
	  virtual void processEvent( LCEvent * evt ) ; 
	  virtual void check( LCEvent * evt ) ; 
	  virtual void end() ;
	 /////////////////CUSTOM/////////////////////////// 
	  void PrintParticle(MCParticle * particle);
	  void PrintChain(std::vector< MCParticle * > * chain);

	  void WriteVertexCollection(LCEvent * evt, std::vector< Vertex * > * bvertexes, std::vector< Vertex * > * bbarvertexes);
	  void Write(std::vector< EVENT::Vertex * > * vertices, int & number);
	  void Write(MCOperator & opera,DecayChain * chain, std::vector< Vertex * > * bvertexes);
	  void AddProngs( VertexMCOperator & vertexOperator, MCOperator & opera, DecayChain * chain, std::vector< Vertex * > * vertices, std::vector<int> & parameters, IMPL::LCCollectionVec * col = NULL);
	  void Write(const std::vector< MCParticle * > particle , int v);
	  double getMissingPt(const std::vector< MCParticle * > & bdaugthers, const  std::vector< MCParticle * > & cdaughters, Vertex * vertex);
	  void WriteQuarksCollection(LCEvent * evt, std::vector< MCParticle * > & quarks);
	  void WriteMisReco(std::vector< MCParticle * > * particles);
	  void GetAsymmetry(std::vector< MCParticle * > & particles);
	  void ClearVariables(); 
	 protected:
	
	  // Input/output collection names
	  
	  std::string _colName ;
	  std::string _outputcolName;
	  std::string _outputquarkcolName;
	  std::string _outputBStarName; 
	  std::string _outputProngsName; 
	  std::string _colRelName;
	  std::string _outputBbarStarName; 

	  // Parameters

	  std::vector<PDGTYPE> _pdgs;
	  IntVec inputPdg;
	  int _tagParameter;
	  int _initialQuarkPDGparameter;
	  float _aParameter;
	  float _bParameter;
	  int _writeBonlyParameter;
	  int _writeROOTparameter;

	  // Root variables
	  TFile * _hfile;
	  TTree * _hTree;
	  TTree * _hVertexTree;
	  TTree * _hTrackTree;
	  TTree * _hBStarTree;
	  TTree * _hMisRecoTree;
	  std::string _hfilename ;
	  
	  // Inner variables
	  int _tag;
	  int _numberOfB0;
	  float _firstVertexDistance[2];
	  float _secondVertexDistance[2];
	  int _totalBcharge;
	  int _ccharge;
	  int _cbarcharge;
	  int _bcharge;
	  int _bbarcharge;
	  float _baccuracy;
	  float _bbaraccuracy;
	  float _bIPdistance;
	  float _bbarIPdistance;
	  float _btracks;
	  float _bbartracks;
	  float _ctracks;
	  float _cbartracks;
	  float _bdistance;
	  float _bbardistance;
	  float _bmomentum;
	  float _bbarmomentum;
	  float _cmomentum;
	  float _cbarmomentum;
	  float _caccuracy;
	  float _cbaraccuracy;
	  int _bnumber;
	  int _bbarnumber;
	  int _cnumber;
	  int _cbarnumber;
	  int _btotalnumber;
	  int _bbartotalnumber;
	  int _bnumber_f;
	  int _bbarnumber_f;
	  int _cnumber_f;
	  int _cbarnumber_f;
	  float _bptmiss;
	  float _bbarptmiss;
	  float _cosquark;
	  float _cosantiquark;

	
	  static const int MAXV = 15;
	  int _numberOfVertexes;
	  float _distanceFromIP[MAXV];
	  float _coordinates[MAXV][3];
	  int _PDG[MAXV];
	  int _generation[MAXV];
	  int _charge[MAXV];
	  int _numberOfParticles[MAXV];
	  float _energyOfParticles[MAXV][MAXV];
	  float _momentumOfParticles[MAXV][MAXV];
	  float _massOfParticles[MAXV][MAXV];
	  int _interactionOfParticles[MAXV][MAXV];
	
	  float _bptrack[MAXV];
	  float _betatrack[MAXV];
	  float _boffsettrack[MAXV];
	  float _bbarptrack[MAXV];
	  float _bbaretatrack[MAXV];
	  float _bbaroffsettrack[MAXV];
	  float _cptrack[MAXV];
	  float _cetatrack[MAXV];
	  float _coffsettrack[MAXV];
	  float _cbarptrack[MAXV];
	  float _cbaretatrack[MAXV];
	  float _cbaroffsettrack[MAXV];

	  static const int MAXVV = 30;
	  int _misreconumber;
	  float _misrecotheta[MAXVV];
	  float _misrecocostheta[MAXVV];
	  float _misrecomomentum[MAXVV];
	  float _misrecopt[MAXVV];

	  int _bstarnumber;
	  float _bstarmomentum[MAXV];
	  float _bstaroffset[MAXV];
	  double ip[3];
	  int _nRun ;
	  int _nEvt ;
	} ;
} /* TTbarAnalysis */
#endif



