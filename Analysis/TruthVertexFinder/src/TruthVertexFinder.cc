#include "TruthVertexFinder.hh"
using std::string;
using std::vector;
using std::map;
namespace TTbarAnalysis
{
	TruthVertexFinder aTruthVertexFinder ;
	TruthVertexFinder::TruthVertexFinder() : Processor("TruthVertexFinder") 
	{

	    _description = "TruthVertexFinder extructs a secondary vertex from generator collection." ;


	    registerInputCollection( LCIO::MCPARTICLE,
        	    "CollectionName" , 
	            "Name of the MCParticle collection"  ,
        	    _colName ,
           	 std::string("MCParticleSkimmed")
	    );
	    registerOutputCollection( LCIO::VERTEX,
        	    "OutputCollectionName" , 
	            "Name of the Vertex collection"  ,
        	    _outputcolName ,
           	 std::string("MCVertex")
	    );
	    registerOutputCollection( LCIO::MCPARTICLE,
        	    "OutputBStarName" , 
	            "Name of the Vertex collection"  ,
        	    _outputBStarName ,
           	 std::string("BStar")
	    );
	    registerOutputCollection( LCIO::MCPARTICLE,
        	    "OutputProngsName" , 
	            "Name of the Prongs collection"  ,
        	    _outputProngsName ,
           	 std::string("EGProngs")
	    );
	    /*registerOutputCollection( LCIO::MCPARTICLE,
        	    "OutputBStarName" , 
	            "Name of the Vertex collection"  ,
        	    _outputBStarName ,
           	 std::string("BStar")
	    );*/
	    registerOutputCollection(LCIO::MCPARTICLE,
	    		"QuarkCollectionName" , 
	            "Name of the b-quark collection"  ,
        	    _outputquarkcolName,
           	 std::string("MCbquarks")
	    );
	    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
	    	"RelCollectionName",
		"Name of the Jet relation collection",
		_colRelName,
	    	std::string("RecoMCTruthLink")
	    );
	    registerProcessorParameter("DecayChainPDGs" , 
	            "PDGs of desired decays"  ,
        	    inputPdg,
           	 vector<int>()
	    );
	    _tagParameter = 6;
	    registerProcessorParameter("tagPDG" , 
	            "PDG of desired particle"  ,
        	    _tagParameter,
           	 _tagParameter
	    );
	    _writeROOTparameter = 0;
	    registerProcessorParameter("writeROOT" , 
	            "write ROOT file?"  ,
        	    _writeROOTparameter,
           	 _writeROOTparameter
	    );
	     _initialQuarkPDGparameter = 5;
	    registerProcessorParameter("initialQuarkPDG" , 
	            "PDG of initial particle"  ,
        	 _initialQuarkPDGparameter,
           	 _initialQuarkPDGparameter
	    );
	    _aParameter = 0.005;
	    registerProcessorParameter("a" , 
	            "a parameter of accuracy in mm"  ,
        	    _aParameter,
           	 _aParameter
	    );
	    _bParameter = 0.01;
	    registerProcessorParameter("b" , 
	            "b parameter of accuracy in mm"  ,
        	    _bParameter,
           	 _bParameter
	    );
	    _writeBonlyParameter = 1;
	    registerProcessorParameter("writeBonly" , 
	            "b parameter"  ,
        	    _writeBonlyParameter,
           	 _writeBonlyParameter
	    );
		_hfilename = "GenVertexTest.root";
	    registerProcessorParameter("ROOTFileName" , 
	            "ROOT File Name"  ,
        	    _hfilename,
           	 _hfilename
	    );
		//_pdgs.push_back(BOTTOM_HADRONS);
		//_pdgs.push_back(CHARMED_HADRONS);
		/*_pdgs.push_back(BOTTOM_MESONS);
		_pdgs.push_back(CHARMED_MESONS);
		_pdgs.push_back(EXCEPTIONAL_PDGTYPE);*/
	    ip[0] = 0.0;
	    ip[1] = 0.0;
	    ip[2] = 0.0;
	    ClearVariables();
	}	



	void TruthVertexFinder::init() 
	{ 
		streamlog_out(DEBUG) << "   init called  " << std::endl;
		printParameters() ;
		if (inputPdg.size() < 2) 
		{
			streamlog_out(ERROR) << " ERROR: size " << inputPdg.size() << " - Running on default settings" << std::endl;
			inputPdg.push_back(500);
			inputPdg.push_back(400);
			inputPdg.push_back(0);
		}
		for (unsigned int i = 0; i < inputPdg.size(); i++) 
		{
			_pdgs.push_back((PDGTYPE)inputPdg[i]);
		}
		_nRun = 0 ;
		_nEvt = 0 ;
		_hfilename = "TrashMCTest.root";
		if (_writeROOTparameter == 0) 
		{
			return;
		}
		_hfile = new TFile( _hfilename.c_str(), "RECREATE", _hfilename.c_str() ) ;
		_hBStarTree = new TTree( "BStar", "My test tree!" );
		_hBStarTree->Branch("bstarnumber",  &_bstarnumber, "bstarnumber/I");
		_hBStarTree->Branch("bstarmomentum",  _bstarmomentum, "bstarmomentum[bstarnumber]/F");
		_hBStarTree->Branch("bstaroffset",  _bstaroffset, "bstaroffset[bstarnumber]/F");
		_hVertexTree = new TTree( "Vertices", "My test tree!" );
		_hVertexTree->Branch("numberOfVertices", &_numberOfVertexes, "numberOfVertexes/I");
		_hVertexTree->Branch("distance", _distanceFromIP, "distance[numberOfVertexes]/F");
		_hVertexTree->Branch("coordinates", _coordinates, "coordinates[numberOfVertexes][3]/F");
		_hVertexTree->Branch("PDG", _PDG, "PDG[numberOfVertexes]/I");
		_hVertexTree->Branch("charge", _charge, "charge[numberOfVertexes]/I");
		_hVertexTree->Branch("generation", _generation, "generation[numberOfVertexes]/I");
		_hVertexTree->Branch("numberOfParticles", _numberOfParticles, "numberOfParticles[numberOfVertexes]/I");
		_hVertexTree->Branch("energyOfParticles", _energyOfParticles, "energyOfParticles[numberOfVertexes][15]/F");
		_hVertexTree->Branch("momentumOfParticles", _momentumOfParticles, "momentumOfParticles[numberOfVertexes][15]/F");
		_hVertexTree->Branch("massOfParticles", _massOfParticles, "massOfParticles[numberOfVertexes][15]/F");
		_hVertexTree->Branch("interactionOfParticles", _interactionOfParticles, "interactionOfParticles[numberOfVertexes][15]/I");
		if (_writeROOTparameter < 2) 
		{
			return;
		}
		_hTree = new TTree( "Stats", "My test tree!" );
		_hTree->Branch("tag", &_tag, "tag/I");
		_hTree->Branch("cosquark", &_cosquark, "cosquark/F");
		_hTree->Branch("cosantiquark", &_cosantiquark, "cosantiquark/F");
		_hTree->Branch("totalBcharge", &_totalBcharge, "totalBcharge/I");
		_hTree->Branch("ccharge", &_ccharge, "ccharge/I");
		_hTree->Branch("cbarcharge", &_cbarcharge, "cbarcharge/I");
		_hTree->Branch("bcharge", &_bcharge, "bcharge/I");
		_hTree->Branch("bbarcharge", &_bbarcharge, "bbarcharge/I");
		_hTree->Branch("baccuracy", &_baccuracy, "baccuracy/F");
		_hTree->Branch("bbaraccuracy", &_bbaraccuracy, "bbaraccuracy/F");
		_hTree->Branch("bIPdistance", &_bIPdistance, "bIPdistance/F");
		_hTree->Branch("bbarIPdistance", &_bbarIPdistance, "bbarIPdistance/F");
		_hTree->Branch("bdistance", &_bdistance, "bdistance/F");
		_hTree->Branch("bbardistance", &_bbardistance, "bbardistance/F");
		_hTree->Branch("bmomentum", &_bmomentum, "bmomentum/F");
		_hTree->Branch("bbarmomentum", &_bbarmomentum, "bbarmomentum/F");
		_hTree->Branch("cmomentum", &_cmomentum, "cmomentum/F");
		_hTree->Branch("cbarmomentum", &_cbarmomentum, "cbarmomentum/F");
		_hTree->Branch("caccuracy", &_caccuracy, "caccuracy/F");
		_hTree->Branch("cbaraccuracy", &_cbaraccuracy, "cbaraccuracy/F");
		_hTree->Branch("bnumber", &_bnumber, "bnumber/I");
		_hTree->Branch("bbarnumber", &_bbarnumber, "bbarnumber/I");
		_hTree->Branch("cnumber", &_cnumber, "cnumber/I");
		_hTree->Branch("cbarnumber", &_cbarnumber, "cbarnumber/I");
		_hTree->Branch("btotalnumber", &_btotalnumber, "btotalnumber/I");
		_hTree->Branch("bbartotalnumber", &_bbartotalnumber, "bbartotalnumber/I");
		_hTree->Branch("bptmiss", &_bptmiss, "bptmiss/F");
		_hTree->Branch("bbarptmiss", &_bbarptmiss, "bbarptmiss/F");
		_hTree->Branch("bnumber_f", &_bnumber_f, "bnumber_f/I");
		_hTree->Branch("bbarnumber_f", &_bbarnumber_f, "bbarnumber_f/I");
		_hTree->Branch("cnumber_f", &_cnumber_f, "cnumber_f/I");
		_hTree->Branch("cbarnumber_f", &_cbarnumber_f, "cbarnumber_f/I");
		//_hTree->Branch("firstVertexDistance", _firstVertexDistance, "firstVertexDistance[numberOfB0]/F");
		//_hTree->Branch("secondVertexDistance", _secondVertexDistance, "secondVertexDistance[numberOfB0]/F");
		//******************************************************************************************************//
		_hTrackTree = new TTree( "Tracks", "My test tree!" );
		_hTrackTree->Branch("bnumber", &_bnumber, "bnumber/I");
		_hTrackTree->Branch("bbarnumber", &_bbarnumber, "bbarnumber/I");
		_hTrackTree->Branch("cnumber", &_cnumber, "cnumber/I");
		_hTrackTree->Branch("cbarnumber", &_cbarnumber, "cbarnumber/I");

		_hTrackTree->Branch("bptrack", _bptrack, "bptrack[bnumber]/F");
		_hTrackTree->Branch("cptrack", _cptrack, "cptrack[cnumber]/F");
		_hTrackTree->Branch("bbarptrack", _bbarptrack, "bbarptrack[bbarnumber]/F");
		_hTrackTree->Branch("cbarptrack", _cbarptrack, "cbarptrack[cbarnumber]/F");

		_hTrackTree->Branch("betatrack", _betatrack, "betatrack[bnumber]/F");
		_hTrackTree->Branch("cetatrack", _cetatrack, "cetatrack[cnumber]/F");
		_hTrackTree->Branch("bbaretatrack", _bbaretatrack, "bbaretatrack[bbarnumber]/F");
		_hTrackTree->Branch("cbaretatrack", _cbaretatrack, "cbaretatrack[cbarnumber]/F");
		
		_hTrackTree->Branch("boffsettrack", _boffsettrack, "boffsettrack[bnumber]/F");
		_hTrackTree->Branch("coffsettrack", _coffsettrack, "coffsettrack[cnumber]/F");
		_hTrackTree->Branch("bbaroffsettrack", _bbaroffsettrack, "bbaroffsettrack[bbarnumber]/F");
		_hTrackTree->Branch("cbaroffsettrack", _cbaroffsettrack, "cbaroffsetttrack[cbarnumber]/F");

		_hMisRecoTree = new TTree( "Misreco", "My test tree!" );
		_hMisRecoTree->Branch("misreconumber",  &_misreconumber, "misreconumber/I");
		_hMisRecoTree->Branch("misrecotheta",  _misrecotheta, "misrecotheta[misreconumber]/F");
		_hMisRecoTree->Branch("misrecocostheta",  _misrecocostheta, "misrecocostheta[misreconumber]/F");
		_hMisRecoTree->Branch("misrecomomentum",  _misrecomomentum, "misrecomomentum[misreconumber]/F");
		_hMisRecoTree->Branch("misrecopt",  _misrecopt, "misrecopt[misreconumber]/F");

	}


	void TruthVertexFinder::processRunHeader( LCRunHeader* run) 
	{ 
		_nRun++ ;
	} 
	
	void TruthVertexFinder::PrintParticle(MCParticle * particle)
	{
		if (!particle) 
		{
			return;
		}
		streamlog_out(DEBUG) << std::fixed << std::setw( 6 ) << std::setprecision( 3 ) << std::setfill( ' ' );
		streamlog_out(DEBUG)<<"|"<<particle->getPDG() <<"\t\t|"<<particle->getMass()<<"\t\t|"<<particle->getCharge()  <<"\t\t|"<<particle->getEnergy()<<"\t\t|"<<particle->getVertex()[0]<<"\t\t|"<<particle->getVertex()[1]<<"\t\t|"<<particle->getVertex()[2] <<"\t\t|\n";
	
	}
	void TruthVertexFinder::Write(vector< Vertex * > * vertices, int & number)
	{
		if (!vertices || vertices->size() < _pdgs.size()-1) 
		{
			return;
		}
		for (unsigned int i = 0; i < vertices->size(); i++) 
		{
			Vertex * vertex = vertices->at(i);
			_PDG[_numberOfVertexes] = vertex->getParameters()[1];
			_generation[_numberOfVertexes] = vertex->getParameters()[2];
			ReconstructedParticle * particle = vertex->getAssociatedParticle();
			_charge[_numberOfVertexes] = particle->getCharge();
			_numberOfParticles[_numberOfVertexes] = particle->getParticles().size();
			_distanceFromIP[_numberOfVertexes] = vertex->getParameters()[0];
			MyVertex * myvertex = static_cast< MyVertex * >(vertex);
			for (unsigned int j = 0; j < particle->getParticles().size(); j++) 
			{
				_momentumOfParticles[_numberOfVertexes][j] = MathOperator::getModule(particle->getParticles()[j]->getMomentum());
				_massOfParticles[_numberOfVertexes][j] = particle->getParticles()[j]->getMass();
				_interactionOfParticles[_numberOfVertexes][j] = myvertex->__GetMCParticles()[j]->isDecayedInCalorimeter();
			}
			_numberOfVertexes++;
		}
	}
	void TruthVertexFinder::PrintChain(vector< MCParticle * > * chain)
	{
		if (!chain) 
		{
			return;
		}
		for (unsigned int i = 0; i < chain->size(); i++) 
		{
			PrintParticle(chain->at(i));
		}
	}
	void TruthVertexFinder::WriteQuarksCollection(LCEvent * evt, std::vector< MCParticle * > & quarks)
	{
		IMPL::LCCollectionVec * mc = new IMPL::LCCollectionVec ( LCIO::MCPARTICLE ) ;
		if (quarks.size() == 2) 
		{
			//streamlog_out(DEBUG) << "PDG: " << quarks[0]->getPDG() << '\n';
			//streamlog_out(DEBUG) << "PDG: " << quarks[1]->getPDG() << '\n';
			mc->addElement(quarks[0]);
			mc->addElement(quarks[1]);
		}
		evt->addCollection( mc , _outputquarkcolName ) ;
	}

	void TruthVertexFinder::AddProngs( VertexMCOperator & vertexOperator, MCOperator & opera, DecayChain * chain, vector< Vertex * > * verticies, std::vector<int> & parameters, IMPL::LCCollectionVec * col)
	{
		if (!verticies || !chain) 
		{
			return;
		}
		int vsize = verticies->size();
		bool useRelation = false;
		if (vsize > 1) 
		{
			vector< MCParticle * > bdaughters = opera.SelectStableCloseDaughters(chain->Get(0), chain->Get(1)->getPDG());
			vertexOperator.AddProngs(verticies->at(0), bdaughters, useRelation);
			vector< MCParticle * > cdaughters = opera.SelectStableCloseDaughters(chain->Get(1));
			vertexOperator.AddProngs(verticies->at(1), cdaughters, useRelation);
			if (col) 
			{
				for (unsigned int i = 0; i < bdaughters.size(); i++) 
				{
					parameters.push_back( (float) chain->GetParentPDG() / 5.0 * 2);
					col->addElement(bdaughters[i]);
				}
				for (unsigned int i = 0; i < cdaughters.size(); i++) 
				{
					parameters.push_back( (float) chain->GetParentPDG() / 5.0 * 3);
					col->addElement(cdaughters[i]);
				}
			}
		}
		else 
		{
			PrintParticle(chain->Get(0));
			vector< MCParticle * > bdaughters = opera.SelectStableCloseDaughters(chain->Get(0));
			streamlog_out(DEBUG)<< "Prongs for quark " << chain->GetParentPDG() << ": \n";
			for (unsigned int i = 0; i < bdaughters.size(); i++) 
			{
				 PrintParticle(bdaughters[i]);
			}
			vertexOperator.AddProngs(verticies->at(0), bdaughters, useRelation);
			if (col) 
			{
				for (unsigned int i = 0; i < bdaughters.size(); i++) 
				{
					col->addElement(bdaughters[i]);
				}
			}
			
		}
	}

	void TruthVertexFinder::processEvent( LCEvent * evt ) 
	{ 
		try
		{
			LCCollection* col = evt->getCollection( _colName );
			streamlog_out(DEBUG)<< "***********TruthVertexFinder*"<<_nEvt<<"***************\n";
			LCCollection* rel = evt->getCollection(_colRelName);
			MCOperator opera(col,rel);
			VertexMCOperator vertexOperator(rel);
			
			_tag = (opera.CheckProcessForPair(_tagParameter))? 1 : 0;
			vector< MCParticle * > bquarks = opera.GetPairParticles(_pdgs[0]);
	 		//GetAsymmetry(bquarks);
			_nEvt ++ ;
			int initQuark = abs(_initialQuarkPDGparameter);
			streamlog_out(DEBUG)<<"\t|PDG\t\t|Mass\t\t|Charge\t\t|Energy\t\t|Vtx X\t\t|Vtx Y\t\t|Vtx Z\t\t|\n";
			DecayChain * bChainRaw = opera.Construct(string("b-quark decay chain"), initQuark, _pdgs);
			DecayChain * bChain = opera.RefineDecayChain(bChainRaw, _pdgs);
			IMPL::LCCollectionVec * mc = new IMPL::LCCollectionVec ( LCIO::MCPARTICLE ) ;
			mc->setSubset();
			if (bChain) 
			{
				vector< MCParticle * > daughters = opera.SelectStableCloseDaughters(bChainRaw->Get(0), bChain->Get(0)->getPDG());	
				streamlog_out(DEBUG)<<"Additional B particles: \n";
				for (unsigned int i = 0; i < daughters.size(); i++) 
				{
					vector< float > direction = MathOperator::getDirection(daughters[i]->getMomentum());
					_bstaroffset[_bstarnumber] = MathOperator::getDistanceTo(ip, direction, bChain->Get(1)->getVertex());
					_bstarmomentum[_bstarnumber++] = MathOperator::getModule(daughters[i]->getMomentum());
					mc->addElement(daughters[i]);
					PrintParticle(daughters[i]);
				}
			}
			DecayChain * bbarChainRaw = opera.Construct(string("bbar-quark decay chain"), 0-initQuark, _pdgs);
			DecayChain * bbarChain = opera.RefineDecayChain(bbarChainRaw, _pdgs);
			if (bbarChain) 
			{
				vector< MCParticle * > daughters = opera.SelectStableCloseDaughters(bbarChainRaw->Get(0), bbarChain->Get(0)->getPDG());	
				streamlog_out(DEBUG)<<"Additional Bbar particles: \n";
				for (unsigned int i = 0; i < daughters.size(); i++) 
				{
					vector< float > direction = MathOperator::getDirection(daughters[i]->getMomentum());
					_bstaroffset[_bstarnumber] = MathOperator::getDistanceTo(ip, direction, bbarChain->Get(1)->getVertex());
					_bstarmomentum[_bstarnumber++] = MathOperator::getModule(daughters[i]->getMomentum());
					mc->addElement(daughters[i]);
					PrintParticle(daughters[i]);
				}
			}
			evt->addCollection( mc , _outputBStarName ) ;

			vector< Vertex * > * bverticies = vertexOperator.Construct(bChain);
			vector< Vertex * > * bbarverticies = vertexOperator.Construct(bbarChain);
			IMPL::LCCollectionVec * prongs = new IMPL::LCCollectionVec ( LCIO::MCPARTICLE ) ;
			vector <int> parameters;
			prongs->setSubset ();
			//prongs->setFlag	(16);
			AddProngs(vertexOperator, opera, bChain, bverticies,parameters, prongs);
			AddProngs(vertexOperator, opera, bbarChain, bbarverticies,parameters, prongs);
			
			
			
			WriteQuarksCollection(evt, bquarks);
			WriteVertexCollection(evt, bverticies, bbarverticies);
			
			evt->addCollection( prongs , _outputProngsName);
			prongs->parameters().setValues("trackIDs", parameters);
			prongs->parameters().setValue("GenVertexTag", _tag);
			//int number = 0;
			if (_writeROOTparameter > 0) 
			{
				_numberOfVertexes = 0;
				Write(bverticies,_numberOfVertexes);
				Write(bbarverticies, _numberOfVertexes);
				if (_writeROOTparameter > 1) 
				{
					Write(opera, bChain,bverticies);
					Write(opera, bbarChain,bbarverticies);
					_hTree->Fill();
					_bnumber = (_bnumber < 0)? 0: _bnumber;
					_bbarnumber = (_bbarnumber < 0)? 0: _bbarnumber;
					_cnumber = (_cnumber < 0)? 0: _cnumber;
					_cbarnumber = (_cbarnumber < 0)? 0: _cbarnumber;
					_hTrackTree->Fill();
						
				}
				_hBStarTree->Fill();
				_hVertexTree->Fill();
				ClearVariables();
			}
			streamlog_out(DEBUG)<<"B cos: " << _cosquark << '\n';
			streamlog_out(DEBUG)<<"Bbar cos: " <<_cosantiquark << '\n';
	
		}
		catch( DataNotAvailableException &e)
		{
			streamlog_out(DEBUG) << "No collection!" << std::endl ;
		}
	}
	void TruthVertexFinder::GetAsymmetry(std::vector< MCParticle * > & particles)
	{
		/*if (!particles) 
		{
			return;
		}*/
		if (particles.size() < 2) 
		{
			return;
		}
		vector< float > direction = MathOperator::getDirection(particles.at(0)->getEndpoint());
		vector< float > antidirection = MathOperator::getDirection(particles.at(1)->getEndpoint());
		_cosquark = std::cos(MathOperator::getAngles(direction)[1]);
		_cosantiquark = std::cos(MathOperator::getAngles(antidirection)[1]);
		//delete particles;
	}
	void TruthVertexFinder::Write(MCOperator & opera, DecayChain * chain, vector< Vertex * > * verticies)
	{
		if (!chain || !chain->Get(0)) 
		{
			return;
		}
		for (int i = 0; i < chain->GetSize(); i++) 
		{
			PrintParticle(chain->Get(i));
		}
		//vector< MCParticle * > * misreco = new vector< MCParticle * >();
		vector< float > direction = MathOperator::getDirection(chain->Get(0)->getMomentum());
		if (chain->GetParentPDG() > 0 && chain->GetSize() > 1) 
		{
			
			const vector< MCParticle * > daughters = static_cast< MyVertex * >(verticies->at(0))->__GetMCParticles(); // opera.SelectStableCloseDaughters(chain->Get(0), chain->Get(1)->getPDG()); //opera.ScanForVertexParticles(bverticies->at(0)->getPosition(), 1e-3);
			_bcharge = (int)chain->Get(0)->getCharge();
			_bnumber = daughters.size();
			_bnumber_f = opera.SelectStableCloseDaughters(chain->Get(0), chain->Get(1)->getPDG(),true).size();//opera.SelectStableCloseDaughters(chain->Get(0), chain->Get(1)->getPDG(),true).size();//opera.CheckDaughterVisibility(daughters).size();
			_bIPdistance = verticies->at(0)->getParameters()[0];
			Write(daughters, 1);
			streamlog_out(DEBUG)<<"Vertex b-quark: " << verticies->at(0)->getParameters()[0]<< " n-tracks: " << _bnumber << '\n';

			const vector< MCParticle * > cdaughters =static_cast< MyVertex * >(verticies->at(1))->__GetMCParticles(); //opera.SelectStableCloseDaughters(chain->Get(1)); //opera.ScanForVertexParticles(bverticies->at(1)->getPosition(), 1e-3);
			_cnumber = cdaughters.size();
			_cnumber_f = opera.SelectStableCloseDaughters(chain->Get(1),0,true).size();//opera.CheckDaughterVisibility(cdaughters).size();
			
			opera.CheckDaughterVisibility(daughters);
			Write(cdaughters, 2);
			
			streamlog_out(DEBUG)<<"Vertex c-quark: " << verticies->at(1)->getParameters()[0] <<" n-tracks: " << _cnumber <<  '\n';
			opera.CheckDaughterVisibility(cdaughters);
			_bdistance = MathOperator::getDistance(verticies->at(1)->getPosition(), verticies->at(0)->getPosition());
			_btotalnumber = _cnumber + _bnumber;
			streamlog_out(DEBUG)<<"Checking b-quark meson...\n";
			bool compatible = opera.CheckCompatibility(daughters, chain->Get(0), chain->Get(1)->getCharge());
				
			streamlog_out(DEBUG)<<"Checking c-quark meson...\n";
			compatible = opera.CheckCompatibility(cdaughters, chain->Get(1));
			_bptmiss = getMissingPt(daughters, cdaughters, verticies->at(0));
			streamlog_out(DEBUG)<<"Missing pt for b-quark hadron: " << _bptmiss << "\n";	
			_ccharge = (int)chain->Get(1)->getCharge();
			_bmomentum = MathOperator::getModule(chain->Get(0)->getMomentum());
			_baccuracy = opera.GetAccuracy(chain->Get(0), _aParameter, _bParameter); 
			_cmomentum = MathOperator::getModule(chain->Get(1)->getMomentum());
			_caccuracy = opera.GetAccuracy(chain->Get(1), _aParameter, _bParameter);
			_cosquark = std::cos(MathOperator::getAngles(direction)[1]);
		}
		if (chain->GetParentPDG() < 0 && chain->GetSize() > 1) 
		{
			const vector< MCParticle * > daughters = static_cast< MyVertex * >(verticies->at(0))->__GetMCParticles();// opera.SelectStableCloseDaughters(chain->Get(0), chain->Get(1)->getPDG()); // opera.ScanForVertexParticles(bbarverticies->at(0)->getPosition(), 1e-3);
			Write(daughters, -1);
			_bbarcharge = (int) chain->Get(0)->getCharge();
			_bbarnumber = daughters.size();
			_bbarIPdistance = verticies->at(0)->getParameters()[0];
			_bbarnumber_f = opera.SelectStableCloseDaughters(chain->Get(0), chain->Get(1)->getPDG(),true).size();//opera.CheckDaughterVisibility(daughters).size();
		        streamlog_out(DEBUG)<<"Vertex bbar-quark"<< 0 <<": " << verticies->at(0)->getParameters()[0] <<" n-tracks: " << _bbarnumber <<  '\n';
			const vector< MCParticle * > cdaughters= static_cast< MyVertex * >(verticies->at(1))->__GetMCParticles(); //opera.SelectStableCloseDaughters(chain->Get(1)); //opera.ScanForVertexParticles(bbarverticies->at(1)->getPosition(), 1e-3);
			_cbarnumber = cdaughters.size();
			opera.CheckDaughterVisibility(daughters);
			streamlog_out(DEBUG)<<"Vertex cbar-quark: " << verticies->at(1)->getParameters()[0] <<" n-tracks: " <<_cbarnumber << '\n';
			opera.CheckDaughterVisibility(cdaughters);
			_cbarnumber_f = opera.SelectStableCloseDaughters(chain->Get(1),0,true).size();//opera.CheckDaughterVisibility(cdaughters).size();
			_bbardistance = MathOperator::getDistance(verticies->at(1)->getPosition(), verticies->at(0)->getPosition());
			_bbartotalnumber = _cbarnumber + _bbarnumber;
			_cbarcharge = (int) chain->Get(1)->getCharge();
			Write(cdaughters, -2);
			_bbarptmiss = getMissingPt(daughters, cdaughters, verticies->at(0));
			streamlog_out(DEBUG)<<"Missing pt for bbar-quark hadron: " << _bbarptmiss << "\n";	
			_bbarmomentum = MathOperator::getModule(chain->Get(0)->getMomentum());
			_bbaraccuracy = opera.GetAccuracy(chain->Get(0), _aParameter, _bParameter); 
			_cbarmomentum = MathOperator::getModule(chain->Get(1)->getMomentum());
			_cbaraccuracy = opera.GetAccuracy(chain->Get(1), _aParameter, _bParameter); 
			_cosantiquark = std::cos(MathOperator::getAngles(direction)[1]);

		}
		//WriteMisReco(misreco);
	}
	void TruthVertexFinder::WriteMisReco(vector< MCParticle * > * particles)
	{

		streamlog_out(DEBUG) << "Misreco: " << _misreconumber << " particles: " << particles->size() << '\n';
		for (unsigned int i = _misreconumber; i < _misreconumber + particles->size(); i++) 
		{
			vector< float > direction = MathOperator::getDirection(particles->at(i-_misreconumber)->getMomentum());
			_misrecotheta[i] = MathOperator::getAngles(direction)[1];
			_misrecocostheta[i] = std::cos(MathOperator::getAngles(direction)[1]);
			_misrecomomentum[i] = MathOperator::getModule(particles->at(i-_misreconumber)->getMomentum());
			_misrecopt[i] = MathOperator::getPt(particles->at(i-_misreconumber)->getMomentum());
		}
		_misreconumber += particles->size();
	}
	double TruthVertexFinder::getMissingPt(const vector< MCParticle * > & bdaugthers, const vector< MCParticle * > & cdaughters, Vertex * vertex)
	{
		const float * position = vertex->getPosition();
		for (int i = 0; i < 3; i++) 
		{
			streamlog_out(DEBUG) << i << ": " << position[i];
		}
		streamlog_out(DEBUG) << '\n';
		vector< const double * > vectors;
		for (unsigned int i = 0; i < bdaugthers.size(); i++) 
		{
			vectors.push_back(bdaugthers[i]->getMomentum());
		}
		for (unsigned int i = 0; i < cdaughters.size(); i++) 
		{
			vectors.push_back(cdaughters[i]->getMomentum());
		}
		double missing =  MathOperator::getMissingPt(vectors, position);
		return missing;
	}
	void TruthVertexFinder::Write(const vector< MCParticle * > daughters, int v)
	{
		float * offset = NULL;
		float * pt = NULL;
		float * p = NULL;
		float * eta = NULL;
		switch(v)
		{
			case 1:
			offset = _boffsettrack;
			p = _bptrack;
			eta = _betatrack;
			break;
			case 2:
			offset = _coffsettrack;
			p = _cptrack;
			eta = _cetatrack;
			break;
			case -2:
			offset = _cbaroffsettrack;
			p = _cbarptrack;
			eta = _cbaretatrack;
			break;
			case -1:
			offset = _bbaroffsettrack;
			p = _bbarptrack;
			eta = _bbaretatrack;
			break;

		}
		int size = daughters.size();
		for (int i = 0; i < size; i++) 
		{
			MCParticle * daughter = daughters[i];
			vector< float > direction = MathOperator::getDirection(daughter->getMomentum());
			offset[i] = MathOperator::getDistanceTo(ip, direction, daughter->getVertex());
			//float pt[i] = MathOperator::getPt(daughter->getMomentum());
			p[i] = MathOperator::getModule(daughter->getMomentum());
			eta[i] = MathOperator::getAngles(direction)[1];
		}
	}
	void TruthVertexFinder::WriteVertexCollection(LCEvent * evt, vector< Vertex * > * bvertexes, vector< Vertex * > * bbarvertexes)
	{
		IMPL::LCCollectionVec * mc = new IMPL::LCCollectionVec ( EVENT::LCIO::VERTEX ) ;
		if (bvertexes) 
		{
			for (unsigned int i = 0; i < bvertexes->size(); i++) 
			{
				
				mc->addElement(bvertexes->at(i));
			}
		}
		if (bbarvertexes) 
		{
			for (unsigned int i = 0; i < bbarvertexes->size(); i++)
			{
			        mc->addElement(bbarvertexes->at(i));
			}
		}
		mc->parameters().setValue("GenVertexTag", _tag);
		evt->addCollection( mc , _outputcolName ) ;
	}

	void TruthVertexFinder::ClearVariables()
	{
		_misreconumber = 0;
		_tag = false;
		_bstarnumber = 0;
	  	_bptmiss = -1.0;
	  	_bbarptmiss = -1.0;
		_cosquark = -2.0;
		_cosantiquark = -2.0;
		_totalBcharge = -3.0;
		_ccharge = -3.0;
		_cbarcharge = -3.0;
		_bcharge = -3.0;
		_bbarcharge = -3.0;
		_bdistance = -1.0;
		_bbardistance = -1.0;
		_bmomentum = -1.0;
		_bbarmomentum = -1.0;
		_bbaraccuracy = -1.0;
		_baccuracy = -1.0;
		_cmomentum = -1.0;
		_cbarmomentum = -1.0;
		_cbaraccuracy = -1.0;
		_caccuracy = -1.0;
		_bnumber = -1;
		_bbarnumber = -1;
		_cnumber = -1;
		_cbarnumber = -1;
		_bbartotalnumber = -1;
		_btotalnumber = -1;
		_bnumber_f = -1;
		_bbarnumber_f = -1;
		_cnumber_f = -1;
		_cbarnumber_f = -1;
		for (int i = 0; i < 2; i++) 
		{
			_firstVertexDistance[i] = 0.0;
			_secondVertexDistance[i] = 0.0;
			_numberOfB0 = 0;
		}
		for (int i = 0; i < MAXV; i++) 
		{
			_bstaroffset[i]  = -1.0;
			_bstarmomentum[i] = -1.0;
			_bptrack[i] = -1.0;
			_betatrack[i] = -1.0;
			_boffsettrack[i] = -1.0;
			_bbarptrack[i] = -1.0;
			_bbaretatrack[i] = -1.0;
			_bbaroffsettrack[i] = -1.0;
			_cptrack[i] = -1.0;
			_cetatrack[i] = -1.0;
			_coffsettrack[i] = -1.0;
			_cbarptrack[i] = -1.0;
			_cbaretatrack[i] = -1.0;
			_cbaroffsettrack[i]  = -1.0;
			_numberOfParticles[i] = -1;
			_charge[i] = 0;
			_PDG[i] = 0;
			_generation[i] = 0;
			for (int j = 0; j < MAXV; j++) 
			{
				_energyOfParticles[i][j] = -1.0;
				_momentumOfParticles[i][j] = -1.0;
				_massOfParticles[i][j] = -1.0;
				_interactionOfParticles[i][j] = -1;
			}
		}
	}
	
	
	void TruthVertexFinder::check( LCEvent * evt ) 
	{ 
		// nothing to check here - could be used to fill checkplots in reconstruction processor
		
	}
	
	
	void TruthVertexFinder::end()
	{ 
		if (_writeROOTparameter == 0) 
		{
			return;
		}
		_hfile->cd();
		_hfile->Write();
		_hfile->Close();
	
	    //   streamlog_out(DEBUG) << "TruthVertexFinder::end()  " << name() 
	    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    // 	    << std::endl ;
	
	}
} /* TTbarAnalysis */
