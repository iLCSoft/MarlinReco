#include "KIT.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;
using namespace IMPL;
using namespace EVENT;
using namespace UTIL;
                         	
KIT aKIT ;
KIT::KIT() : Processor("KIT") {
  
  // modify processor description
  _description = "KIT does whatever it does ..." ;
  
  registerProcessorParameter( "ECAL_Collection",
			      "Name of the CalorimeterHit Collection for ECAL ",
			      _Ecal_col,
			      std::string("ECAL"));

  registerProcessorParameter( "Core_Collection",
			      "Name of the Cluster Collection for Cores ",
			      _CoreCollection,
			      std::string("CORE"));

  registerProcessorParameter( "Cleaning",
			      "To do the cleaning on hits or not ",
			      _ToClean,
			      std::string("YES"));
  registerProcessorParameter( "TopologicalCut",
			      "At which number of neighbors to put the threshold, condition is < so you need to put N+1 ",
			       _CleanCut,
			      (int)5);

  registerProcessorParameter( "NumberOfLevels",
			      "Number of levels for central loop ",
			       _N,
			      (int)10);
  vector<float> miipstep;
  miipstep.push_back(0.1);
  miipstep.push_back(1.5);
  miipstep.push_back(2.5);
  miipstep.push_back(4.0);
  miipstep.push_back(6.0);
  miipstep.push_back(9.0);
  miipstep.push_back(16.0);
  miipstep.push_back(26.0);
  miipstep.push_back(41.0);
  miipstep.push_back(65.0);

  registerProcessorParameter( "Levels",
			     "Levels for central loop in MIP ",
			      _miipstep,
			       miipstep);

  registerProcessorParameter( "MinHit0",
			     "Minimal Number of hits for ground level cluster ",
			      _MinHit0,
			      (int)4);

  registerProcessorParameter( "MinHitSplit",
			      "Minimal Number of hits for i-th level cluster ",
			      _MinHitSplit,
			      (int)2);

  registerProcessorParameter( "Rcut",
			      "Fluctuation suprresion cut",
			      _Rcut,
			      (double)0.4);
  registerProcessorParameter( "Distcut",
			      "Square of distance cut for merging ",
			      _Distcut,
			      (double)35.0);
  registerProcessorParameter( "Coscut",
			      "Cosine of the angle for merging ",
			      _Coscut,
			      (double)0.95);
 
}


void KIT::init() { 
     
   printParameters() ;
  _nRun = 0 ;
  _nEvt = 0 ;   
}

void KIT::processRunHeader( LCRunHeader*  /*run*/) { 
      
  _nRun++ ;
} 

void KIT::processEvent( LCEvent * evt ) {  
 
  try{   
        LCCollection* colt = evt->getCollection(_Ecal_col.c_str()) ;
	CellIDDecoder<CalorimeterHit> CDECAL(colt);
	LCCollectionVec * clscol = new LCCollectionVec(LCIO::CLUSTER);
      if( colt!=0)
	{
	  unsigned int nelem=colt->getNumberOfElements(); 
	  // MAIN CONTAINER OF SHITS 
	  vector<Superhit2*> calo[10]; 
                        
	  // creating all superhits 
	  CreateAllShits2(colt,CDECAL,calo); 

	  // precalculation
	  TotalPrecalc2(calo,CDECAL,nelem,_CleanCut);  

	  // setting the parameters of the alghorithm
	  vector <PROTSEED2> prs2;
	  CoreCut2 Ccut;
	  Ccut.Rcut=_Rcut;
	  Ccut.Distcut=_Distcut;
	  Ccut.Coscut=_Coscut;
	  Ccut.MinHit0=(unsigned int) _MinHit0;
	  Ccut.MinHitSplit=(unsigned int) _MinHitSplit;
          const unsigned int N=_N;
	  vector< vector<Tmpcl2*> > bbb(N);
	 
	  // finding cores.
	  if( _ToClean=="YES" || _ToClean=="yes")
	    {
	      // cout << " da koristim cat " << endl;
	       FindCores2(&(calo[4]), bbb , &prs2,_N,_miipstep,Ccut);     
	    }else{
	       FindCores2(&(calo[0]), bbb , &prs2,_N,_miipstep,Ccut);   
	    }


	  // writting out in LCIO
	 
	  for(unsigned int i=0;i<prs2.size();i++)
	    {
	      if(prs2[i].active==true)
		{
		
		  ClusterImpl * cluster = new ClusterImpl();
		  
		  for(unsigned int j=0;j<prs2[i].cl->hits.size();j++)
		    {
		      cluster->addHit( prs2[i].cl->hits[j]->chit,(float)1.0);
		    }		
		  
		  cluster->setTypeBit(prs2[i].level);
	       
		  float position[3];
		  position[0]=(float)prs2[i].cl->getCenter()[0];
		  position[1]=(float)prs2[i].cl->getCenter()[1];
		  position[2]=(float)prs2[i].cl->getCenter()[2];

		  float energy=(float)prs2[i].cl->getEnergy();

		  cluster->setPosition(position); 
		  cluster->setEnergy(energy);
		  clscol->addElement(cluster);
	        }
	    }
	

	  // for strong memory and nice dreams ..
	  for(unsigned int i=0;i<N;i++)
	    {
	      if( bbb[i].size()!=0)
		for(unsigned int im=0;im<bbb[i].size();im++)
		  {
		    delete bbb[i][im];
		  }
	    }

	  for(unsigned int im=0;im<2;im++)
	    {        
	      if(calo[im].size()!=0)
		for( unsigned int iij=0;iij<calo[im].size();iij++)
		  delete (calo[im])[iij];
	    }

    }      
    
   evt->addCollection(clscol,_CoreCollection.c_str());

  }catch(DataNotAvailableException &e) {}

  _nEvt ++ ;
}  



void KIT::check( LCEvent *  /*evt*/ ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void KIT::end(){}



