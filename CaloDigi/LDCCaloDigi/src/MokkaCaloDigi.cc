#include "MokkaCaloDigi.h"
#include <iostream>
#include <cmath>
#include <vector>
#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif
#include <IMPL/CalorimeterHitImpl.h>
#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <lcio.h>
#include <iterator>
#define SHIFT_M 0
#define SHIFT_S 3
#define SHIFT_I 6
#define SHIFT_J 15
#define SHIFT_K 24
#define SHIFT_2 30
#define SHIFT_1 31
#include "ced_cli.h"
//   Layer types: should be included
#define PHOTON_LAYER (1<<CED_LAYER_SHIFT)
#define NHADR_LAYER  (2<<CED_LAYER_SHIFT)
#define CHADR_LAYER  (3<<CED_LAYER_SHIFT)
#define TPC_LAYER    (4<<CED_LAYER_SHIFT)
#define ECAL_LAYER   (5<<CED_LAYER_SHIFT)
#define HCAL_LAYER   (6<<CED_LAYER_SHIFT)

#define MASK_M (unsigned int) 0x00000007
#define MASK_S (unsigned int) 0x00000038
#define MASK_I (unsigned int) 0x00007FC0
#define MASK_J (unsigned int) 0x00FF8000
#define MASK_K (unsigned int) 0x3F000000
#define MASK_2 (unsigned int) 0x40000000
#define MASK_1 (unsigned int) 0x80000000

using namespace lcio ;
using namespace marlin ;
using namespace std;
using namespace IMPL;

MokkaCaloDigi aMokkaCaloDigi ;


MokkaCaloDigi::MokkaCaloDigi() : Processor("MokkaCaloDigi") {
  
  // modify processor description
  _description = "Mokka digitizer..." ;
  

  // register steering parameters: name, description, class-variable, default value
 
  std::vector<std::string> ecalCollections;

  ecalCollections.push_back(std::string("ecal02_EcalBarrel"));
  ecalCollections.push_back(std::string("ecal02_EcalEndcap"));

  registerProcessorParameter( "ECALCollections" , 
			      "ECAL Collection Names" ,
			      _ecalCollections ,
			       ecalCollections);

  std::vector<std::string> hcalCollections;

  hcalCollections.push_back(std::string("hcalFeScintillator_HcalBarrelEnd"));
  hcalCollections.push_back(std::string("hcalFeScintillator_HcalBarrelReg"));
  hcalCollections.push_back(std::string("hcalFeScintillator_HcalEndCaps"));

  registerProcessorParameter("HCALCollections" , 
			     "HCAL Collection Names" , 
			     _hcalCollections , 
			     hcalCollections);

  registerProcessorParameter( "NewHCALCellSize" , 
			      "size of the new cell (integer) "  ,
			      cell_size ,
			      int(3) ) ;
  registerProcessorParameter( "NewHCALCollName" , 
			     "name for the new collection "  ,
			      _newCollNameHCAL ,
			      std::string("HCAL")) ;

  registerProcessorParameter( "NewECALCollName" , 
			      "name for the new collection "  ,
			      _newCollNameECAL ,
			      std::string("ECAL")) ;

  registerProcessorParameter("ECALThreshold" , 
			     "Threshold for ECAL Hits in GeV" ,
			     _thresholdEcal,
			     (float)1.0e-4);

  registerProcessorParameter("HCALThreshold" , 
			     "Threshold for HCAL Hits in GeV" ,
			     _thresholdHcal,
			     (float)4.0e-4);

  std::vector<int> ecalLayers;
  ecalLayers.push_back(30);
  ecalLayers.push_back(100);


  registerProcessorParameter("ECALLayers" , 
			     "Index of ECal Layers" ,
			     _ecalLayers,
			     ecalLayers);

  

  std::vector<int> hcalLayers;
  hcalLayers.push_back(100);

  registerProcessorParameter("HCALLayers" , 
			     "Index of HCal Layers" ,
			     _hcalLayers,
			     hcalLayers);


  std::vector<float> calibrEcal;
  calibrEcal.push_back(31.3);
  calibrEcal.push_back(83.0);


  registerProcessorParameter("CalibrECAL" , 
			     "Calibration coefficients for ECAL" ,
			     _calibrCoeffEcal,
			     calibrEcal);
  

  std::vector<float> calibrHcal;
  calibrHcal.push_back(27.3);

  registerProcessorParameter("CalibrHCAL" , 
			     "Calibration coefficients for HCAL" ,
			     _calibrCoeffHcal,
			     calibrHcal);


  registerProcessorParameter("IfDigitalEcal" ,
			     "Digital Ecal" , 
			     _digitalEcal , 
			     0);


  registerProcessorParameter("IfDigitalHcal" ,
			     "Digital Hcal" , 
			     _digitalHcal , 
			     0);


}


void MokkaCaloDigi::init() { 

  // usually a good idea to
   printParameters() ;
   //ced_client_init("localhost",7286);
   //  ced_register_elements();

  _nRun = -1 ;
  _nEvt = 0 ;
  //  cell_size=4;
}

void MokkaCaloDigi::processRunHeader( LCRunHeader* run) { 

  _nRun++  ;
} 

void MokkaCaloDigi::processEvent( LCEvent * evt ) { 
  //ced_new_event();
    vector <CalorimeterHitImpl*> hitall;
    vector <CalorimeterHitImpl*> hitnew;
    SimCalorimeterHit* tmp_hit;
    int Asell;
    int zamena;
    int ID_k;
    int ID_i;
    int ID_j;
    int ID_M;
    int ID_S;   
    int ID_iprime;  
    int ID_kprime;
    int ID_jprime;
    int ID_Mprime;
    int ID_Sprime;
    unsigned int N_coll;  
    float Eold;
    float Eadd;
    float Esum;
    float point[3];
    float point_new[3];
    float teta_r;
    int Nmax_z;int z_last_cell_size;
    Asell=cell_size*cell_size;
    int xl[40]={ 159, 161, 163, 165, 167, 169, 171, 173, 175,177, 179,
                     181,183,185,187,189,191,193,195,197,200,202,
                     204,206,208,210,212,214,216,218,220,222,211,
                     200,188,176,164,152,140,129};
    int  zoff[40]={0,0,0,0,0,0,3,6,10,10,10,10,10,10,10,10,10,10,10,10,
                     10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10};

    Nmax_z=106/cell_size;
    z_last_cell_size=106%cell_size;

    int Nxmax[40];
    int xlcs[40];
    for(int i=0;i<40;i++)
      {
	Nxmax[i]=xl[i]/cell_size;
	xlcs[i]=xl[i]%cell_size;
      }           
    int  Nzmax[40];
    int zlcs[40];
    for(int i=0;i<40;i++)
      {
	Nzmax[i]=(106+zoff[i])/cell_size;
	zlcs[i]=(106+zoff[i])%cell_size;
      }  
   int cellid ;
   int z_id_swap[40];   
   for (int i=0;i<40;i++)
       { z_id_swap[i]=105+zoff[i];  }            
   /********************* initial calculation done *********/
    int loop_off;
    hitnew.erase(hitall.begin(),hitall.end());
    hitall.erase(hitall.begin(),hitall.end());

     for ( N_coll=0;N_coll<_hcalCollections.size();++N_coll)
      {   
	loop_off=hitall.size();
        //  hitall.erase(hitall.begin(),hitall.end());                    
	try{

             LCCollection* col = evt->getCollection( (_hcalCollections[N_coll]).c_str() ) ;
             int N = col->getNumberOfElements();

        for (int i=0;i<N;i++)
        {
	tmp_hit=(dynamic_cast<SimCalorimeterHit*>( col->getElementAt(i) ));   

	CalorimeterHitImpl * tmp_fill = new  CalorimeterHitImpl();

	  
	tmp_fill->setEnergy(tmp_hit->getEnergy());
	tmp_fill->setPosition(tmp_hit->getPosition());
	tmp_fill->setCellID0(tmp_hit->getCellID0());
	

	point[0] =  tmp_hit->getPosition()[0];
	point[1] =  tmp_hit->getPosition()[1];
	point[2] = tmp_hit->getPosition()[2];
	      //  unsigned int kolor=0xf2290e;
	      // ced_hit(point[0],point[1],point[2],0,3,kolor);

	hitall.push_back(tmp_fill);           
	}

       for(int i=0; i<N;i++)     
             {
	       cellid = (hitall[i+loop_off])->getCellID0();
	       ID_M=(cellid & MASK_M) >> SHIFT_M; // reed module number on it depends further calculation
	       ID_j=(cellid & MASK_J) >> SHIFT_J;
	       ID_S=(cellid & MASK_S) >> SHIFT_S;	
	       ID_k=(cellid & MASK_K) >> SHIFT_K;       
	       ID_i=(cellid & MASK_I) >> SHIFT_I;
	       //	         cout << ID_M << " "<< ID_S << " "<< ID_i << " "<< ID_j<< " "<< ID_k << endl;
              
	       switch ( ID_M )
               {
	     
               
	      case 2: case 3: case 4: case 5:		
		 zamena=ID_j;
                 ID_j=ID_k;
                 ID_k=zamena;
                 ID_k=ID_k/cell_size;
                 ID_i=ID_i/cell_size; 
               break;
               case 1:
		
		 zamena=ID_j;
                 ID_j=ID_k;
                 ID_k=zamena;
                  ID_i=ID_i/cell_size;                  
                  ID_k=z_id_swap[ID_j]-ID_k;
	          ID_k=ID_k/cell_size; 
		 
               break;
       
               case 0: case 6:
                 zamena=ID_i;
                 ID_i=ID_j;
	         ID_j=zamena;
		 switch (ID_S) 
		   {
                   case 0:case 1:case 2:
                     if ( ID_j >= 30 )
		       { ID_j=ID_j-30;
                         ID_i=ID_i+30;
		       }
                     else 
                       {
			 ID_S=ID_S+1;  
                         zamena=ID_i-30;
                         ID_i=29-ID_j;
                         ID_j=zamena;
		       }
                   break;

                    case 3:
                    if ( ID_j >= 30 )
		       { ID_j=ID_j-30;
                         ID_i=ID_i+30;
		       }
                     else 
                       {
			 ID_S=0;  
                         zamena=ID_i-30;
                         ID_i=29-ID_j;
                         ID_j=zamena;
		       }
		   break;                              
	           } //swith ID_S

                 ID_j=ID_j/cell_size;
                 ID_i=ID_i/cell_size; 
                
	       break;  
	     
               }//switch

               cellid=0;
               cellid=((ID_M<<SHIFT_M)&MASK_M)|((ID_S<<SHIFT_S)&MASK_S)|((ID_j<<SHIFT_J)&MASK_J)|((ID_k<<SHIFT_K)&MASK_K)|((ID_i<<SHIFT_I)&MASK_I);  // calculate new ID 
	       (hitall[i+loop_off])->setCellID0(cellid); //assign new ID to hits
              (hitall[i+loop_off])->setCellID1(1);
             } // loop for calculating new index 
       //    cout <<"before"<< endl;

       for (int i=0; i<N;i++) // loop now and merge cells with same ID
       {	       
        if ( (hitall[i+loop_off])->getCellID1())  
	      {
		cellid=(hitall[i+loop_off])->getCellID0();
                ID_k=(cellid & MASK_K) >> SHIFT_K;
		ID_i=(cellid & MASK_I) >> SHIFT_I;
		ID_j=(cellid & MASK_J) >> SHIFT_J;  
		ID_M=(cellid & MASK_M) >> SHIFT_M;      
		ID_S=(cellid & MASK_S) >> SHIFT_S;
	
	 for (int kk=i+1;kk<N;kk++)
	   {

              cellid=(hitall[kk+loop_off])->getCellID0();
	       	ID_kprime=(cellid & MASK_K) >> SHIFT_K;
		ID_iprime=(cellid & MASK_I) >> SHIFT_I;
		ID_jprime=(cellid & MASK_J) >> SHIFT_J;  
		ID_Mprime=(cellid & MASK_M) >> SHIFT_M;      
		ID_Sprime=(cellid & MASK_S) >> SHIFT_S;


          if ( ( ID_kprime == ID_k) && (ID_iprime == ID_i)
	     &&( ID_jprime == ID_j) && (ID_Sprime == ID_S ) && (ID_Mprime == ID_M))
	    {
		  Eold=float((hitall[i+loop_off])->getEnergy());
		  Eadd=float((hitall[kk+loop_off])->getEnergy());
                  Esum=Eold+Eadd;  
                                
                  (hitall[i+loop_off])->setEnergy(Esum) ;
                  (hitall[kk+loop_off])->setCellID1(0); // do not count this any more !
             }
	   } //loop over kk   
	 //now calculate new coordinates of cell
             
	 switch(ID_M)
           {
           case 2:case 3: case 4:
             point[1]=1931.25+24.5*ID_j;	   
            if ((ID_i<Nxmax[ID_j])&&(ID_k<Nmax_z))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+cell_size*0.5;
                point[2]=(ID_M-3)*108.20-53.0+ID_k*cell_size+cell_size*0.5;
               (hitall[i+loop_off])->setCellID1(Asell); 
              }
	    if ((ID_i==Nxmax[ID_j])&&(ID_k<Nmax_z))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+xlcs[ID_j]*0.5;
                point[2]=(ID_M-3)*108.20-53.0+ID_k*cell_size+cell_size*0.5;
                (hitall[i+loop_off])->setCellID1(xlcs[ID_j]*cell_size); 
              }
	    if ((ID_i<Nxmax[ID_j])&&(ID_k==Nmax_z))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+cell_size*0.5;
                point[2]=(ID_M-3)*108.20-53.0+ID_k*cell_size+z_last_cell_size*0.5;
                (hitall[i+loop_off])->setCellID1(cell_size*z_last_cell_size); 
              
              }
	    if ((ID_i==Nxmax[ID_j])&&(ID_k==Nmax_z))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+xlcs[ID_j]*0.5;
                point[2]=(ID_M-3)*108.20-53.0+ID_k*cell_size+z_last_cell_size*0.5;
                 (hitall[i+loop_off])->setCellID1(xlcs[ID_j]*z_last_cell_size); 
              } 
            teta_r=-ID_S*0.7853981633974483;
	    point_new[2]=point[2]*10;
	    point[0]=point[0]*10;
	    point_new[0]=point[0]*cos(teta_r)+point[1]*sin(teta_r);
	    point_new[1]=-point[0]*sin(teta_r)+point[1]*cos(teta_r);             
	    (hitall[i+loop_off])->setPosition(point_new);
	   break;
          
           case 1:
	     //  cout << ID_k << " ";
	    point[1]=1931.25+24.5*ID_j;
            if ((ID_i<Nxmax[ID_j])&&(ID_k<Nzmax[ID_j]))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+cell_size*0.5;
                point[2]=(ID_M-3)*108.20+53.0-ID_k*cell_size-cell_size*0.5;
		(hitall[i+loop_off])->setCellID1(Asell); 
              }
	    if ((ID_i==Nxmax[ID_j])&&(ID_k<Nzmax[ID_j]))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+xlcs[ID_j]*0.5;
                point[2]=(ID_M-3)*108.20+53.0-ID_k*cell_size-cell_size*0.5;
		(hitall[i+loop_off])->setCellID1(xlcs[ID_j]*cell_size);
              }
	    if ((ID_i<Nxmax[ID_j])&&(ID_k==Nzmax[ID_j]))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+cell_size*0.5;
                point[2]=(ID_M-3)*108.20+53.0-ID_k*cell_size-zlcs[ID_j]*0.5;
              (hitall[i+loop_off])->setCellID1(cell_size*zlcs[ID_j]); 
              
              }
	    if ((ID_i==Nxmax[ID_j])&&(ID_k==Nzmax[ID_j]))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+xlcs[ID_j]*0.5;
                point[2]=(ID_M-3)*108.2+53.0-ID_k*cell_size-zlcs[ID_j]*0.5;
		(hitall[i+loop_off])->setCellID1(xlcs[ID_j]*zlcs[ID_j]);
              } 
            teta_r=-ID_S*0.7853981633974483;
	    point_new[2]=point[2]*10;
	    point[0]=point[0]*10;
	    point_new[0]=point[0]*cos(teta_r)+point[1]*sin(teta_r);
	    point_new[1]=-point[0]*sin(teta_r)+point[1]*cos(teta_r);             
	    (hitall[i+loop_off])->setPosition(point_new);
	   break;
           
           case 5:

	    point[1]=1931.25+24.5*ID_j;
            if ((ID_i<Nxmax[ID_j])&&(ID_k<Nzmax[ID_j]))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+cell_size*0.5;
                point[2]=(ID_M-3)*108.20-53.0+ID_k*cell_size+cell_size*0.5-zoff[ID_j];
		(hitall[i+loop_off])->setCellID1(Asell); 
              }
	    if ((ID_i==Nxmax[ID_j])&&(ID_k<Nzmax[ID_j]))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+xlcs[ID_j]*0.5;
                point[2]=(ID_M-3)*108.20-53.0+ID_k*cell_size+cell_size*0.5-zoff[ID_j];
		(hitall[i+loop_off])->setCellID1(xlcs[ID_j]*cell_size);
              }
	    if ((ID_i<Nxmax[ID_j])&&(ID_k==Nzmax[ID_j]))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+cell_size*0.5;
                point[2]=(ID_M-3)*108.20-53.0+ID_k*cell_size+zlcs[ID_j]*0.5-zoff[ID_j];
		(hitall[i+loop_off])->setCellID1(cell_size*zlcs[ID_j]);
              }
	    if ((ID_i==Nxmax[ID_j])&&(ID_k==Nzmax[ID_j]))
              {
		point[0]=(-xl[ID_j]/2.0)+ID_i*cell_size+xlcs[ID_j]*0.5;
                point[2]=(ID_M-3)*108.20-53.0+ID_k*cell_size+zlcs[ID_j]*0.5-zoff[ID_j];
		(hitall[i+loop_off])->setCellID1(xlcs[ID_j]*zlcs[ID_j]);
              } 
	    teta_r=-ID_S*0.7853981633974483;
	    point_new[2]=point[2]*10;
	    point[0]=point[0]*10;
	    point_new[0]=point[0]*cos(teta_r)+point[1]*sin(teta_r);
	    point_new[1]=-point[0]*sin(teta_r)+point[1]*cos(teta_r);             
	    (hitall[i+loop_off])->setPosition(point_new);
           break;
           
           case 6:
	     point[0]=(ID_i*cell_size+cell_size*0.5)*10;
	     point[1]=(ID_j*cell_size+cell_size*0.5)*10;
	    
	        if (ID_M == 0) 
		  {
		    point_new[2]=-3039.25-24.5*ID_k;
                  }
                else 
                  {
                    point_new[2]=3039.25+24.5*ID_k; 
		  }

              switch(ID_S ) {
                  case 0:
                   point_new[0]=point[0]-300;
                   point_new[1]=point[1]+300;
		  break;
		  case 1:
		   point_new[0]=point[1]+300;
                   point_new[1]=-point[0]+300;
		  break;
		  case 2:
                   point_new[0]=-point[0]+300;
                   point_new[1]=-point[1]-300;
		  break;
		  case 3:
		   point_new[0]=-point[1]-300;
                   point_new[1]=point[0]-300;
                  break;
	          }   
              // (hitall[i])->setCellID1(); 	
	      (hitall[i+loop_off])->setCellID1(Asell); 
	      (hitall[i+loop_off])->setPosition(point_new);         
	   break;
	    case 0:
	     point[0]=(ID_i*cell_size+cell_size*0.5)*10;
	     point[1]=(ID_j*cell_size+cell_size*0.5)*10;
	    
	        if (ID_M == 0) 
		  {
		    point_new[2]=-3039.25-24.5*ID_k;
                  }
                else 
                  {
                    point_new[2]=3039.25+24.5*ID_k; 
		  }

              switch(ID_S ) {
                  case 0:
                   point_new[0]=-point[0]+300;
                   point_new[1]=point[1]+300;
		  break;
		  case 1:
		   point_new[0]=-point[1]-300;
                   point_new[1]=-point[0]+300;
		  break;
		  case 2:
                   point_new[0]=point[0]-300;
                   point_new[1]=-point[1]-300;
		  break;
		  case 3:
		   point_new[0]=point[1]+300;
                   point_new[1]=point[0]-300;
                  break;
	          }    
             	(hitall[i+loop_off])->setCellID1(Asell); 
	      (hitall[i+loop_off])->setPosition(point_new);         
	   break;        
           }//switch
             
             // now we have a new coordinate and new energy and we store the new hit 
	 //unsigned int kolor=0x25cb17;

	  //    ced_hit(point_new[0],point_new[1],point_new[2],0,3,kolor);
             hitnew.push_back(hitall[i+loop_off]);

	   }//if
          }//loop over all hits
       //  ced_send_event();
    

	     }//try 
         catch(DataNotAvailableException &e){}

      } // loop over all the collections 
     //ced_draw_event();
   //  cout <<endl;
           LCCollectionVec* CalorimeterHitgang = new LCCollectionVec( LCIO::CALORIMETERHIT);
           LCFlagImpl flag;
	   flag.setBit(LCIO::CHBIT_LONG);
	   CalorimeterHitgang->setFlag(flag.getFlag());



      for(unsigned int i=0;i<hitnew.size();i++)
         {
             cellid = (hitnew[i])->getCellID0();
	     ID_M=(cellid & MASK_M) >> SHIFT_M;
	       ID_j=(cellid & MASK_J) >> SHIFT_J;
	       ID_S=(cellid & MASK_S) >> SHIFT_S;	
	       ID_k=(cellid & MASK_K) >> SHIFT_K;       
	       ID_i=(cellid & MASK_I) >> SHIFT_I;
	       
               cellid=0;
	       switch ( ID_M )
               {
	     
	
	       case 1:case 2: case 3: case 4: case 5:		
		 zamena=ID_k;
                 ID_k=ID_j;
                 ID_j=zamena;

               cellid=((ID_M<<SHIFT_M)&MASK_M)|((ID_S<<SHIFT_S)&MASK_S)|((ID_j<<SHIFT_J)&MASK_J)|((ID_k<<SHIFT_K)&MASK_K)|((ID_i<<SHIFT_I)&MASK_I);  
	       (hitnew[i])->setCellID0(cellid);
	       // cout<< ID_k << " "; 	
                  break;
               default:
                cellid=((ID_M<<SHIFT_M)&MASK_M)|((ID_S<<SHIFT_S)&MASK_S)|((ID_j<<SHIFT_J)&MASK_J)|((ID_k<<SHIFT_K)&MASK_K)|((ID_i<<SHIFT_I)&MASK_I);  
	       (hitnew[i])->setCellID0(cellid);
	       // cout<< ID_k << " "; 
		 break;                  
	       }
	       float energy = hitnew[i]->getEnergy();

	       if (energy > _thresholdHcal) {

		 int cellid = hitnew[i]->getCellID0();
		 float calibr_coeff(1.);
		 int layer = cellid >> 24 ;
		 for (unsigned int k(0); k < _hcalLayers.size(); ++k) {
		   int min,max;
		   if (k == 0) 
		     min = 0;
		   else 
		     min = _hcalLayers[k-1];
		   max = _hcalLayers[k];
		   if (layer >= min && layer < max) {
		     calibr_coeff = _calibrCoeffHcal[k];
		     break;
		   }
		 } 
		 if (_digitalHcal) {
		   hitnew[i]->setEnergy(calibr_coeff); 
		 }
		 else {
		   hitnew[i]->setEnergy(calibr_coeff*energy);
		 }


		 CalorimeterHitgang->addElement(hitnew[i]);
	       }
         }
       evt->addCollection( CalorimeterHitgang , _newCollNameHCAL.c_str() )  ;

      


  LCCollectionVec *ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
 
  ecalcol->setFlag(flag.getFlag());
 


// 
// * Reading Collections of ECAL Simulated Hits * 
// 

  for (unsigned int i(0); i < _ecalCollections.size(); ++i) {
      try{
	  LCCollection * col = evt->getCollection( _ecalCollections[i].c_str() ) ;
	  int numElements = col->getNumberOfElements();

	  for (int j(0); j < numElements; ++j) {
	      SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	     
                    
	      float energy = hit->getEnergy();
	      if (energy > _thresholdEcal) {
		CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
		int Cellid = hit->getCellID0();
		float calibr_coeff(1.);
		int layer = Cellid >> 24;
		for (unsigned int k(0); k < _ecalLayers.size(); ++k) {
		  int min,max;
		  if (k == 0) 
		    min = 0;		      
		  else 
		    min = _ecalLayers[k-1];		      
		  max = _ecalLayers[k];
		  if (layer >= min && layer < max) {
		    calibr_coeff = _calibrCoeffEcal[k];
		    break;
		  }
		} 
		calhit->setCellID0(Cellid);
		if (_digitalEcal) {
		  calhit->setEnergy(calibr_coeff); 
		}
		else {
		  calhit->setEnergy(calibr_coeff*energy);
		}

		calhit->setPosition(hit->getPosition());
		ecalcol->addElement(calhit);
	      }
	  }
      }
      catch(DataNotAvailableException &e){ 
      }
  }

  evt->addCollection(ecalcol, _newCollNameECAL.c_str());

  _nEvt ++ ;

} // end of event processor 


void MokkaCaloDigi::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void MokkaCaloDigi::end(){ 
  
//   std::cout << "MyProcessor::end()  " << name() 
// 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
// 	    << std::endl ;

}

