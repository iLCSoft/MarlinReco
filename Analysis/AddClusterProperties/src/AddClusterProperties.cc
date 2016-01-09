#include "AddClusterProperties.h"
#include "WeightedPoints3D.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Cluster.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include "marlin/Exceptions.h"
#include <vector>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA


using namespace lcio ;
using namespace marlin ;


AddClusterProperties aAddClusterProperties ;


AddClusterProperties::AddClusterProperties() : Processor("AddClusterProperties") {

    // modify processor description
    _description = "AddClusterProperties does whatever it does ..." ;


    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "PFOCollectionName" ,
                           "Name of the input PFO collection"  ,
                           _PFOName ,
                           std::string("PandoraPFOs") ) ;
  
    registerInputCollection( LCIO::CLUSTER,
                           "ClusterCollection" , 
                           "Name of the Clusters input collection"  ,
                           _clusterCollectionName ,
                           std::string("PandoraClusters") ) ;


}







void AddClusterProperties::init() { 

    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    // usually a good idea to
    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}


void AddClusterProperties::processRunHeader( LCRunHeader* run) { 

    _nRun++ ;
} 



void AddClusterProperties::processEvent( LCEvent * evt ) { 


    int npoints_index=0;
    int sum_wgt_index=0 ;
    int sum_wgt2_index=0 ;
    int sum_wgt4_index=0 ;
    int ecal_index=0 ;
    int hcal_index=0 ;
    int yoke_index=0 ;
    int lcal_index=0 ;
    int lhcal_index=0 ;
    int bcal_index=0 ;
    StringVec shape_keys ;
    LCCollection* clucol = NULL;
    try{
        clucol = evt->getCollection( _clusterCollectionName  );
        StringVec shapeParams ;
        shapeParams = clucol->getParameters().getStringVals ("ClusterShapeParameters", shapeParams);
        shapeParams.push_back("npoints") ;
        shapeParams.push_back("sum_wgt") ;
        shapeParams.push_back("sum_wgt^2") ;
        shapeParams.push_back("sum_wgt^4") ;
        clucol->parameters().setValues( "ClusterShapeParameters" , shapeParams ) ;
        shape_keys = clucol->getParameters().getStringVals( std::string("ClusterShapeParameters"),shape_keys);
        for ( unsigned kkk=0 ; kkk < shape_keys.size() ; kkk++ ) {
	  if ( shape_keys[kkk] == "npoints" )   { npoints_index  = kkk ; }
	  if ( shape_keys[kkk] == "sum_wgt" ) { sum_wgt_index = kkk ; }
	  if ( shape_keys[kkk] == "sum_wgt^2" ) { sum_wgt2_index = kkk ; }
          if ( shape_keys[kkk] == "sum_wgt^4" ) { sum_wgt4_index = kkk ; }
        }
        StringVec subDetectorNames;
        subDetectorNames = clucol->getParameters().getStringVals ("ClusterSubdetectorNames", subDetectorNames);
        for ( unsigned kkk=0 ; kkk <  subDetectorNames.size() ; kkk++ ) {
	  if ( subDetectorNames[kkk] == "ecal"  )   { ecal_index  = kkk ;  }
	  if ( subDetectorNames[kkk] == "hcal"  )   { hcal_index  = kkk ;  }
	  if ( subDetectorNames[kkk] == "yoke"  )   { yoke_index  = kkk ;  }
	  if ( subDetectorNames[kkk] == "lcal"  )   { lcal_index  = kkk ;  }
	  if ( subDetectorNames[kkk] == "lhcal" )   { lhcal_index = kkk ;  }
	  if ( subDetectorNames[kkk] == "bcal"  )   { bcal_index  = kkk ;  }
        }
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) <<  _clusterCollectionName   << " collection not available" << std::endl;
        clucol = NULL;
    }

    if( clucol != NULL ){
      int nclu = clucol->getNumberOfElements()  ;
      for(int i=0; i< nclu ; i++){
        streamlog_out(DEBUG3) << "i cluster " << i << std::endl;
        ClusterImpl* clu = dynamic_cast<ClusterImpl*>( clucol->getElementAt( i ) ) ;
        CalorimeterHitVec hits=clu->getCalorimeterHits();
        int np=clu->getCalorimeterHits().size();
        if ( np == 0 ) {
          streamlog_out(ERROR) << " Cluster with no hits. Input probably isnt a RECO-file, so just stop here ! " << std::endl;
          throw StopProcessingException( dynamic_cast<marlin::Processor*>(this) );  
        } 

        // get the hits: 

	std::vector<double> ehit,xhit,yhit,zhit;
        for (int ihit=0; ihit< np ; ihit++){
          CalorimeterHit* calo_hit  = dynamic_cast<CalorimeterHit*>(hits[ihit]);
          ehit.push_back(calo_hit->getEnergy());
          xhit.push_back(calo_hit->getPosition()[0]);
          yhit.push_back(calo_hit->getPosition()[1]);
          zhit.push_back(calo_hit->getPosition()[2]);
        }

        // Analayse the cluster:

	WeightedPoints3D wgtp = WeightedPoints3D( int(hits.size()) , &ehit[0] , &xhit[0], &yhit[0], &zhit[0] ); 

        double* cog=wgtp.getCentreOfGravity();

        double *covv=wgtp.getCentreOfGravityErrors();
        double cov[3][3]; for( int iii=0 ; iii<3 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { cov[iii][jjj]=covv[jjj+iii*3] ; } } ; 

        double *eval=wgtp.getEigenVal(); 
        double *eval_err=wgtp.getEigenValErrors(); 

        double *evpv=wgtp.getEigenVecPolar(); 
	double evp[2][3]; for( int iii=0 ; iii<2 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { evp[iii][jjj]=evpv[jjj+iii*3] ; } } ;

        double *evpc=wgtp.getEigenVecCartesian(); 
	double evc[2][3]; for( int iii=0 ; iii<2 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { evc[iii][jjj]=evpc[jjj+iii*3] ; } } ;

        double* evpev=wgtp.getEigenVecPolarErrors();              
	double evpe[2][2][3]; for( int iii=0 ; iii<2 ; iii++ ) { for ( int jjj=0 ; jjj<2 ; jjj++ ) {for ( int kkk=0 ; kkk<3 ; kkk++ ) { evpe[iii][jjj][kkk]=evpev[kkk+jjj*3+iii*3*2] ; } } };

        double sum_wgtsqr = wgtp.getTotalSquaredWeight();
        double sum_wgt4 = wgtp.getTotalQuarticWeight();
        double sum_wgt = wgtp.getTotalWeight();
        streamlog_out(DEBUG3) << "i Totals : "  << sum_wgt << " " << clu->getEnergy() << std::endl;  
	//        if ( wgtp ) delete wgtp ;

        // Pack/translate/cast/... :

        FloatVec PositionError ;
        PositionError.push_back(cov[0][0]*sum_wgtsqr/(sum_wgt*sum_wgt));
	PositionError.push_back(cov[0][1]*sum_wgtsqr/(sum_wgt*sum_wgt));
	PositionError.push_back(cov[1][1]*sum_wgtsqr/(sum_wgt*sum_wgt));
	PositionError.push_back(cov[0][2]*sum_wgtsqr/(sum_wgt*sum_wgt));
	PositionError.push_back(cov[1][2]*sum_wgtsqr/(sum_wgt*sum_wgt));
	PositionError.push_back(cov[2][2]*sum_wgtsqr/(sum_wgt*sum_wgt));
        streamlog_out(DEBUG3) << "i PositionError : "  << std::endl;  streamlog_out(DEBUG3) << "i " ;
        for ( int iii=0 ; iii < 6 ; iii++ ) { streamlog_out(DEBUG3) << PositionError[iii] << " " ;}  ;streamlog_out(DEBUG3) << std::endl;

        float Position[3];
        Position[0] = cog[0];
        Position[1] = cog[1];
        Position[2] = cog[2];
        streamlog_out(DEBUG3) << "i Position : "  << std::endl;  streamlog_out(DEBUG3) << "i " ;
        for ( int iii=0 ; iii < 3 ; iii++ ) { streamlog_out(DEBUG3) << Position[iii] << " " ;}  ;streamlog_out(DEBUG3) << std::endl;

        float theta = evp[0][2];
        float phi = evp[1][2];
        streamlog_out(DEBUG3) << "i theta/phi : "  << std::endl;  streamlog_out(DEBUG3) << "i " ;
        streamlog_out(DEBUG3) << theta << " "  << phi <<  std::endl;

        FloatVec DirectionError ;
	DirectionError.push_back(evpe[0][0][2]);
	DirectionError.push_back(evpe[1][0][2]);
	DirectionError.push_back(evpe[1][1][2]);
        streamlog_out(DEBUG3) << "i DirectionError : "  << std::endl;  streamlog_out(DEBUG3) << "i " ;
        for ( int iii=0 ; iii < 3 ; iii++ ) { streamlog_out(DEBUG3) << DirectionError[iii] << " " ;}  ;streamlog_out(DEBUG3) << std::endl;
 
        FloatVec shape = clu->getShape() ; 
	shape.resize(shape_keys.size());
	shape[npoints_index]=1.0*hits.size() ; 
	shape[sum_wgt_index]=sum_wgt; 
	shape[sum_wgt2_index]=sum_wgtsqr; 
        shape[sum_wgt4_index]=sum_wgt4; 
 

        float Eerror=clu->getEnergyError();
        if ( Eerror == 0.0 ) {
            // not set, so as per HLRWS:
          float E=clu->getEnergy();
          const FloatVec pec  = clu->getSubdetectorEnergies();
          float Eem = pec[ecal_index]+pec[lcal_index]+pec[bcal_index];
          float Ehad = pec[hcal_index]+pec[yoke_index]+pec[lhcal_index];
          float Eerror=0.0;
          if ( Eem/E < 0.95 ) {
            Eerror=sqrt((0.6*sqrt(E))*(0.6*sqrt(E)) + (0.03*E)*(0.03*E));
          } else {
            Eerror=sqrt((0.17*sqrt(E))*(0.17*sqrt(E)) + (0.01*E)*(0.01*E));
          }
          clu->setEnergyError(Eerror);
	  streamlog_out(DEBUG3) << " energies : " << E << " " << Eem << " " << Ehad << " " << Eerror << " " << Eerror/sqrt(E) << std::endl;
        }

        for ( int kkk = npoints_index ; kkk <=sum_wgt4_index ; kkk++ ) {if (shape_keys[kkk] !=  shape_keys[kkk] ) { streamlog_out(WARNING) << " shape_keys " << kkk << " is NaN " << std::endl ; }}
        for ( int kkk = 0 ; kkk < 2 ; kkk++ ) {if ( Position[kkk] !=  Position[kkk] ) { streamlog_out(WARNING)<< " Position " << kkk << " is NaN " << std::endl ; }}
        for ( int kkk = 0 ; kkk < 6 ; kkk++ ) {if ( PositionError[kkk] !=  PositionError[kkk] ) { streamlog_out(WARNING) << " PositionError " << kkk << " is NaN " << std::endl ; }}
        for ( int kkk = 0 ; kkk < 3 ; kkk++ ) {if ( DirectionError[kkk] !=  DirectionError[kkk] ) { streamlog_out(WARNING) << " DirectionError " << kkk << " is NaN " << std::endl ; }}
        if ( theta !=    theta )    {  streamlog_out(WARNING)  << " theta is Nan " << std::endl ; }
        if ( phi   !=    phi   )    {  streamlog_out(WARNING)  << " phi is Nan " << std::endl ; }

        // add to cluster:

        clu->setShape( shape ) ;
        clu->setPosition(Position);
        clu->setPositionError(&PositionError[0]);
        clu->setITheta(theta);
        clu->setIPhi(phi);
        clu->setDirectionError(&DirectionError[0]);

        if ( streamlog::out.write<streamlog::DEBUG2>() ){
          debuging(clucol,clu,cog,&cov[0][0],eval,eval_err,&evp[0][0],&evpe[0][0][0],&evc[0][0],np,sum_wgt,sum_wgtsqr,sum_wgt4);
        }

      }
    }

    // Now loop PFOs to add the error on P for neutral ones, from the cluster
    // NB: not enough to just use the cluster COG + errors + error on E: Also must
    // check if the neutral actually is a V0 or gamma-conversion, in which case the
    // p should be from the fitted values (or at least the sum on the tracks)

    LCCollection* pfocol = NULL;
    try{
      pfocol = evt->getCollection( _PFOName );
    }
    catch( lcio::DataNotAvailableException e )
    {
        streamlog_out(WARNING) << _PFOName    << " collection not available" << std::endl;
        pfocol = NULL;
    }
    if( pfocol != NULL ){
      int npfo = pfocol->getNumberOfElements()  ;
      for(int i=0; i< npfo ; i++){
             //  get the PFO...
        ReconstructedParticleImpl* part = dynamic_cast<ReconstructedParticleImpl*>( pfocol->getElementAt( i ) ) ;
             //   and its clusters ...
        ClusterVec clusters=part->getClusters();
        int nclu =  part->getClusters().size();
        TrackVec tracks=part->getTracks();
        int ntrk =  part->getTracks().size();
        if ( part->getCharge() == 0 && nclu == 1 && ntrk == 0 ) {
          // for now, only look at the clear-cut cases - see above
          streamlog_out(DEBUG3) << "i normal, neutral pfo " << i << std::endl;
          Cluster* clu= dynamic_cast<Cluster*>(clusters[0]) ;
          const float* cogv = clu->getPosition();
          float dist = 0.0;
          for ( int iii=0 ; iii < 3 ; iii++ ) {
            dist+=cogv[iii]*cogv[iii];
          }
          dist = sqrt(dist);
          float Eclu = clu->getEnergy();
          float E = part->getEnergy();
	  streamlog_out(DEBUG3) << " pfo/cluster E " << E << " / " << Eclu << std::endl;
          float Eerror = clu->getEnergyError();
          const double* mompfo=part->getMomentum();
          double mom[4] ;
          for ( int kkk=0 ; kkk<3 ; kkk++) {mom[kkk] = E*cogv[kkk]/dist;}
          mom[3] = E ;
          streamlog_out(DEBUG3)  << "i clu-mom: " ; for (int kkk=0 ; kkk<4 ; kkk++ ) { streamlog_out(DEBUG3)  << mom[kkk] << " " ;} streamlog_out(DEBUG3)  << std::endl;
          streamlog_out(DEBUG3)  << "i pfo-mom: " ; for (int kkk=0 ; kkk<4 ; kkk++ ) { streamlog_out(DEBUG3)  << mompfo[kkk] << " " ;} streamlog_out(DEBUG3)  << std::endl;
          part->setMomentum(mom);
          FloatVec PositionError = clu->getPositionError();
          double covmat[4][4];
          int nnn = 0;
          for ( int kkk=0 ; kkk < 3 ; kkk++ ) {
            for ( int lll=0 ; lll<=kkk ; lll++) {
              covmat[lll][kkk] =  PositionError[nnn];
  	      covmat[kkk][lll] =  PositionError[nnn];
              nnn++ ;
            }
          }
          for ( int kkk=0 ; kkk < 4 ; kkk++ ) {
            covmat[kkk][3] = 0. ;  
            covmat[3][kkk] = 0. ;  
          }  
          covmat[3][3] = Eerror*Eerror;
	  streamlog_out(DEBUG3)  << "i covmat : " << std::endl;
	  for ( int kkk=0 ; kkk <4 ; kkk++ ) {
  	    for ( int iii=0 ; iii <4 ; iii++ ) { streamlog_out(DEBUG3) << covmat[kkk][iii] << " " ;}
            streamlog_out(DEBUG3) << std::endl; 
          }  

          double dp_drE[4][4] ;
          double prefact= E/(dist*dist*dist);
          for ( int kkk=0 ; kkk < 4 ; kkk++ ) {
            if ( kkk < 3 ) {
              for ( int lll=0 ; lll< 3 ; lll++) {
                if ( lll == kkk ) {
                  dp_drE[lll][kkk] = prefact*(dist*dist+cogv[kkk]*cogv[kkk]);
                } else {
                  dp_drE[lll][kkk] = prefact*cogv[kkk]*cogv[lll];
                }
              }
            } else {
              for ( int lll=0 ; lll< 3 ; lll++) {
                dp_drE[lll][3] = 0 ;
              }
              for ( int lll=0 ; lll< 3 ; lll++) {
                dp_drE[3][lll] = prefact*cogv[lll]*dist*dist/E ;
              }
	      dp_drE[3][3] = prefact*dist*dist*dist/E;
            }
          }
	  streamlog_out(DEBUG3)  << "i dp_drE : " << std::endl;
	  for ( int kkk=0 ; kkk <4 ; kkk++ ) {
  	    for ( int iii=0 ; iii <4 ; iii++ ) { streamlog_out(DEBUG3) << dp_drE[kkk][iii] << " " ;}
            streamlog_out(DEBUG3) << std::endl; 
          }
           // propagate
          double p_cov[4][4];
          double tmp_mat [4][4];
          for ( int iii =0 ; iii < 4 ; iii++ ) {
            for ( int jjj=0 ; jjj < 4 ; jjj++ ) {
              tmp_mat [iii][jjj] = 0. ;
              for (int kkk=0 ; kkk < 4 ; kkk++ ) {
                tmp_mat [iii][jjj] += covmat[iii][kkk]* dp_drE[kkk][jjj];
              }
            }
          }
          for ( int iii =0 ; iii < 4 ; iii++ ) {
            for ( int jjj=0 ; jjj < 4 ; jjj++ ) {
              p_cov[iii][jjj] = 0. ;
              for (int kkk=0 ; kkk < 4 ; kkk++ ) {
	        p_cov[iii][jjj] += dp_drE[kkk][iii]*tmp_mat[kkk][jjj] ;
              }
            }
          }
	  streamlog_out(DEBUG3)  << "i pcov : " << std::endl;
	  for ( int kkk=0 ; kkk <4 ; kkk++ ) {
  	    for ( int iii=0 ; iii <4 ; iii++ ) { streamlog_out(DEBUG3) << "i " << p_cov[kkk][iii] << " " ;}
            streamlog_out(DEBUG3) << std::endl; 
          }  
          nnn=0;
          float   p_cov_v[10];
          for ( int iii=0 ; iii < 4 ; iii++ ) {
            for ( int jjj=0 ; jjj <= iii ; jjj++ ) {
              p_cov_v[nnn] = p_cov[iii][jjj] ; 
              if ( p_cov_v[nnn] != p_cov_v[nnn] ) { streamlog_out(WARNING) << " p_cov_v " << nnn << " is NaN " <<  p_cov_v[nnn] << std::endl ; }
              nnn++ ;
	    }
          }

          part->setCovMatrix(p_cov_v);

        } else {
          if ( part->getCharge() == 0 ) {
	    streamlog_out(DEBUG3) << " atypical neutral : nclu = " << nclu << " ,  ntrk = " << ntrk << std::endl;  
          }
        }
      }   
    }
    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

    streamlog_out(DEBUG4) << "   processing event: " << _nEvt << " which is event " << evt->getEventNumber() 
        << "   in run:  " << evt->getRunNumber() << std::endl ;



    _nEvt ++ ;
}
void AddClusterProperties::debuging(LCCollection* clucol ,ClusterImpl* clu,double* cog,double* cov,double* eval,double* eval_err,double* evp,double* evpe,double* evc,int np,double sum_wgt,double sum_wgtsqr,double sum_wgt4) {
 
          streamlog_out(DEBUG2) << " input: cog                    " << cog[0] << " " << cog[1] << " " << cog[2] <<  std::endl; 
          streamlog_out(DEBUG2) << " input: covariance mat         " << cov[0*3+0] << " "  << cov[0*3+1] << " " << cov[0*3+2] <<   std::endl;  
          streamlog_out(DEBUG2) << " input: covariance mat         " << cov[1*3+0] << " "  << cov[1*3+1] << " " << cov[1*3+2] <<   std::endl;  
          streamlog_out(DEBUG2) << " input: covariance mat         " << cov[2*3+0] << " "  << cov[2*3+1] << " " << cov[2*3+2] <<   std::endl;  
          streamlog_out(DEBUG2) << " input: eigenvals              " << eval[0] << " "  << eval[1] << " "  << eval[2] << std::endl; 
          streamlog_out(DEBUG2) << " input: V(eigenvals)           " << eval_err[0] << " "  << eval_err[1] << " "  << eval_err[2] << std::endl; 
          streamlog_out(DEBUG2) << " input: eigen vec 1 (p)        " << evp[0*3+0] << " "  << evp[1*3+0] <<  std::endl; 
          streamlog_out(DEBUG2) << " input: eigen vec 2 (p)        " << evp[0*3+1] << " "  << evp[1*3+1] <<  std::endl; 
          streamlog_out(DEBUG2) << " input: eigen vec 3 (p)        " << evp[0*3+2] << " "  << evp[1*3+2] <<  std::endl; 
          streamlog_out(DEBUG2) << " input: cov(eigen vec) 1 (p)   " << evpe[0*2*3+0*3+0] << " "  << evpe[1*2*3+1*3+0] << " "  << evpe[1*2*3+0*3+0] <<  std::endl;
          streamlog_out(DEBUG2) << " input: cov(eigen vec) 2 (p)   " << evpe[0*2*3+0*3+1] << " "  << evpe[1*2*3+1*3+1] << " "  << evpe[1*2*3+0*3+1] <<  std::endl;
          streamlog_out(DEBUG2) << " input: cov(eigen vec) 3 (p)   " << evpe[0*2*3+0*3+2] << " "  << evpe[1*2*3+1*3+2] << " "  << evpe[1*2*3+0*3+2] <<  std::endl;
          streamlog_out(DEBUG2) << " input: eigen vec 1 (c)        " << evc[0*3+0] << " "  << evc[1*3+0] << " "  << evc[1*3+0] << std::endl; 
          streamlog_out(DEBUG2) << " input: eigen vec 2 (c)        " << evc[0*3+1] << " "  << evc[1*3+1] << " "  << evc[1*3+1] << std::endl; 
          streamlog_out(DEBUG2) << " input: eigen vec 3 (c)        " << evc[0*3+2] << " "  << evc[1*3+2] << " "  << evc[1*3+2] << std::endl; 
          streamlog_out(DEBUG2) << " input: E, np, sum, sum_sq, sum_4 " << clu->getEnergy() << " " << np << " " << sum_wgt << " " << sum_wgtsqr << " " << sum_wgt4 <<  std::endl; 


         // check reading back: get values now stored in clu
          const float* cogv = clu->getPosition();
          vector<double> cog_fr_clu ;
          for (int iii=0; iii<3 ; iii++) { cog_fr_clu.push_back(cogv[iii]);}
          FloatVec PositionError = clu->getPositionError();
          StringVec shape_keys ;
          shape_keys = clucol->getParameters().getStringVals( std::string("ClusterShapeParameters"),shape_keys);
          FloatVec wgts_sumv = clu->getShape() ; 
          int npnt = 0 ;
          double wgt_sum= 0.0;
          double wgt_sq_sum= 0.0;
          double wgt_4_sum= 0.0;
          for ( unsigned kkk=0 ; kkk < shape_keys.size() ; kkk++ ) {
	    if ( shape_keys[kkk] == "npoints" )   {  npnt = int(wgts_sumv[kkk])  ; }
	    if ( shape_keys[kkk] == "sum_wgt" )   {  wgt_sum = wgts_sumv[kkk]  ; }
	    if ( shape_keys[kkk] == "sum_wgt^2" ) {  wgt_sq_sum = wgts_sumv[kkk]  ; }
            if ( shape_keys[kkk] == "sum_wgt^4" ) {  wgt_4_sum = wgts_sumv[kkk]  ; }
          }
 
          double seen_covmat[3][3];
          int nnn = 0;
          for ( int kkk=0 ; kkk < 3 ; kkk++ ) {
            for ( int lll=0 ; lll<=kkk ; lll++) {
              seen_covmat[lll][kkk] =  PositionError[nnn]*wgt_sum*wgt_sum/wgt_sq_sum;
  	      seen_covmat[kkk][lll] =  PositionError[nnn]*wgt_sum*wgt_sum/wgt_sq_sum;
              nnn++ ;
            }
          }
          vector<double>  cov_fr_clu ; 
          for( int jjj=0 ; jjj<3 ; jjj++) {for (int iii=0; iii<3 ; iii++) { cov_fr_clu.push_back(seen_covmat[jjj][iii]);}}
	  //          for (int iii=0; iii<6 ; iii++) { cov.push_back(PositionError[iii]);}
          FloatVec DirectionError = clu->getDirectionError();
          vector<double> thphcov;
          for (int iii=0; iii<3 ; iii++) { thphcov.push_back(DirectionError[iii]);}
          // a new WeightedPoints3D object, ctored with the cluster-shape params
          WeightedPoints3D wgtp_readback = WeightedPoints3D(  cog_fr_clu, cov_fr_clu , thphcov ,npnt, wgt_sum , wgt_sq_sum , wgt_4_sum);
          // read back stuff just entered
          double* cog_readback=wgtp_readback.getCentreOfGravity();
          double* covv=wgtp_readback.getCentreOfGravityErrors();
          int     np_readback=wgtp_readback.getNumberOfPoints();
          double  sum_wgt_readback=wgtp_readback.getTotalWeight();
          double  sum_wgtsqr_readback=wgtp_readback.getTotalSquaredWeight();
          double  sum_wgt4_readback=wgtp_readback.getTotalQuarticWeight();
          double  cov_readback[3][3] ; 
          for( int iii=0 ; iii<3 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { cov_readback[iii][jjj]=covv[jjj+iii*3] ; } } ; 
          //for( int iii=0 ; iii<3 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { cov_readback[iii][jjj]=covv[jjj+iii*3]*sum_wgtsqr_readback/(sum_wgt_readback*sum_wgt_readback) ; } } ; 
          // now get things not actually entered, but calculated.  
          double* evpv=wgtp_readback.getEigenVecPolar(); 
          double  evp_readback[2][3] ;
          for( int iii=0 ; iii<2 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { evp_readback[iii][jjj]=evpv[jjj+iii*3] ; } } ;
          double* evcv=wgtp_readback.getEigenVecCartesian(); 
          double  evc_readback[3][3] ;
          for( int iii=0 ; iii<3 ; iii++ ) { for ( int jjj=0 ; jjj<3 ; jjj++ ) { evc_readback[iii][jjj]=evcv[jjj+iii*3] ; } } ;
          double* eval_readback=wgtp_readback.getEigenVal();
          double* eval_err_readback=wgtp_readback.getEigenValErrors();
          // this should be a mix of input and calculated: [][][3] is entry, the others calcuated. The latter should be
          // different from the exact values from the points, since the fourth moments are missing.
          double* evpev=wgtp_readback.getEigenVecPolarErrors();              
          double  evpe_readback[2][2][3] ; 
          for( int iii=0 ; iii<2 ; iii++ ) { for ( int jjj=0 ; jjj<2 ; jjj++ ) {for ( int kkk=0 ; kkk<3 ; kkk++ ) { evpe_readback[iii][jjj][kkk]=evpev[kkk+jjj*3+iii*3*2] ; } } };
          streamlog_out(DEBUG2) << " read back: cog                    " << cog_readback[0] << " " << cog_readback[1] << " " << cog_readback[2] <<  std::endl; 
          streamlog_out(DEBUG2) << " read back: covariance mat         " << cov_readback[0][0] << " "  << cov_readback[0][1] << " " << cov_readback[0][2] <<   std::endl;  
          streamlog_out(DEBUG2) << " read back: covariance mat         " << cov_readback[1][0] << " "  << cov_readback[1][1] << " " << cov_readback[1][2] <<   std::endl;  
          streamlog_out(DEBUG2) << " read back: covariance mat         " << cov_readback[2][0] << " "  << cov_readback[2][1] << " " << cov_readback[2][2] <<   std::endl;  
          streamlog_out(DEBUG2) << " read back: eigenvals              " << eval_readback[0] << " "  << eval_readback[1] << " "  << eval_readback[2] << std::endl; 
          streamlog_out(DEBUG2) << " read back: V(eigenvals)           " << eval_err_readback[0] << " "  << eval_err_readback[1] << " "  << eval_err_readback[2] << std::endl; 
          streamlog_out(DEBUG2) << " read back: eigen vec 1 (p)        " << evp_readback[0][0] << " "  << evp_readback[1][0] <<  std::endl; 
          streamlog_out(DEBUG2) << " read back: eigen vec 2 (p)        " << evp_readback[0][1] << " "  << evp_readback[1][1] <<  std::endl; 
          streamlog_out(DEBUG2) << " read back: eigen vec 3 (p)        " << evp_readback[0][2] << " "  << evp_readback[1][2] <<  std::endl; 
          streamlog_out(DEBUG2) << " read back: cov(eigen vec) 1 (p)   " << evpe_readback[0][0][0] << " "  << evpe_readback[1][1][0] << " "  << evpe_readback[1][0][0] <<  std::endl;
          streamlog_out(DEBUG2) << " read back: cov(eigen vec) 2 (p)   " << evpe_readback[0][0][1] << " "  << evpe_readback[1][1][1] << " "  << evpe_readback[1][0][1] <<  std::endl;
          streamlog_out(DEBUG2) << " read back: cov(eigen vec) 3 (p)   " << evpe_readback[0][0][2] << " "  << evpe_readback[1][1][2] << " "  << evpe_readback[1][0][2] <<  std::endl;
          streamlog_out(DEBUG2) << " read back: eigen vec 1 (c)        " << evc_readback[0][0] << " "  << evc_readback[1][0] << " "  << evc_readback[1][0] << std::endl; 
          streamlog_out(DEBUG2) << " read back: eigen vec 2 (c)        " << evc_readback[0][1] << " "  << evc_readback[1][1] << " "  << evc_readback[1][1] << std::endl; 
          streamlog_out(DEBUG2) << " read back: eigen vec 3 (c)        " << evc_readback[0][2] << " "  << evc_readback[1][2] << " "  << evc_readback[1][2] << std::endl; 
          streamlog_out(DEBUG2) << " read back: E, np, sum, sum_sq, sum_4 " << clu->getEnergy() << " " << np_readback << " " << sum_wgt_readback << " " << sum_wgtsqr_readback << " " << sum_wgt4_readback <<  std::endl; 
}

void AddClusterProperties::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AddClusterProperties::end(){ 

    //   std::cout << "AddClusterProperties::end()  " << name() 
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;

}

