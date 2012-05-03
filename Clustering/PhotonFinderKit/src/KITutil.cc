
#include <algorithm>


#include "KITutil.h"

void CreateAllShits2(LCCollection* colt,CellIDDecoder<CalorimeterHit>& id,vector<Superhit2*>* calo)
{

  unsigned int nelem=colt->getNumberOfElements();
  

  for(unsigned int ll=0;ll<nelem;ll++)
    {     
      CalorimeterHit* ch=dynamic_cast<CalorimeterHit*>(colt->getElementAt(ll)); 
      unsigned int layer  =(unsigned int)id(ch)["K-1"];
      Superhit2 * stmp;

       if( layer<=PGDB[PGdb::ECAL1_BAR].max_lay)
	    stmp= new Superhit2(PGDB[PGdb::ECAL1_BAR].mip_whole,ch);
	  else
	    stmp= new Superhit2(PGDB[PGdb::ECAL2_BAR].mip_whole,ch);

       calo[0].push_back(stmp);
	  
    }
}

void TotalPrecalc2(vector<Superhit2*>* calo,CellIDDecoder<CalorimeterHit>& id,
		  unsigned int nelem, int Ncut)
{

     vector<Superhit2*> tmpsh[9];

     for(unsigned int ll=0;ll<nelem;ll++)
       {
	 CalorimeterHit* chh=(calo[0])[ll]->chit;
	 int stove  =id(chh)["S-1"];
	 int module =id(chh)["M"];
	 int layer  =id(chh)["K-1"];
	 
	 (calo[0])[ll]->S=stove;
	 (calo[0])[ll]->M=module;
	 (calo[0])[ll]->K=layer;
	 (calo[0])[ll]->tip=1;

        if( module !=0 && module!=6)        
	  {
	    if(layer!=0)
	      {
		tmpsh[stove].push_back((calo[0])[ll]);	      
	      }else{
	      
	   double fi=atan2(-(calo[0])[ll]->chit->getPosition()[0],(calo[0])[ll]->chit->getPosition()[1]);
	      if ( fi< 0.0) 
		{
		  if( stove!=0)
		  fi=2.0*M_PI+fi;
		  else
		  fi=fabs(fi);           
		}
	    
	      double diff=fi-M_PI_4*stove;
	      tmpsh[stove].push_back((calo[0])[ll]);
              if ( fabs(diff)<M_PI_4/2.0)
		{
		  tmpsh[stove].push_back((calo[0])[ll]);
		}else{

	        (calo[0])[ll]->connect=true;
		    
		tmpsh[stove].push_back((calo[0])[ll]);
		
		int stovem1=stove-1;
		
		if( stovem1>=0)
                   tmpsh[stovem1].push_back((calo[0])[ll]);
		else
		   tmpsh[7].push_back((calo[0])[ll]);
	        }
	      }                       
	  }else{
	  tmpsh[8].push_back((calo[0])[ll]);
	  }
     }

      

   double celx=PGDB[PGdb::ECAL1_BAR].cell_size;
   double distcut=0.98*(2.0*celx)*(2.0*celx);
  
   for(unsigned int ll=0;ll<8;ll++)
      {
	if(tmpsh[ll].size()!=0)
	  Precalc2(tmpsh[ll],PGDB[PGdb::ECAL1_BAR].r_inner,
		            PGDB[PGdb::ECAL1_CAP].z_inner,
	                     PGDB[PGdb::ECAL1_BAR].cell_size,distcut,true,ll,1,id); 
      }


   if(tmpsh[8].size()!=0)
     Precalc2(tmpsh[8],PGDB[PGdb::ECAL1_BAR].r_inner,
	              PGDB[PGdb::ECAL1_CAP].z_inner,
	              PGDB[PGdb::ECAL1_BAR].cell_size,distcut,false,1,1,id); 



   for(unsigned int ii=0;ii<calo[0].size();ii++)
       calo[0][ii]->top=calo[0][ii]->neighbours.size();

 
     for(unsigned int ii=0;ii<(calo[0]).size();ii++)
       { 
        if ((calo[0])[ii]->top !=0)
	  {
	     if((calo[0])[ii]->top >= 1 && (calo[0])[ii]->top <Ncut)
	       {
	         (calo[2]).push_back((calo[0])[ii]);  
	       }else{
	         (calo[4]).push_back((calo[0])[ii]);  
	       }
          }else{  // zero or not 
	     (calo[6]).push_back((calo[0])[ii]);
	  }
       }
}













void Precalc2(vector< Superhit2* >& shvec,double r, double z, double cell, double dist,bool tip,int stove,int module,CellIDDecoder<CalorimeterHit>& id){

           unsigned int nPts=shvec.size();
           ANNpointArray dataPts;
           ANNpoint queryPt;
           ANNidxArray nnIdx;
           ANNdistArray dists;
           ANNkd_tree* kdTree;

	   queryPt=annAllocPt(3); // 3d point 
	   dataPts=annAllocPts(nPts,3);
           nnIdx=new ANNidx[36];
	   dists=new ANNdist[36];
	   float xxx[3];
           float Ecalradius=r;
	   float Ecalhalfz=z;
	   float Ecalcellsize=cell;

	   for ( unsigned int ii=0; ii<nPts;ii++)
	     {  
	       GridTransform2(shvec[ii]->chit,Ecalradius,Ecalhalfz,Ecalcellsize,xxx,tip,stove,module,id);
	       dataPts[ii][0]=xxx[0];  
	       dataPts[ii][1]=xxx[1];  
	       dataPts[ii][2]=xxx[2];  
	       shvec[ii]->point[0]=xxx[0];
	       shvec[ii]->point[1]=xxx[1];
	       shvec[ii]->point[2]=xxx[2];
	     }


	   kdTree=new ANNkd_tree(dataPts,nPts,3);                                                       
	   ANNdist trd=dist;
	  
           for (unsigned int ii=0; ii<nPts;ii++)
	     {  	      
	       unsigned int koliko = (unsigned int)kdTree->annkFRSearch(dataPts[ii],trd,26,nnIdx,dists,0.0);
	      if (koliko>26) 
 	          koliko=26;
	    if( tip )  
	      {
	       if( shvec[ii]->connect )  
		 {
		   
		  for(unsigned int ij=0;ij<koliko;ij++)
		   {
		     if( shvec[ii]->neighbours.end() == 
		       find(shvec[ii]->neighbours.begin(),shvec[ii]->neighbours.end(), shvec[nnIdx[ij]]))
		      {
			shvec[ii]->neighbours.push_back(shvec[nnIdx[ij]]);
			shvec[ii]->top++;
		
		      }
	             if( shvec[nnIdx[ij]]->neighbours.end() == 
	        find(shvec[nnIdx[ij]]->neighbours.begin(),shvec[nnIdx[ij]]->neighbours.end(), shvec[ii]))
		      {		
			shvec[nnIdx[ij]]->neighbours.push_back(shvec[ii]);
			shvec[nnIdx[ij]]->top++;
 		      }

                   }

		 }else{
		 if( koliko !=0)
		   {
		     shvec[ii]->top=koliko;		   
		     for(unsigned int ij=0;ij<koliko ;ij++)
		       shvec[ii]->neighbours.push_back( shvec[nnIdx[ij]]);
		   }
	         }
	      }else{
	      if( koliko !=0)
		   {
	         shvec[ii]->top=koliko;
		for(unsigned int ij=0;ij<koliko;ij++)
		 shvec[ii]->neighbours.push_back( shvec[nnIdx[ij]]);
		   }	        
	      }
	     }
	   delete [] nnIdx;			     
	   delete [] dists;
	   annDeallocPts(dataPts);
	   delete kdTree;
	   annClose();

}



void GridTransform2( CalorimeterHit* clh,float& radius, float& halfz, float& cellsize,float*X,
		   bool tip,int stove,int module,CellIDDecoder<CalorimeterHit>& id)
{

            unsigned int layer;
	    double coss;
	    double sinn;

	     if( !tip)
		{
		   stove  = id(clh)["S-1"];
		   module = id(clh)["M"];
		   layer  = id(clh)["K-1"];
      
		}else{
		int stove2= id(clh)["S-1"];               
		if ( stove2!=stove)
		  {

		    coss= cos(M_PI_4*stove);
		    sinn= sin(M_PI_4*stove);
		    float tmp_y=coss*clh->getPosition()[1]-sinn*clh->getPosition()[0];
		    layer=(unsigned int) ((tmp_y-radius)/PGDB[PGdb::ECAL1_BAR].sampling);

                  if ( layer> PGDB[PGdb::ECAL1_BAR].n_sampl-1 ) 
		    {
	              layer=PGDB[PGdb::ECAL1_BAR].n_sampl-1+ 
			(unsigned int)(( tmp_y-radius- PGDB[PGdb::ECAL1_BAR].n_sampl* 
				PGDB[PGdb::ECAL1_BAR].sampling)/PGDB[PGdb::ECAL2_BAR].sampling); 
		    }
                 
		  }else{
		   layer = id(clh)["K-1"]; 
		  }

	        }
	   
	            coss= cos(M_PI_4*stove);
		    sinn= sin(M_PI_4*stove);

	       if( module !=0 && module!=6)
		 {
                  float tmp_x=coss*clh->getPosition()[0]+sinn*clh->getPosition()[1];
		  float tmp_y=radius+(layer+1)*cellsize;
		 
		  X[0]= coss*tmp_x-sinn*tmp_y;
		  X[1]= coss*tmp_y+sinn*tmp_x;
		  X[2]=clh->getPosition()[2];
		 
                 }else{ // endcap
	           X[0]=clh->getPosition()[0];
		   X[1]=clh->getPosition()[1];
		   X[2]= (-1+module/3)*(halfz+layer*cellsize);
	        }

}


void FindCores2(Shitvec2* secal1, vector<Tmpclvec2>& bbb , vector <PROTSEED2> * prs2,
		unsigned int N, vector<float> miipstep, CoreCut2 Ccut)
{       
          typedef pair<int,int>   test;
	  int psl_plun_global=0;
	  
	  double Diagby2=PGDB[PGdb::ECAL1_BAR].cell_size*sqrt(2.0)/2.0;
	  double d=PGDB[PGdb::ECAL1_BAR].r_inner;
	  double z=PGDB[PGdb::ECAL1_CAP].z_inner;
	  double Xatf[3];
	  vector< vector<Superhit2*> > blk(N) ;
	 
	  for( unsigned int iim=0;iim<secal1->size();iim++)
	    {	      
	      int koji=0;

	      for(unsigned int ib=0;ib<N;ib++)
		{
		  if( (*secal1)[iim]->mip > miipstep[ib])
		      koji=ib;
		  else 
		    break;
	         }  
	      for( int im=koji;im>=0 ;im--)
	        blk[im].push_back((*secal1)[iim]);
	    }

	  vector<Tmpcl2*>  testvektor;

	  for( unsigned int iim=0;iim<N;iim++)
	    { 
	      cluster5(&blk[iim],&bbb[iim]); 
	    }

	  for(unsigned int iim=0;iim<N;iim++)
	    {

	      if(bbb[iim].size()!=0)
		{
		  for(unsigned int im=0;im<bbb[iim].size();im++)
		    {
		      bbb[iim][im]->calcCenter();
		      bbb[iim][im]->findInertia();
		     
		    }
		}
	    }
	
	  int poslednjipun=N-1;  
	  for( unsigned int im=0;im<N;im++)
	    {
	      if(bbb[im].size()==0)
		{
		  poslednjipun=im-1;
		  break;
		}
	      poslednjipun=im;
	    }
	  psl_plun_global=poslednjipun;
	 
	  vector < PROTSEED2 > prs;

		for(unsigned int im=0;im<bbb[0].size();im++)
		  {      
		    if(bbb[0][im]->hits.size()>Ccut.MinHit0)
		      {
			PROTSEED2 km;
			km.cl=bbb[0][im];
			km.level=0;
			km.X[0]=0.0;
			km.X[1]=0.0;
			km.X[2]=0.0;
			km.active=true;
			prs.push_back(km);
		      }
		  }
	       
	    
	      for( unsigned int iim=0;iim<prs.size();iim++)
		{ 		 
		  Xatf[0]=0.0;
		  Xatf[1]=0.0;
		  Xatf[2]=0.0;
				
		  LineCaloIntersectD2(prs[iim].cl->getCenter(), prs[iim].cl->direction,d,z,Xatf);

		  prs[iim].X[0]=Xatf[0];
		  prs[iim].X[1]=Xatf[1];
		  prs[iim].X[2]=Xatf[2];
		  prs2->push_back(prs[iim]);
		 
		}
 	        
	  for(int i=1;i<=poslednjipun;i++)
	   {
	   
 	      for( unsigned int im=0;im<prs2->size();im++)
 		{   
		  if( (*prs2)[im].active)
		    {
 		       Tmpclvec2 izv;
		       ClusterInCluster2((*prs2)[im].cl, bbb[i],izv);

		  if( izv.size()<=1) 
		    {
		      // 1-in-1 nothing to do 
		    }else{
		    		  
		    vector < int > gssf;
		    vector < int > gssfa;
		    for( unsigned int imj=0;imj<izv.size();imj++)
		        {			
				  
			  for(unsigned int imm=0;imm<izv.size();imm++)
			    {
			      if(imm!=imj)
				{
				  
				  if( izv[imj]->hits.size()>=Ccut.MinHitSplit 
				   && izv[imm]->hits.size()>=Ccut.MinHitSplit)
				  {
				    double XIP[3];
				    XIP[0]=0.0;  XIP[1]=0.0;  XIP[2]=0.0;

				  // Split the split in two categories 
				  // clear one (larger distance) and not so clear
				  
				  if((LinePointDistance2(XIP,(*prs2)[im].cl->center,izv[imm]->center)>Diagby2 || 
				      LinePointDistance2(XIP,(*prs2)[im].cl->center,izv[imj]->center)>Diagby2 ))
				    {
				       if(find(gssf.begin(),gssf.end(),(int)imm)==gssf.end())
					 gssf.push_back(imm);
				       if(find(gssf.begin(),gssf.end(),(int)imj)==gssf.end())
					 gssf.push_back(imj);

				    }else{
				         if(find(gssfa.begin(),gssfa.end(),(int)imm)==gssfa.end())
					    gssfa.push_back(imm);
				         if(find(gssfa.begin(),gssfa.end(),(int)imj)==gssfa.end())
					    gssfa.push_back(imj);
				    }
				   }
				}
			    }
		         }

		    if(gssf.size()!=0) 
		      {
					    
		      int nh1=0;
		      
		      int nh3=(*prs2)[im].cl->hits.size();
		      for( unsigned int j=0;j<gssf.size();j++)
			   nh1+=izv[gssf[j]]->hits.size();
		   
		      double ratio = ((double)nh1)/((double)nh3);
		      (*prs2)[im].active=false;

		      if( ratio > Ccut.Rcut )
			{
			  for( unsigned int j=0;j<gssfa.size();j++)
			   {			     
			      Xatf[0]=0.0;
			      Xatf[1]=0.0;
			      Xatf[2]=0.0;

			      PROTSEED2 a;
		
			      a.cl=izv[gssfa[j]];
			      a.level=i;
			      LineCaloIntersectD2(izv[gssfa[j]]->getCenter(), 
						  izv[gssfa[j]]->direction,d,z,Xatf);
			      a.X[0]=Xatf[0];
			      a.X[1]=Xatf[1];
			      a.X[2]=Xatf[2];
			      a.active=true;
			      bool stavi=true;
			      for( unsigned int im=0;im<prs2->size();im++)
				{
				  if((*prs2)[im].active)
				    {
				      if( (*prs2)[im].cl==a.cl)
					stavi=false;
				    }
				}
			      if(stavi)
				prs2->push_back(a);
			      
			   }
			}
		      
		      
                      for( unsigned int j=0;j<gssf.size();j++)
			{
			  
			  Xatf[0]=0.0;
			  Xatf[1]=0.0;
			  Xatf[2]=0.0;

			  PROTSEED2 a;
			 
			  a.cl=izv[gssf[j]];
			  a.level=i;
			  LineCaloIntersectD2(izv[gssf[j]]->getCenter(), 
					     izv[gssf[j]]->direction,d,z,Xatf);
			  a.X[0]=Xatf[0];
			  a.X[1]=Xatf[1];
			  a.X[2]=Xatf[2];
			  a.active=true;
		
			  bool stavi=true;
			  for( unsigned int im=0;im<prs2->size();im++)
			    {
			      if((*prs2)[im].active)
				{
				  if( (*prs2)[im].cl==a.cl)
				    stavi=false;
				}
			    }
			  if(stavi)
			  prs2->push_back(a);
			}
                     
		      
		      } 
		  
		    }
		    }// if active 
 		}
	      
	   }


	 vector <test> kojsv;

	  for( unsigned int im=0;im<prs2->size();im++)
 	    {
	      if((*prs2)[im].active)
		    {		      

		    double Xa[3];
		    ModuleNormal2((*prs2)[im].cl->center,z, Xa);
	

		      for( unsigned int ij=0;ij<prs2->size();ij++)
			{
			  if( ij!=im && (*prs2)[ij].active)
			    {
			     
			      double dircc[3];
			      dircc[0]=(*prs2)[im].cl->center[0]-(*prs2)[ij].cl->center[0];
			      dircc[1]=(*prs2)[im].cl->center[1]-(*prs2)[ij].cl->center[1];
			      dircc[2]=(*prs2)[im].cl->center[2]-(*prs2)[ij].cl->center[2];    
			    
			      //MERGE CONDITION
			      if ( D_cl_cl2((*prs2)[im].cl,(*prs2)[ij].cl)<Ccut.Distcut
				   && fabs( Dot2((*prs2)[ij].cl->center,dircc))>Ccut.Coscut)
				 {
				   test a,b;
				   a.first=im;
				   a.second=ij;
				   b.first=ij;
				   b.second=im;
				   if(find(kojsv.begin(),kojsv.end(),a)==kojsv.end() && 
				      find(kojsv.begin(),kojsv.end(),b)==kojsv.end())
				     {
				     kojsv.push_back(a);
				  
				     }
				 }
			      
			    }
			}
		    }
	    }


	  if(kojsv.size()!=0)
	    {
           
	      int iskoristio[kojsv.size()];
	      for(unsigned int i=0;i<kojsv.size();i++)
		{
		  iskoristio[i]=0;
		}
	      for(unsigned int i=0;i<kojsv.size();i++)
		{
		  if( iskoristio[i]==0)
		    {
		  vector <int> trojka;
		  trojka.push_back(kojsv[i].first);
		  trojka.push_back(kojsv[i].second);
		  iskoristio[i]=1;
                  for(unsigned int j=0;j<kojsv.size();j++)
		    {
		     
		      if(
			 (find(trojka.begin(),trojka.end(),kojsv[j].first)!=trojka.end() && 
			  find(trojka.begin(),trojka.end(),kojsv[j].second)==trojka.end()  ) || 
			 (find(trojka.begin(),trojka.end(),kojsv[j].second)!=trojka.end() && 
			  find(trojka.begin(),trojka.end(),kojsv[j].first)==trojka.end()  ))
			{
			
			  if(find(trojka.begin(),trojka.end(),kojsv[j].first)==trojka.end())
			    {
			      if(iskoristio[j]!=1)
				{ 
				 
				  iskoristio[j]=1;
				  if( (*prs2)[kojsv[j].first].active)
 				  trojka.push_back(kojsv[j].first);
				}

			    } else{
			       if(iskoristio[j]!=1)
				{ 
				  if( (*prs2)[kojsv[j].second].active)
				     trojka.push_back(kojsv[j].second);
				  iskoristio[j]=1;
				}
			    }  
			} 
		    }
	
		 
                  int kakvi=1;
		  for( unsigned int  j=1;j<trojka.size();j++)
		    {
		      if((*prs2)[trojka[0]].level!=(*prs2)[trojka[j]].level)
			kakvi++;
		    }
		 
	   if(kakvi==1)
	     {
		    
		   Tmpcl2* cl = new Tmpcl2();
		        for(unsigned int j=0;j<trojka.size();j++)
			  {
			    if( (*prs2)[trojka[j]].active)
			      {
				for(unsigned int jj=0;jj<(*prs2)[trojka[j]].cl->hits.size();jj++)
				  {	
				    cl->hits.push_back((*prs2)[trojka[j]].cl->hits[jj]);
				  }
				(*prs2)[trojka[j]].active=false;
			      }
			  }
		     if( cl->hits.size()!=0)
		       {
			
		       cl->calcCenter();
		       cl->findInertia();
		       bbb[(*prs2)[trojka[0]].level].push_back(cl);
		     
			Xatf[0]=0.0;
			Xatf[1]=0.0;
			Xatf[2]=0.0;
		
			PROTSEED2 a; 
			a.cl=cl;
			a.level=(*prs2)[trojka[0]].level;
			LineCaloIntersectD2(cl->getCenter(), 
					   cl->direction,d,z,Xatf);
			a.X[0]=Xatf[0];
			a.X[1]=Xatf[1];
			a.X[2]=Xatf[2];
			a.active=true;
		
			bool stavi=true;
			  for( unsigned int im=0;im<prs2->size();im++)
			    {
			      if((*prs2)[im].active)
				{
				  if( (*prs2)[im].cl==a.cl)
				    stavi=false;
				}
			    }
			  if(stavi)
			  prs2->push_back(a);
			    
		       }else{
		       delete cl;
		       }
		    }else{
		   
                       
		    vector <int > nivoi;
		    for(unsigned int j=0;j<trojka.size();j++)
		        {
			 if( find(nivoi.begin(),nivoi.end(),(*prs2)[trojka[j]].level)==nivoi.end())
			   nivoi.push_back((*prs2)[trojka[j]].level);
		        }

		     if(trojka.size()==2) 
		       {
			 if( (*prs2)[trojka[0]].cl->getEnergy() >=(*prs2)[trojka[1]].cl->getEnergy())
			   {
			     (*prs2)[trojka[1]].active=false;
			   }else{
			     (*prs2)[trojka[0]].active=false;
		           }
		       }else{
                          
		       vector <int> izabrani;
		       int isko[trojka.size()];
		       for(unsigned int jj=0;jj<trojka.size();jj++)
			     {
			       isko[jj]=0;
			     }
		       for(unsigned int j=0;j<nivoi.size();j++)   
			 {
			       vector<int> sabr;
			       for(unsigned int jj=0;jj<trojka.size();jj++)
				 {
				   if(isko[jj]==0)
				     {
				       if( nivoi[j]==(*prs2)[trojka[jj]].level)
					 {
					   sabr.push_back(trojka[jj]);
					   isko[jj]=1;
					 }
				 } 
			 }

			if(sabr.size()>1)
			  {

			    Tmpcl2* cl = new Tmpcl2();
			         for(unsigned int jk=0;jk<sabr.size();jk++)
				   {
				       (*prs2)[sabr[jk]].active=false;
				     for(unsigned int jj=0;jj<(*prs2)[sabr[jk]].cl->hits.size();jj++)
				       cl->hits.push_back((*prs2)[sabr[jk]].cl->hits[jj]);
				    } 

				  cl->calcCenter();
				  cl->findInertia();
				  bbb[nivoi[j]].push_back(cl);
				
				  Xatf[0]=0.0;
				  Xatf[1]=0.0;
				  Xatf[2]=0.0;
				 
				  PROTSEED2 a;
		    
				  a.cl=cl;
				  a.level=(*prs2)[nivoi[j]].level;
				  LineCaloIntersectD2(cl->getCenter(), 
					   cl->direction,d,z,Xatf);
				  a.X[0]=Xatf[0];
				  a.X[1]=Xatf[1];
				  a.X[2]=Xatf[2];
				  a.active=true;
			
				  bool stavi=true;
				  for( unsigned int im=0;im<prs2->size();im++)
				    {
				      if((*prs2)[im].active)
					{
					  if( (*prs2)[im].cl==a.cl)
					    stavi=false;
					}
				    }
				  if(stavi)
				    {
				    prs2->push_back(a);
				    izabrani.push_back(prs2->size()-1);
				    }else{
				      delete cl;
				    }
			     }else{
			     izabrani.push_back(sabr[0]);
			     }			   
			 }
		       if( izabrani.size()>1)
			 {
			   typedef pair <double,int> DpI;
			   vector <DpI> bfjl;
			   for( unsigned int j=0;j<izabrani.size();j++)
			     {
			       DpI tmp;
			       tmp.first=(*prs2)[izabrani[j]].cl->getEnergy();			    
			       tmp.second=izabrani[j];
			       bfjl.push_back(tmp);
			     }
			  
			     sort(bfjl.begin(),bfjl.end());

		             for( unsigned int j=0;j<izabrani.size()-1;j++)
			     {
			       (*prs2)[bfjl[j].second].active=false;
			     }
			    
			    
			 }       
		     }
  
		    }
		  }   
	    } // koj =0.0	   
        }
}




void cluster5( vector<Superhit2*>* shv, vector<Tmpcl2*>* clv) 
{
  
 if( (*shv).size()!=0)
  {

    for(unsigned int i=0;i<(*shv).size();i++)
       {	 
	 (*shv)[i]->is_assigned=false;
       }

    unsigned int shvsz=(*shv).size();
 
    for(unsigned int i=0;i<shvsz;i++)  
      {
	
	   if( (*shv)[i]->top==0 )
	      {
		Tmpcl2* cl =new Tmpcl2();
		cl->hits.push_back((*shv)[i]);
		(*shv)[i]->is_assigned=true;
		clv->push_back(cl);
		continue;
	      }else{

	     if( (*shv)[i]->is_assigned==false)
	       {
		 vector<Superhit2*> sshv;
		 stack <Superhit2*> shstack;
		 sshv.push_back((*shv)[i]);
		 
		 for(unsigned int j=0;j<(*shv)[i]->neighbours.size();j++)
		   {   
		     if( !((*shv)[i]->neighbours[j]->is_assigned) && 
			 shv->end() !=find(shv->begin(),shv->end(), (*shv)[i]->neighbours[j]))
		       {
		     sshv.push_back((*shv)[i]->neighbours[j]);
		     (*shv)[i]->neighbours[j]->is_assigned=true;
		     shstack.push((*shv)[i]->neighbours[j]);
		       }
		   }

		 while(!shstack.empty())
		     {
		       Superhit2* sh=shstack.top();
		       shstack.pop();
		       //sshv.push_back(sh);
	
		       for(unsigned int j=0;j<sh->neighbours.size();j++)
			   {  
			     if( sh->neighbours[j]!=0 )
			       {
			     if( ! sh->neighbours[j]->is_assigned &&   				 
				 shv->end() !=find(shv->begin(),shv->end(), sh->neighbours[j]))
			       {
			       if(sshv.end()==find(sshv.begin(),sshv.end(),sh->neighbours[j])) 
				 {
				   sshv.push_back(sh->neighbours[j]);
				   sh->neighbours[j]->is_assigned=true;
				 
				   shstack.push(sh->neighbours[j]);
				 }
			        }
			       }
			   }  
		     }
		 if(sshv.size()!=0)
		   {
		     Tmpcl2* cl=new Tmpcl2();
		     for(unsigned int im=0;im<sshv.size();im++)
		       {
			 cl->hits.push_back(sshv[im]);
		       }
		     clv->push_back(cl);
		   }
	       }    
              }

      }// over all hits

  for(unsigned int i=0;i<shv->size();i++)
     if((*shv)[i]->is_assigned==false)
       {
	 Tmpcl2* cl=new Tmpcl2();
	 cl->hits.push_back((*shv)[i]);
	 clv->push_back(cl);
       }

  }
}


Superhit2::~Superhit2(){}

Superhit2::Superhit2(float E, CalorimeterHit* chitin){
  chit=chitin;
   for(unsigned int i=0;i<3;i++)
    {
     point[i]=0.0;   
     pointt[i]=chitin->getPosition()[i];
    }
   mip=chit->getEnergy()/E;
   connect=false;
   is_assigned=false;
   mipE=E;
   top=0;
   cl=0; 
   S=0;
   M=0;
   K=0;
   tip=0;
     
}
double* Tmpcl2:: getCenter()
{
  return center;
}


Tmpcl2::Tmpcl2(){
  energy=0.0;
}
Tmpcl2::~Tmpcl2(){}

double Tmpcl2::getEnergy(){

      energy =0.0;
     for( unsigned int i=0; i<hits.size();i++)
     {
      energy+=hits[i]->chit->getEnergy();
     }
      return energy;
}
void Tmpcl2::calcCenter()
{
  center[0]=0.0;
  center[1]=0.0;
  center[2]=0.0;
  double norm=0.0;
  double en=0.0;
  for( unsigned int i=0; i<hits.size();i++)
    {
     en=hits[i]->chit->getEnergy();
     center[0]+=hits[i]->pointt[0]*en;
     center[1]+=hits[i]->pointt[1]*en;
     center[2]+=hits[i]->pointt[2]*en;
     norm+=en;
    } 
  center[0]=center[0]/norm;
  center[1]=center[1]/norm;
  center[2]=center[2]/norm;
}



void Tmpcl2::findInertia()
{

  double aIne[3][3];  
  for (int i(0); i < 3; ++i) {
      for (int j(0); j < 3; ++j) {
	  aIne[i][j] = 0.0;
      }
  }

 
 for( unsigned int i=0; i<hits.size();i++)
    {
      double dX = hits[i]->pointt[0] - center[0];
      double dY = hits[i]->pointt[1] - center[1];
      double dZ = hits[i]->pointt[2] - center[2];
      double eneh=hits[i]->chit->getEnergy();
      eneh=1.0;
      aIne[0][0] += eneh*(dY*dY+dZ*dZ);
      aIne[1][1] += eneh*(dX*dX+dZ*dZ);
      aIne[2][2] += eneh*(dX*dX+dY*dY);
      aIne[0][1] -= eneh*dX*dY;
      aIne[0][2] -= eneh*dX*dZ;
      aIne[1][2] -= eneh*dY*dZ;
  }

  for (int i(0); i < 2; ++i) {
      for (int j = i+1; j < 3; ++j) {
	  aIne[j][i] = aIne[i][j];
      }
  }
 

  gsl_matrix_view aMatrix = gsl_matrix_view_array((double*)aIne,3,3);
  gsl_vector* aVector = gsl_vector_alloc(3);
  gsl_matrix* aEigenVec = gsl_matrix_alloc(3,3);
  gsl_eigen_symmv_workspace* wa = gsl_eigen_symmv_alloc(3);
  gsl_eigen_symmv(&aMatrix.matrix,aVector,aEigenVec,wa);
  gsl_eigen_symmv_free(wa);
  gsl_eigen_symmv_sort(aVector,aEigenVec,GSL_EIGEN_SORT_ABS_ASC);
 double norm=0.0;
  for (int i(0); i < 3; i++) {
    inteigen[i] = gsl_vector_get(aVector,i);
    norm+=inteigen[i]* inteigen[i];
    for (int j(0); j < 3; j++) {
      inteigenvec[i+3*j] = gsl_matrix_get(aEigenVec,i,j);
    }
  }
   norm=sqrt(norm);
 

  for (int i(0); i < 3; i++) 
    inteigen[i]=inteigen[i]/norm;

  direction[0]=inteigenvec[0];
  direction[1]=inteigenvec[1];
  direction[2]=inteigenvec[2];

  gsl_vector_free(aVector);
  gsl_matrix_free(aEigenVec);


}

void LineCaloIntersectD2( double* X1, double* dir,double&d,double&zmax, double*X)
{
  
  double X2[3];
  X2[0]=X1[0]+dir[0]*400.0;
  X2[1]=X1[1]+dir[1]*400.0;
  X2[2]=X1[2]+dir[2]*400.0;

  LineCaloIntersect2( X2,X1,d,zmax,X);

  if( X[0]==0.0) 
    {

      X2[0]=X1[0]-dir[0]*400.0;
      X2[1]=X1[1]-dir[1]*400.0;
      X2[2]=X1[2]-dir[2]*400.0;

      LineCaloIntersect2(X2, X1,d,zmax,X);
      
    }
}

void LineCaloIntersect2(double* X1, double* X2,double&d,double&zmax,  double*X)
{
        double  n[10][3];

	double koren2=M_SQRT2/2.0;
   
	n[0][0]=0.0; 
	n[0][1]=1.0;
	n[0][2]=0.0;

	n[1][0]=koren2; 
	n[1][1]=koren2;
	n[1][2]=0.0;

	n[2][0]=1.0; 
	n[2][1]=0.0;
	n[2][2]=0.0;

	n[3][0]=-1.0; 
	n[3][1]=0.0;
	n[3][2]=0.0;

	n[4][0]=0.0; 
	n[4][1]=-1.0;
	n[4][2]=0.0;

	n[5][0]=-koren2; 
	n[5][1]=-koren2;
	n[5][2]=0.0;
 
	n[6][0]=koren2; 
	n[6][1]=-koren2;
	n[6][2]=0.0;

	n[7][0]=-koren2; 
	n[7][1]=koren2;
	n[7][2]=0.0;

	
	n[8][0]=0.0; 
	n[8][1]=0.0;
	n[8][2]=1.0;

	n[9][0]=0.0; 
	n[9][1]=0.0;
	n[9][2]=-1.0;
 
 

  if ( abs(X2[2])<zmax)
   {

    
         vector < vector <double> > test;
	 for ( unsigned int i=0; i<8;i++)
	   {
	     double tmp=-(n[i][0]*n[i][0]+n[i][1]*n[i][1]+n[i][2]*n[i][2])*d
	                +n[i][0]*X1[0]+n[i][1]*X1[1]+n[i][2]*X1[2];

	     double tmp2=n[i][0]*(X1[0]-X2[0])+n[i][1]*(X1[1]-X2[1])+n[i][2]*(X1[2]-X2[2]);
	    
	     double t=-2.0;     
	     if(tmp2!=0.0)
	       {
		 t=tmp/tmp2;
	
	       }

	     if (t>0.0 && t<=1.0)
	       {
		 vector <double> tmpi;
		 tmpi.push_back(X1[0]+(X2[0]-X1[0])*t);
		 tmpi.push_back(X1[1]+(X2[1]-X1[1])*t);
		 tmpi.push_back(X1[2]+(X2[2]-X1[2])*t);
		 test.push_back(tmpi);
	       }

	   }
             
  if (test.size()==1)
      {
	 
      X[0]=test[0][0];
      X[1]=test[0][1];
      X[2]=test[0][2];
      double rts=sqrt(X[0]*X[0]+X[1]*X[1]);
      if( rts< ((d+2.0)*1.0823922002924))
	{
      return;
	}else{
	 X[0]=0.0;
	 X[1]=0.0;
	 X[2]=0.0;
	 return;
        } 
      }else{

	double rmin=1e+22;
	unsigned int imin=999;
  
	for( unsigned int i=0; i<test.size();i++)
	  {
	    double Xa[3];
	    Xa[0]=test[i][0]; Xa[1]=test[i][1]; Xa[2]=test[i][2];

	    double tmpr=sqrt(Xa[0]*Xa[0]+Xa[1]*Xa[1]);
	    if( tmpr<rmin)
	      {
		rmin=tmpr;
		imin=i;
	      }
	  }

	if ( imin< test.size())
	  {
	    X[0]=test[imin][0];
	    X[1]=test[imin][1];
	    X[2]=test[imin][2];
	    return;
	  }else{
	    X[0]=0.0;
	    X[1]=0.0;
	    X[2]=0.0;
	    return;
	  }
      }
   
    }else { // endcap 
                double n[8][3];
                if(X2[2]>0.0)
                 {
		   n[0][0]=0.0; 
		   n[0][1]=0.0;
		   n[0][2]=1.0;
		 }else{
		  n[0][0]=0.0; 
		  n[0][1]=0.0;
		  n[0][2]=-1.0;
		}


		double tmp=-(n[0][0]*n[0][0]+n[0][1]*n[0][1]+n[0][2]*n[0][2])*zmax
		            +n[0][0]*X1[0]+n[0][1]*X1[1]+n[0][2]*X1[2];
		double tmp2=n[0][0]*(X1[0]-X2[0])+n[0][1]*(X1[1]-X2[1])+n[0][2]*(X1[2]-X2[2]);
		double t=0.0;
     
		if(tmp2!=0.0)
		  {
		    t=tmp/tmp2;
		  }
		if(t>0.0 && t<=1.0)
		  {
		    X[0]=X1[0]+(X2[0]-X1[0])*t ;
		    X[1]=X1[1]+(X2[1]-X1[1])*t ;
		    X[2]=X1[2]+(X2[2]-X1[2])*t ; 
		    return;
		  }else{
		   X[0]=0.0;
		   X[1]=0.0;
		   X[2]=0.0;
		   return;
		  }
   }
}


void ClusterInCluster2(Tmpcl2* cl, vector<Tmpcl2*>& clv)
{
 
   sort(cl->hits.begin(),cl->hits.end());

   for(unsigned int i=0;i<clv.size();i++)
      {
	sort(clv[i]->hits.begin(),clv[i]->hits.end());
        if( includes(cl->hits.begin(),cl->hits.end(),clv[i]->hits.begin(),clv[i]->hits.end()))
	  { 
	     cl->daughters.push_back(clv[i]);
	     clv[i]->parents.push_back(cl);
	  }
     }

}

double LinePointDistance2( double* X1, double* X2, double* X0){

double tmp1=X1[1]*X2[0]-X1[0]*X2[1]-X1[1]*X0[0]+X2[1]*X0[0]+X1[0]*X0[1]-X2[0]*X0[1]; 
double tmp2=X1[2]*X2[0]-X1[0]*X2[2]-X1[2]*X0[0]+X2[2]*X0[0]+X1[0]*X0[2]-X2[0]*X0[2]; 
double tmp3=X1[2]*X2[1]-X1[1]*X2[2]-X1[2]*X0[1]+X2[2]*X0[1]+X1[1]*X0[2]-X2[1]*X0[2]; 
double tmp4=(X1[0]-X2[0])*(X1[0]-X2[0])+(X1[1]-X2[1])*(X1[1]-X2[1])+(X1[2]-X2[2])*(X1[2]-X2[2]);
 
 double tmp5=tmp1*tmp1+tmp2*tmp2+tmp3*tmp3;
 tmp5=tmp5/tmp4;
 return sqrt(tmp5);

}

void ModuleNormal2(double* X1,double& zmax, double* X0)
{


  if ( abs(X1[2])<zmax)
    {
  
    double fi=atan2(X1[0],X1[1]);
    double fi0=M_PI_4/2.0;
   
    int i;
    if( fi>=0.)
      i=(int)((fi+fi0)/M_PI_4);
    else
      i=(int)((fi-fi0)/M_PI_4);


    double fip=i*M_PI_4;
    if( i!=0)
      {
	X0[0]=sin(fip);
	X0[1]=cos(fip);
	X0[2]=0.0;
      }else{
      X0[0]=0.0;
      X0[1]=1.0;
      X0[2]=0.0;

      }
    }else{

    if( X1[2]>0.)
      {
	X0[0]=0.0;
	X0[1]=0.0;
	X0[2]=1.0;
      }else{
	X0[0]=0.0;
	X0[1]=0.0;
	X0[2]=-1.0;
      }
  }
}
double D_cl_cl2(Tmpcl2* cl1,Tmpcl2* cl2) 
{

 vector <double> dist;
 
 
    for( unsigned int i=0; i<cl1->hits.size();i++)
    {
      double x1=cl1->hits[i]->pointt[0];
      double y1=cl1->hits[i]->pointt[1];
      double z1=cl1->hits[i]->pointt[2]; 
 
   for( unsigned int j=0; j<cl2->hits.size();j++)
    {

      double tmp=((x1-cl2->hits[j]->pointt[0])*(x1-cl2->hits[j]->pointt[0])+
		  (y1-cl2->hits[j]->pointt[1])*(y1-cl2->hits[j]->pointt[1])+
		  (z1-cl2->hits[j]->pointt[2])*(z1-cl2->hits[j]->pointt[2]) ) ;     
      dist.push_back(tmp);
    }
    }
   sort(dist.begin(),dist.end());
   return sqrt(dist[0]);

}
inline double Dot2(double* X1,double* X2)
{
  double n1=sqrt(X1[0]*X1[0]+X1[1]*X1[1]+X1[2]*X1[2]);
  double n2=sqrt(X2[0]*X2[0]+X2[1]*X2[1]+X2[2]*X2[2]);

  return (X1[0]*X2[0]+X1[1]*X2[1]+X1[2]*X2[2])/(n1*n2);
}
void ClusterInCluster2(Tmpcl2* cl, vector<Tmpcl2*>& clv,vector<Tmpcl2*>& clout)
{

     
    sort(cl->hits.begin(),cl->hits.end());

    for(unsigned int i=0;i<clv.size();i++)
      {
	sort(clv[i]->hits.begin(),clv[i]->hits.end());
        if( includes(cl->hits.begin(),cl->hits.end(),clv[i]->hits.begin(),clv[i]->hits.end()))
	  clout.push_back(clv[i]);
      }
}
Photon2::~Photon2(){}

Photon2:: Photon2(double Ein,double* pravac, double* pocetak)
{
 Ee=Ein;
 dir[0]=pravac[0];
 dir[1]=pravac[1];
 dir[2]=pravac[2];
 start[0]=pocetak[0];
 start[1]=pocetak[1];
 start[2]=pocetak[2];

 PGdb::ZONE zons ;
 zons=PGdb::ECAL1_BAR;
  
 Z        =PGDB[zons].Zeff;
 x0eff    =PGDB[zons].x0eff;
 sampling =PGDB[zons].sampling;
 Eceff    =PGDB[zons].Eceff;
 eprime   =PGDB[zons].eprime;
 Rm       =PGDB[zons].Rmeff;

 z1=0.0251+0.00319*log(Ee*1000.0);
 z2=0.1162-0.000381*Z;
 k1=0.659-0.00309*Z;
 k2=0.645;
 k3=-2.59;
 k4=0.3585+0.0421*log(Ee*1000.0);
 p1=2.632-0.00094*Z; 
 p2=0.401+0.00187*Z;
 p3=1.313-0.0686*log(Ee*1000.0);
 y=Ee*1000.0/Eceff;

 Fs=x0eff/sampling;

 Thom=log(y)-0.858;
 Tsam=Thom-0.59/Fs-0.53*(1.0-eprime);
 alfahom=0.21+(0.492+2.38/Z)*log(y);
 alfasam=alfahom-0.444/Fs;
 betasam=0.5;

}


void Photon2::Prob(CalorimeterHit* ch,double cut,double* out)
{
 
 double X[3];
 PointOnLine22(start, dir,ch->getPosition(),X);

 double TP[3];  TP[0]=ch->getPosition()[0]-start[0];
                TP[1]=ch->getPosition()[1]-start[1];
		TP[2]=ch->getPosition()[2]-start[2];
 // distance must have direction  
 if( Dot2(TP,dir)<0.0)
   {
     out[0]=0.0;
     out[1]=0.0;
     return ;
   }

 double t=sqrt( (start[0]-X[0])*(start[0]-X[0])+
		(start[1]-X[1])*(start[1]-X[1])+
		(start[2]-X[2])*(start[2]-X[2]))/x0eff;
   double tau= t/Tsam; 
   double pos[3];  pos[0]=ch->getPosition()[0];pos[1]=ch->getPosition()[1];pos[2]=ch->getPosition()[2];

   double r= LinePointDistance2(start,dir,pos);

// average radial profiles 
   double RChom=z1+z2*tau;
   double RThom=k1*(exp(k3*(tau-k2))+exp(k4*(tau-k2)));
   double phom=p1*exp((p2-tau)/p3-exp((p2-tau)/p3));

// average radial profiles sampling
 double RCsam=RChom-0.0203*(1.0-eprime)+0.0397*exp(-tau)/Fs;
 double RTsam=RThom-0.14*(1.0-eprime)-0.495*exp(-tau)/Fs;
 double psam =phom+(1.0-eprime)*(0.348-0.642*exp(-pow(tau-1.0,2.0))/Fs);
 
 double rb=r/Rm; 
 
 if ( rb < 20.0 && t< 35.0)  
 {
   double tmp= pow(betasam*t,alfasam-1.0)*betasam*exp(-betasam*t)/gsl_sf_gamma(alfasam);
   tmp=Ee*tmp*(psam*2.0*rb*RCsam*RCsam/pow(rb*rb+RCsam*RCsam,2.0) +
	       (1-psam)*(2.0*rb*RTsam*RTsam/pow(rb*rb+RTsam*RTsam,2.0)))/rb;
     if(tmp>cut)
       {
	 out[0]= tmp;
	 out[1]= r;
       }else{
         out[0]=0.0;
	 out[1]=0.0;
       }
 }else{
   out[0]=0.0;
   out[1]=0.0;
 }
 

}



void PointOnLine3(const double* X1,const double* X2,const float* X0,double* Xline)
{

double tmp1=(X1[0]-X0[0])*(X2[0]-X1[0])+(X1[1]-X0[1])*(X2[1]-X1[1])+(X1[2]-X0[2])*(X2[2]-X1[2]);
double tmp4=(X1[0]-X2[0])*(X1[0]-X2[0])+(X1[1]-X2[1])*(X1[1]-X2[1])+(X1[2]-X2[2])*(X1[2]-X2[2]);
double t=-tmp1/tmp4;


 Xline[0]=X1[0]+(X2[0]-X1[0])*t;
 Xline[1]=X1[1]+(X2[1]-X1[1])*t;
 Xline[2]=X1[2]+(X2[2]-X1[2])*t;

}
void PointOnLine22(const double* Xstart,const double* dir,const float* X0,double* Xline)
{
  double X1[3];
  double X2[3];
  X1[0]=Xstart[0]-dir[0]*100.0;
  X1[1]=Xstart[1]-dir[1]*100.0;
  X1[2]=Xstart[2]-dir[2]*100.0;
  X2[0]=Xstart[0]+dir[0]*100.0;
  X2[1]=Xstart[1]+dir[1]*100.0;
  X2[2]=Xstart[2]+dir[2]*100.0;


  PointOnLine3(X1,X2,X0,Xline);

}

double  giveMeEEstimate2(int nivo,double Ecore, vector<CoreCalib2> cc)
{
   for(unsigned int i=0;i<cc.size();i++)
    {
      if( cc[i].level==nivo &&  cc[i].Emax>Ecore)  
	{
	  return cc[i].a+cc[i].b*Ecore;
	}
    }
    return 0.0;
}


void CreateCalibrationLDC00(vector<CoreCalib2>* cc)
{

  
  double aEnom[9]={0.5,1.0,2.0,3.0,5.0,8.0,12.0,20.0,40.0};
  double aa[9][10]={{0.283,0.357,0.337,0.398,0.000,0.000,0.000,0.000,0.000,0.000},
		    {0.703,0.346,0.814,0.785,0.939,0.000,0.000,0.000,0.000,0.000},
		    {0.505,0.821,1.000,1.112,1.527,1.860,0.000,0.000,0.000,0.000},
		    {0.809,0.892,1.999,1.554,2.170,2.870,0.000,0.000,0.000,0.000},
		    {1.075,0.827,1.290,2.350,2.123,3.958,4.814,0.000,0.000,0.000},
		    {1.140,1.310,1.765,3.652,4.030,5.600,6.670,0.000,0.000,0.000},
		    {0.760,1.845,3.444,2.850,5.150,6.760,10.72,0.000,0.000,0.000},
		    {3.770,4.460,2.270,5.930,8.730,11.00,17.40,17.77,0.000,0.000},
		    {1.890,3.560,5.120,7.320,9.650,11.71,21.74,0.000,0.000,0.000}};

  double bb[9][10]={{0.730,0.595,0.772,0.683,0.000,0.000,0.000,0.000,0.000,0.000},
		    {0.493,1.079,0.458,0.592,0.447,0.000,0.000,0.000,0.000,0.000},
                    {0.941,0.813,0.766,0.782,0.559,0.405,0.000,0.000,0.000,0.000},
                    {0.892,0.910,0.550,0.809,0.593,0.347,0.000,0.000,0.000,0.000},
                    {0.918,1.015,0.967,0.800,0.954,0.510,0.356,0.000,0.000,0.000},
		    {0.974,0.986,0.970,0.778,0.784,0.593,0.509,0.000,0.000,0.000},
	            {1.032,0.971,0.869,0.986,0.835,0.742,0.339,0.000,0.000,0.000},
	            {0.891,0.888,1.025,0.884,0.782,0.620,0.314,0.361,0.000,0.000},
	            {1.020,0.992,0.977,0.962,0.945,0.943,0.724,0.000,0.000,0.000}};

double aEmin[9][10]={{0.250,0.200,0.150,0.060,0.000,0.000,0.000,0.000,0.000,0.000},
		     {0.600,0.500,0.350,0.200,0.100,0.000,0.000,0.000,0.000,0.000},
		     {1.500,1.200,1.100,0.900,0.650,0.300,0.000,0.000,0.000,0.000},
		     {2.400,2.100,1.900,1.600,1.300,0.600,0.000,0.000,0.000,0.000},
		     {4.000,4.100,3.900,3.200,2.900,2.400,1.100,0.000,0.000,0.000},
		     {7.000,6.800,6.200,5.500,5.000,4.200,2.600,0.000,0.000,0.000},
		     {10.50,10.00,9.500,8.600,8.400,8.300,5.200,0.000,0.000,0.000},
		     {18.50,18.20,17.40,16.00,15.50,13.00,11.50,9.000,0.000,0.000},
		     {36.00,35.00,34.00,33.50,32.50,31.00,28.00,0.000,0.000,0.000}};

double aEmax[9][10]={{0.600,0.500,0.350,0.300,0.000,0.000,0.000,0.000,0.000,0.000},
		     {1.200,1.100,1.000,0.900,0.700,0.000,0.000,0.000,0.000,0.000},
		     {2.500,2.200,1.900,1.800,1.400,1.200,0.000,0.000,0.000,0.000},
		     {3.300,3.100,3.000,2.700,2.400,2.000,0.000,0.000,0.000,0.000},
		     {5.600,5.500,5.200,4.800,4.300,3.900,2.900,0.000,0.000,0.000},
		     {9.000,8.800,8.300,7.600,7.200,6.400,5.200,0.000,0.000,0.000},
		     {13.50,13.00,12.50,11.70,11.10,11.00,8.300,0.000,0.000,0.000},
		     {22.20,22.80,21.50,20.50,19.00,16.30,15.50,13.50,0.000,0.000},
		     {44.00,43.50,42.50,41.50,39.20,37.50,33.50,0.000,0.000,0.000}};

 for(unsigned int ibl=0;ibl<10;ibl++)
   {
     for(unsigned int ibk=0;ibk<9;ibk++)
       {
	 if( aEmin[ibk][ibl]!=0.0 && aEmax[ibk][ibl]!=0.0)
	   {
	     CoreCalib2 cck;
	     cck.level=ibl;
	     cck.Enom=aEnom[ibk];
	     cck.a=aa[ibk][ibl];
	     cck.b=bb[ibk][ibl];
	     cck.Emin=aEmin[ibk][ibl];
	     cck.Emax=aEmax[ibk][ibl];
	     
	     cc->push_back(cck);
	   }
       }	      
   }
 
}
