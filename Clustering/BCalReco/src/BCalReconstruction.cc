#include <unistd.h>
#include <vector>

#include <string>
#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TROOT.h>

#include <math.h>

#include "BCalReconstruction.h"

using namespace std;

TROOT root("root", "root");

void BCalReconstruction::Init(){

        BcCells = new BCalReconstruction::CellType** [maxlayers]; 
                for (int i=1; i<maxlayers; ++i) {
                        BcCells[i]=new BCalReconstruction::CellType* [maxrings];
                                for (int j=0; j<maxrings; ++j) {
                                        BcCells[i][j]=new BCalReconstruction::CellType [maxphis];
                                                for (int k=0; k<maxphis; ++k){
//                                                         BcCells[i][j][k]=0.;
                                BcCells[i][j][k].sRin=0.;
                                BcCells[i][j][k].sRout=0.;
                                BcCells[i][j][k].sZstart=0.;
                                BcCells[i][j][k].sZend=0.;
                                BcCells[i][j][k].sEdepNeg=0.;
                                BcCells[i][j][k].sEdepPos=0.;
                                BcCells[i][j][k].sSphi=0.;
                                BcCells[i][j][k].sSphi=0.;
                                BcCells[i][j][k].sDphi=0.;
                                BcCells[i][j][k].sPos[0]=i;
                                BcCells[i][j][k].sPos[1]=j;   
                                BcCells[i][j][k].sPos[2]=k;

                                                  }
                               }
                }


}

void BCalReconstruction::Destroy(){

Free3DArray(BcCells);
//Free3DArray(nBcCells);

}

BCalReconstruction::RecCorr BCalReconstruction::GetReconstrCoordinates(int number_layers, int number_rings, int number_pads[], BCalReconstruction::CellType ***info_detector){
 
        the_Layers = number_layers;
        the_Rings = number_rings;
          for(int j = 0; j<the_Rings;j++){
		the_Phis[j] = number_pads[j];
          }


        //read structure


                for (int i=1; i<the_Layers; ++i) {
                                for (int j=0; j<the_Rings; ++j) {
                                                for (int k=0; k<the_Phis[j]; ++k){
                                                         BcCells[i][j][k]=info_detector[i][j][k];
                                                  }
                               }
                }      


        int RingFW[maxrings][maxphis][maxlayers];
        int RingBW[maxrings][maxphis][maxlayers];
        double RenFW[maxrings][maxphis];
        double RenBW[maxrings][maxphis];

        vector<vector<int> > nRinFW = getVector(maxrings,maxphis);
        vector<vector<int> > nRinBW = getVector(maxrings,maxphis);

        //nulls in arays
        for(int r=0; r<the_Rings; r++){
          for(int p=0; p<the_Phis[r]; p++){
            RenFW[r][p]=0.;
            RenBW[r][p]=0.;
            for(int l=1; l<the_Layers; l++){  
              RingFW[r][p][l]=0;
              RingBW[r][p][l]=0;
            }
          }
        }
                  
        int index = 0;
        for(int l=1; l<the_Layers; l++){
          for(int r=0; r<the_Rings; r++){
            for(int p=0; p<the_Phis[r]; p++){
                  index++;
              if (BcCells[l][r][p].sEdepPos>0.){
////////////////                if (BcCells[l][r][p].sPos[0]>0){
                  RingFW[BcCells[l][r][p].sPos[1]][BcCells[l][r][p].sPos[2]][BcCells[l][r][p].sPos[0]] = 1;
                  RenFW[BcCells[l][r][p].sPos[1]][BcCells[l][r][p].sPos[2]] += BcCells[l][r][p].sEdepPos;  
               // cout<<" SUNT IN FW BeamCal"<<endl;
//                cout<<""<<Ren[BcCells[l][r][p].sPos[1]][BcCells[l][r][p].sPos[2]]<<endl;
//                cout<<""<<Ren[BcCells[l][r][p].sPos[1]][BcCells[l][r][p].sPos[2]]<<endl;
//                  cout<<index<< " "<< BcCells[l][r][p].sPos[1] <<" " << BcCells[l][r][p].sPos[2]<< " " <<Ren[BcCells[l][r][p].sPos[1]]$
////////////////                }
              }

              if (BcCells[l][r][p].sEdepNeg>0.){
////////////////                if (BcCells[l][r][p].sPos[0]>0){
                  RingBW[BcCells[l][r][p].sPos[1]][BcCells[l][r][p].sPos[2]][BcCells[l][r][p].sPos[0]] = 1;
                  RenBW[BcCells[l][r][p].sPos[1]][BcCells[l][r][p].sPos[2]] += BcCells[l][r][p].sEdepNeg;
               // cout<<" SUNT IN BW BeamCal"<<endl;
//                cout<<""<<Ren[BcCells[l][r][p].sPos[1]][BcCells[l][r][p].sPos[2]]<<endl;
//                  cout<<index<< " "<< BcCells[l][r][p].sPos[1] <<" " << BcCells[l][r][p].sPos[2]<< " " <<Ren[BcCells[l][r][p].sPos[1]]$
////////////////                }
              }
        
            }
          }
        }
          

      nRinFW = SearchTowers(RingFW);
      nRinBW = SearchTowers(RingBW);

        
        for(int l=1; l<the_Layers; l++){
          for(int r=0; r<the_Rings; r++){
            for(int p=0; p<the_Phis[r]; p++){

                if(nRinFW[BcCells[l][r][p].sPos[1]][BcCells[l][r][p].sPos[2]] == 1)
                  BcCells[l][r][p].sEdepPos = BcCells[l][r][p].sEdepPos;
                else 
                  BcCells[l][r][p].sEdepPos = 0;
                
                if(nRinBW[BcCells[l][r][p].sPos[1]][BcCells[l][r][p].sPos[2]] == 1)
                  BcCells[l][r][p].sEdepNeg = BcCells[l][r][p].sEdepNeg;
                else 
                  BcCells[l][r][p].sEdepNeg = 0;
          
                              
            }
          }
        }
      
     BCalReconstruction::RecCorr reco_obj = {{0},{0.}};
     BCalReconstruction::RecCorr reco_obj_fw = {{0},{0.}},reco_obj_bw = {{0},{0.}};

      reco_obj_fw = SearchClustersFW(BcCells);
     //cout << "Reconstructed side FW BeamCal = " << " " << reco_obj_fw.side[0] << " " << reco_obj_fw.side[1] << endl;

      reco_obj_bw = SearchClustersBW(BcCells);

     //cout << "Reconstructed side BW BeamCal = " << " " << reco_obj_bw.side[0] << " " << reco_obj_bw.side[1] << endl;
 
     if((reco_obj_fw.side[0] == 1 && reco_obj_fw.side[1] == 0) && (reco_obj_bw.side[0] == 0 && reco_obj_bw.side[1] == 0))
       reco_obj = reco_obj_fw;
     else if((reco_obj_bw.side[0] == 0 && reco_obj_bw.side[1] == 1) && (reco_obj_fw.side[0] == 0 && reco_obj_fw.side[1] == 0))
       reco_obj = reco_obj_bw;
     else if((reco_obj_fw.side[0] == 1 && reco_obj_fw.side[1] == 0) && (reco_obj_bw.side[0] == 0 && reco_obj_bw.side[1] == 1)){
       reco_obj.RecEne[0] = reco_obj_fw.RecEne[0];
       reco_obj.RecEne[1] = reco_obj_bw.RecEne[1];
       reco_obj.ErrEne[0] = reco_obj_fw.ErrEne[0];
       reco_obj.ErrEne[1] = reco_obj_bw.ErrEne[1];
       reco_obj.CoordX[0] = reco_obj_fw.CoordX[0];
       reco_obj.CoordX[1] = reco_obj_bw.CoordX[1];
       reco_obj.CoordY[0] = reco_obj_fw.CoordY[0];
       reco_obj.CoordY[1] = reco_obj_bw.CoordY[1];
       reco_obj.CoordZ[0] = reco_obj_fw.CoordZ[0];
       reco_obj.CoordZ[1] = reco_obj_bw.CoordZ[1];
       reco_obj.RecRad[0] = reco_obj_fw.RecRad[0];
       reco_obj.RecRad[1] = reco_obj_bw.RecRad[1];
       reco_obj.RecPhi[0] = reco_obj_fw.RecPhi[0];
       reco_obj.RecPhi[1] = reco_obj_bw.RecPhi[1];
       reco_obj.side[0] = reco_obj_fw.side[0];
       reco_obj.side[1] = reco_obj_bw.side[1];
     }  

    // cout << "Reconstructed side BeamCal = " << " " << reco_obj.side[0] << " " << reco_obj.side[1] << endl;

//     Print(reco_obj); 

//      Free2DArray(nRinFW);
//      Free2DArray(nRinBW);

      return reco_obj;
          
}

vector<vector<int> > BCalReconstruction::SearchTowers(int the_Chains[maxrings][maxphis][maxlayers]){


        //-------- count segments in stack------------
        
        //const  double  Pi    =  3.14159165;
	const  double  Pi    =  4.*atan(1.);
	   // cout << "Pi = " << Pi << endl;
        double rInner = 20.;
        double rOuter = 150.;
        double SgScl = 8.;
        int split = (int)((rOuter-rInner)/SgScl); //number of rings
        // cout << " split = "<<split<<endl;
        double dR = (rOuter-rInner)/((double)split+1);
        //cout << "dR = "<<dR<<endl;
        int nWafers = 8; //number of general segments
//        double sphi = 200.*Pi/180.; //? 
        double sdphi = 320.*Pi/180.;//?
//        double DAsphi = fmod(sphi+sdphi, 2*Pi); //mean betwee 320 and 200
        double Wafer = (320.*Pi/180.)/nWafers;//angle of general sectors
        double DAWafer = (360.*Pi/180.)-sdphi;
        double calc;
        
        //int nRs;
        int DAstart = 7;//ring after there are no daed area
        int nPhis[100], DAnPhis[100];
        double dPhis[100], DAdPhis[100];
        //creating array
        for (int f =0; f<split; f++){
          calc = dR/(rInner+dR/2+dR*f);
          dPhis[f] = Wafer/(double(int(Wafer/calc)+1));
          nPhis[f] = (int) ceil(Wafer/(dR/(rInner+dR/2+dR*f)));
          DAdPhis[f] = 0.;
          DAnPhis[f] = 0;
          if (DAstart != 0 && f>=DAstart){
            DAdPhis[f] = DAWafer/(double(int(Wafer/calc)+1));
            DAnPhis[f] =  (int) ceil(DAWafer/(dR/(rInner+dR/2+dR*f)));  
          }
          //cout<<" DAdPhis = "<<DAdPhis[f]<<" ; DAnPhis = "<<DAnPhis[f]<<endl;
        
        }

        int nline=0;
          
        int maxline=0;
        int nl=0;
        int flag=0;
        int cl=0; //count lines
          
        int Stack[100][2];

        vector<vector<int> > myVectorRin = getVector(maxrings,maxphis);

        for(int t=0; t<100; t++){
          Stack[t][0]=0;
          Stack[t][1]=0;
        }


        //------towers searching!------start-----------

        for(int r=0; r<the_Rings; r++){
          for(int p=0; p<the_Phis[r]; p++){
            nline   = 0;
            maxline = 0;
            flag    = 0;
//            Rin[r][p]=0;
            myVectorRin[r][p] = 0;
        
            for(int l=1; l<the_Layers; l++){
              if (the_Chains[r][p][l]==1){
                nline++;
              }
            }
            nl=nline;
        
        
          
            if (nline>10) {
              // cout << "r= " << r << " ; p = " << p << " ; nline = " << nline << endl;
              nl=0;
              if(the_Chains[r][p][0]==1) nl=1;
              //we are searching events after 5-th layer
              for(int l=6; l<30; l++){
                if (the_Chains[r][p][l+1]==1) {
                  nl+=1;
                  if (nl>maxline) maxline=nl;
                }
                else
                  nl=0;
              }
             maxline=maxline-1; //as first cell will be 0 always!
              //cheking if maximum is more then genetal amount of "1" cells
              if (maxline>nline) {
                //cout << "Warning!"<<endl;
                //cout<<"max -> r= " << r << " ; p = " << p << " ;  maxn = " << maxline << " ;  nline = " << nline <<endl;
                //cout<<" array : " <<endl ;
                for(int l=1; l<34; l++)  //cout <<"l -> "<<l<<" -> "<< the_Chains[r][p][l]<<endl;
                    ;
              }
              //printing all over 10 in line
              if (maxline>=10) {
                //cout<<"max -> r= " << r << " ; p = " << p << " ;  maxn = " << maxline << " ;  nline = " << nline <<endl;
                //cout<<" array : " <<endl ; 
                for(int l=1; l<34; l++)  //cout << the_Chains[r][p][l];
                //cout << endl;
                flag=1;
//                Rin[r][p]=1;
                myVectorRin[r][p] = 1;
//            cout<<"Rin = " <<Rin[r][p]<<endl;
//            cout<<"myVector = " <<myVector[r][p]<<endl;
                Stack[cl][0]=r;
                Stack[cl][1]=p;   
                cl+=1;
                
              }
              else {
//                Rin[r][p]=0;
                myVectorRin[r][p] = 0;
//                the_ChainEn[r][p]=0;
              }
            }
            //cout << " n lines -> " << nl << endl;
          }
        }
        //cerr << " n lines -> " << nl << endl;
                
                
                
        //------towers searching!-------end----------


//         return Rin;
         return myVectorRin;

}

BCalReconstruction::RecCorr BCalReconstruction::SearchClustersFW(BCalReconstruction::CellType ***info_detector){

                std::vector<BCalReconstruction::CellType> allInput;
                std::vector<cluster> clu; 
                std::vector<cluster> emptyclu;   
                std::vector<std::vector<cluster> > allClu; // vector of clusters
                
                BCalReconstruction::RecCorr reconstructed_object = {{0},{0.}};
                BCalReconstruction::CellType ***nBcCells;
                float angle;

                int IP = 3595;
                float xAngle = 0.014;

                angle=-TMath::Abs(xAngle/2);
                IP = TMath::Abs(IP);


 // std::cout << "Check clu " << clu.size() << std::endl;


	nBcCells = new BCalReconstruction::CellType** [the_Layers];
		for (int i=1; i<the_Layers; ++i) {
			nBcCells[i]=new BCalReconstruction::CellType* [the_Rings];
				for (int j=0; j<the_Rings; ++j) {
					nBcCells[i][j]=new BCalReconstruction::CellType [the_Phis[j]];
						for (int k=0; k<the_Phis[j]; ++k){
                               				nBcCells[i][j][k]=info_detector[i][j][k];
                                                  }
                               }
                }


        std::vector<std::vector<int> > Rinput;
        Rinput.resize(22);
        for (int i = 0; i < 22; ++i) Rinput.at(i).resize(160,0);
          
         for(int l=1;l<the_Layers;l++){
                for(int r=0;r<the_Rings;r++){
                        for(int p=0;p<the_Phis[r];p++){
                                ////////////// if(sPos[0]>0 <31)
                                allInput.push_back(nBcCells[l][r][p]);
                           if(nBcCells[l][r][p].sEdepPos > 0)
                              ++(Rinput.at(nBcCells[l][r][p].sPos[1]).at(nBcCells[l][r][p].sPos[2]) );
                        }
                }
        }
            
  //-------------------------------------------
        
        
  int count=0;
  for (unsigned int r = 0; r < Rinput.size(); ++r)
    for (unsigned int p = 0; p < Rinput[r].size(); ++p) {
      if(Rinput[r][p]>=1) count++;
    }
                                
    // cout<< "!!! count of towers after testline ->" << count << endl;   


        for(int l=1; l<the_Layers; l++){
          for(int r=0; r<the_Rings; r++){
            for(int p=0; p<the_Phis[r]; p++){
        
        
    if ( nBcCells[l][r][p].sEdepPos > 0. && Rinput.at(nBcCells[l][r][p].sPos[1]).at(nBcCells[l][r][p].sPos[2]) >= 1 ) {
  
      int flag = 0;
    
      for (unsigned int it = 0; it < clu.size(); ++it ) {
     
        if ( nBcCells[l][r][p].sPos[1] == clu[it].sPos[1] && nBcCells[l][r][p].sPos[2] == clu[it].sPos[2]) {
          clu[it].nclu++;
          flag = 1;           
          // sum up energies after 5-th layer
          if(nBcCells[l][r][p].sPos[0]>5 && nBcCells[l][r][p].sPos[0]<20)
            clu[it].sEdepPos += nBcCells[l][r][p].sEdepPos;
        }
      }

      if ( flag == 0 ) {
        cluster newTower;
        newTower.sPos[0] = nBcCells[l][r][p].sPos[0];
        newTower.sPos[1] = nBcCells[l][r][p].sPos[1];
        newTower.sPos[2] = nBcCells[l][r][p].sPos[2];
        newTower.sEdepPos   = nBcCells[l][r][p].sEdepPos;
        newTower.sRin    = nBcCells[l][r][p].sRin;
        newTower.sRout   = nBcCells[l][r][p].sRout;
        newTower.sZstart = nBcCells[l][r][p].sZstart;
        newTower.sZend   = nBcCells[l][r][p].sZend;
        newTower.sSphi   = nBcCells[l][r][p].sSphi;
        newTower.sDphi   = nBcCells[l][r][p].sDphi;
        newTower.in      = 1;
        clu.push_back(newTower);
        //std::cout<<" cell number = "<<clu.size()<<" R ="<<newTower.sPos[1]<<" Phi = "<<newTower.sPos[2]<<" sDphi ="<<newTower.s$
       
      }
        
    }
        
            }
          }
        }
        
  //std::cout << "found " << clu.size() << " towers in this cluster" << std::endl;
        
        
 if ( clu.size() > 1 ) {
        
    //for all peaks to test for neighbours every cell should have one neighbour as minimum
    for(unsigned int nc=0; nc< clu.size(); ++nc) {
       
      double maxEn = 0.;
      unsigned int numCell = 0;
     
      //search line with maximum energy in line
      for (unsigned int np = 0; np < clu.size(); ++np ) {
        if ( clu[np].sEdepPos > maxEn && clu[np].in==1) {
          maxEn = clu[np].sEdepPos;
          numCell = np;
        }
      } 
      clu[numCell].in=0;
      int nnei=0; //number of neighbours
        
      //search 8 neighbours
      for(unsigned int np=0; np< clu.size(); np++) {
        //compare all cells with neighbour positions
        //for the same radius
        if((clu[np].sEdepPos>0.)&&(clu[np].sPos[0]>0)&&(clu[np].sPos[0]<31)&&(clu[np].sPos[2]<1000)&&(clu[np].sPos[2]>0)) {
     
          if((clu[np].sPos[1]==clu[numCell].sPos[1] && clu[np].sPos[2]==clu[numCell].sPos[2]+1)
             ||(clu[np].sPos[1]==clu[numCell].sPos[1] && clu[np].sPos[2]==clu[numCell].sPos[2]-1)
             ||(clu[np].sPos[1]==clu[numCell].sPos[1]+1 && fabs(clu[numCell].sSphi-clu[np].sSphi)<=1.5*clu[numCell].sDphi) //sDp$
             ||(clu[np].sPos[1]==clu[numCell].sPos[1]-1 && fabs(clu[numCell].sSphi-clu[np].sSphi)<=1.5*clu[numCell].sDphi)) //sD$
            {
              ++nnei;
            }
      
        } // end if
        
      } // end loop np over clu
      
      if (nnei == 0) clu[numCell].sEdepPos = 0.;
        
    } // end loop nc over clu
     
    //print out all rest cells
    for(unsigned int t=0; t < clu.size(); ++t) {
      if((clu[t].sEdepPos>0.)&&(clu[t].sPos[0]>0)&&(clu[t].sPos[0]<31) && clu[t].sPos[2] > 0 && clu[t].sPos[2] < 1000)
        {         // std::cout<<" cell number = "<<t<<" R ="<<clu[t].sPos[1]<<" Phi = "<<clu[t].sPos[2]<<" sDphi ="<<clu[t].sDphi <<" sSphi = "<<clu[t].sSphi<< std::endl;
             
        }//if(clu[t].sEdepPos>0.)
    }//for(int t=0; t<numclu; t++)

  } // end if clu.size() > 1   

      
  //--------Searching for clusters--------------
  //start loop over all clu and search max energetic tower
  //maxEn tower mv to allClu as a beginning of new array clu[].sEdepPos = 0.; not to count this clu anymore
  //searching for neighbous -> loop overall allClu -> inside loop ove clu if any in neighbour ->
  //if yes clu[].sEdepPos = 0; and add this clu to the allClu in this loop.
      
         
//  std::vector< std::vector< float > > exampleVec;
//  std::vector<float> first;
//  exampleVec.push_back(first); 
//  for () {
//    exampleVec.back()[nnnn];
  
//  }
      
  if ( clu.size() > 1 ) {
  
    double maxEn = 0.;
    int numCell = -1;
  
    //search line with maximum energy in line
    for (unsigned int np = 0; np < clu.size(); ++np ) {
      if ( clu[np].sEdepPos > maxEn ) {
        maxEn = clu[np].sEdepPos;
        numCell = np;
      }
    }//for np
    if( numCell != -1){
      allClu.push_back(emptyclu);//just to put some onedimentional vector
      allClu.back().push_back(clu[numCell]);
      clu[numCell].sEdepPos = 0.;
    }
    
    
    //loop over allClu
    //for (std::vector<scoeff>::const_iterator firstLoop = allClu.begin(); firstLoop != allClu.end(); ++firstLoop) {
    //Loop over clusters
    int maxennei = 0;
    double EnMaxNei=0.;
    for (unsigned int fLoop =0; fLoop < allClu.size(); ++fLoop) {//loop over first coeff in allClu
       
      //Loop to find neighbours
      for (unsigned int np = 0; np < clu.size(); ++np ) {//loop over all clu to find neigb
        if (( clu[np].sEdepPos > 0)&&(clu[np].sPos[0]>0)&&(clu[np].sPos[0]<31)
            &&(clu[np].sPos[2]<1000)&&(clu[np].sPos[2]>0)) {
          if ((clu[np].sPos[1]==allClu[fLoop][0].sPos[1] && clu[np].sPos[2]==allClu[fLoop][0].sPos[2]+1)
              ||(clu[np].sPos[1]==allClu[fLoop][0].sPos[1] && clu[np].sPos[2]==allClu[fLoop][0].sPos[2]-1)
              ||(clu[np].sPos[1]==allClu[fLoop][0].sPos[1]+1
                 && fabs(allClu[fLoop][0].sSphi-clu[np].sSphi)<=1.5*allClu[fLoop][0].sDphi)
              ||(clu[np].sPos[1]==allClu[fLoop][0].sPos[1]-1
                 && fabs(allClu[fLoop][0].sSphi-clu[np].sSphi)<=1.5*allClu[fLoop][0].sDphi))
            {
              allClu[fLoop].push_back(clu[np]);//cope of element from clu - is a neighbour
              //search max energetic tower and remember it
              if (clu[np].sEdepPos > EnMaxNei){
                maxennei=allClu[fLoop].size()-1;
              }
      
              clu[np].sEdepPos = 0.;//energy set 0
            
          
            }// if looooong stuff
        } // end if
      }//loop over np
      //std::cout << " next max " << allClu[fLoop]
                 
      if(fabs( - EnMaxNei + allClu[fLoop][0].sEdepPos) > 0.9*allClu[fLoop][0].sEdepPos) {
        for (unsigned int np = 0; np < clu.size(); ++np ) {//loop over all clu to find neigb
          if (( clu[np].sEdepPos > 0)&&(clu[np].sPos[0]>0)&&(clu[np].sPos[0]<31)
              &&(clu[np].sPos[2]<1000)&&(clu[np].sPos[2]>0)) {
            if ((clu[np].sPos[1]==allClu[fLoop][maxennei].sPos[1] && clu[np].sPos[2]==allClu[fLoop][maxennei].sPos[2]+1)
                ||(clu[np].sPos[1]==allClu[fLoop][maxennei].sPos[1] && clu[np].sPos[2]==allClu[fLoop][maxennei].sPos[2]-1)
                ||(clu[np].sPos[1]==allClu[fLoop][maxennei].sPos[1]+1
                   && fabs(allClu[fLoop][maxennei].sSphi-clu[np].sSphi)<=1.5*allClu[fLoop][maxennei].sDphi)
                ||(clu[np].sPos[1]==allClu[fLoop][maxennei].sPos[1]-1
                   && fabs(allClu[fLoop][maxennei].sSphi-clu[np].sSphi)<=1.5*allClu[fLoop][maxennei].sDphi))
              {
                //      ++nnei;
                //std::cout<< "R ->" << clu[np].sPos[1] << " ; Phi ->" << clu[np].sPos[2] << std::endl;
              allClu[fLoop].push_back(clu[np]);//cope of element from clu - is a neighbour
              clu[np].sEdepPos = 0.;//energy set 0
      
        
              }// if looooong stuff
              
          } // end if
                
        }//loop over np
      }
                
                   
      maxEn = 0.;
      numCell = -1;
                
      //search line with maximum energy in line
      for (unsigned int np = 0; np < clu.size(); ++np ) {
        if ( clu[np].sEdepPos > maxEn ) {
          maxEn = clu[np].sEdepPos;
          numCell = np;
        }     
      }//for np
      if( numCell != -1){
        allClu.push_back(emptyclu);//just to put some onedimentional vector
        allClu.back().push_back(clu[numCell]);
        clu[numCell].sEdepPos = 0.;
      }            
      
      
    }//loop over firstLoop
  }//if clu.size>0
  else
    std::cerr << "no towers in calorimeter " << std::endl;
          
          
  // ------- Work with clusters ----- 2D vector -----------
  double En_cluster = 0.;
  double MaxEn_cluster = 0.;
  int clusterSize   = 0 ;
  int MaxCluSize    = 0 ;
  double CluPos        = 0. ;
  double CluPositionSearch = 0.;
  double CluPos2        = 0. ;
  double CluPositionSearch2 = 0.;
  double CluPhi        = 0. ;

  double CluPhiSearch = 0.;
  double CluPhi2        = 0. ;
  double CluPhiSearch2 = 0.;
          
  unsigned int clusterNum    = 0 ;
  //std::cout << "number of clusters" << allClu.size() << std::endl;
  //std::cout << "-------------------------------" <<  std::endl;
  
  for (unsigned int fLoop =0; fLoop < allClu.size(); ++fLoop) {//loop over first coeff in allClu
    //std::cout << "number of cluster " << fLoop << std::endl;
    En_cluster = 0.;
    clusterSize= 0;
    //std::cout << "--------Cluster Information-----------------------" <<  std::endl;
  
    for (unsigned int secLoop =0; secLoop < allClu[fLoop].size(); ++secLoop) {//loop over sec coeff
      En_cluster += allClu[fLoop][secLoop].sEdepPos;
      //std::cout << "R=" <<allClu[fLoop][secLoop].sPos[1] <<"-> Phi=" << allClu[fLoop][secLoop].sPos[2] << "-> En=" << allClu[fL$
      clusterSize += 1;
      CluPositionSearch2 += (( allClu[fLoop][secLoop].sPos[1] + 1. ) * allClu[fLoop][secLoop].sEdepPos);
      CluPositionSearch  += (  allClu[fLoop][secLoop].sPos[1] + 1. );
  
      CluPhiSearch2 += (( allClu[fLoop][secLoop].sPos[2] + 1. ) * allClu[fLoop][secLoop].sEdepPos);
      CluPhiSearch  += (  allClu[fLoop][secLoop].sPos[2] + 1. );
    }//sec loop
    
    //std::cout << "Cluster energy = " << En_cluster << "; Cluster size = " << clusterSize <<std::endl;
    //    if ( MaxCluSize < clusterSize && En_cluster < 5.) {
    if ( MaxCluSize < clusterSize ) {
      MaxCluSize = clusterSize;
      clusterNum = fLoop;
      MaxEn_cluster = En_cluster;
      CluPos2 = ( (CluPositionSearch/(En_cluster)) - 1. );
      CluPos  = ( (CluPositionSearch/clusterSize)  - 1. );
      CluPhi2 = ( (CluPhiSearch/(En_cluster)) - 1. );
      CluPhi  = ( (CluPhiSearch/clusterSize)  - 1. );
  
    }
  }//first loop
    
  //std::cout << "chosen cluster ->" << clusterNum << "-> clu En =" << MaxEn_cluster << std::endl;
  // sort vertion for output in file for later analizis

  if(
      ((int) CluPos == 0 && CluPhi >= 24) ||
      ((int) CluPos == 1 && CluPhi >= 24) ||
      ((int) CluPos == 2 && CluPhi >= 32) ||
      ((int) CluPos == 3 && CluPhi >= 40) ||
      ((int) CluPos == 4 && CluPhi >= 40) ||
      ((int) CluPos == 5 && CluPhi >= 48) ||
      ((int) CluPos == 6 && CluPhi >= 56) ||
      ((int) CluPos == 7 && CluPhi >= 72) ||
      ((int) CluPos == 8 && CluPhi >= 72) ||
      ((int) CluPos == 9 && CluPhi >= 81) ||
      ((int) CluPos == 10 && CluPhi >= 90) ||
      ((int) CluPos == 11 && CluPhi >= 90) ||
      ((int) CluPos == 12 && CluPhi >= 99) ||
      ((int) CluPos == 13 && CluPhi >= 108) ||
      ((int) CluPos == 14 && CluPhi >= 108) ||
      ((int) CluPos == 15 && CluPhi >= 117) ||
      ((int) CluPos == 16 && CluPhi >= 125))
              std::cout << "Attention!!" << " " << "Rad = "<<CluPos << "; Phi = " << CluPhi << std::endl;

  if(MaxCluSize >=2 &&
    (
      ((int) CluPos == 0 && CluPhi < 24) ||
      ((int) CluPos == 1 && CluPhi < 24) ||
      ((int) CluPos == 2 && CluPhi < 32) ||
      ((int) CluPos == 3 && CluPhi < 40) ||
      ((int) CluPos == 4 && CluPhi < 40) ||
      ((int) CluPos == 5 && CluPhi < 48) ||
      ((int) CluPos == 6 && CluPhi < 56) ||
      ((int) CluPos == 7 && CluPhi < 72) ||
      ((int) CluPos == 8 && CluPhi < 72) ||
      ((int) CluPos == 9 && CluPhi < 81) ||
      ((int) CluPos == 10 && CluPhi < 90) ||
      ((int) CluPos == 11 && CluPhi < 90) ||
      ((int) CluPos == 12 && CluPhi < 99) ||
      ((int) CluPos == 13 && CluPhi < 108) ||
      ((int) CluPos == 14 && CluPhi < 108) ||
      ((int) CluPos == 15 && CluPhi < 117) ||
      ((int) CluPos == 16 && CluPhi < 125))

    ){
  //  std::cout <<"En="<< MaxEn_cluster<<"; Rad="<<CluPos << "; Phi=" << CluPhi << std::endl;
    reconstructed_object.RecEne[0] = GetEnergyCalib(MaxEn_cluster);
    reconstructed_object.ErrEne[0] = GetEnergyErr(CluPos,CluPhi);
    reconstructed_object.RecRad[0] = (double) CluPos;
    reconstructed_object.RecPhi[0] = (double) CluPhi;
    reconstructed_object.CoordX[0] = GetCoordRotX(CluPos,CluPhi,IP,angle);
    reconstructed_object.CoordY[0] = GetCoordY(CluPos,CluPhi);
    reconstructed_object.CoordZ[0] = GetCoordRotZ(CluPos,CluPhi,IP,angle);
    reconstructed_object.side[0] = 1;
  }  
  else {
 //   std::cout << "no cluster found >=2" << std::endl;
 //   std::cout <<"En=0"<<"; Rad=0"<< "; Phi=0" << std::endl;

    reconstructed_object.RecEne[0] = 0.;
    reconstructed_object.ErrEne[0] = 0.;
    reconstructed_object.RecRad[0] = 0.;
    reconstructed_object.RecPhi[0] = 0.;
    reconstructed_object.CoordX[0] = 0.;
    reconstructed_object.CoordY[0] = 0.;
    reconstructed_object.CoordZ[0] = 0.;
    reconstructed_object.side[0] = 0;
  }
 

      allInput.clear(); 

  for (int i = 1; i < the_Layers; ++i) {
    for (int j = 0; j < the_Rings; ++j)
      delete [] nBcCells[i][j];
   
    delete [] nBcCells[i];
  }
  delete [] nBcCells;

    return reconstructed_object;

}


BCalReconstruction::RecCorr BCalReconstruction::SearchClustersBW(BCalReconstruction::CellType ***info_detector){

                std::vector<BCalReconstruction::CellType> allInput;
                std::vector<cluster> clu; 
                std::vector<cluster> emptyclu;   
                std::vector<std::vector<cluster> > allClu; // vector of clusters
                
                BCalReconstruction::RecCorr reconstructed_object = {{0},{0.}};
                BCalReconstruction::CellType ***nBcCells;
                float angle;

                int IP = 3595;
                float xAngle = 0.014;

                angle = TMath::Abs(xAngle/2);
                IP = -TMath::Abs(IP);


  //std::cout << "Check clu " << clu.size() << std::endl;



	nBcCells = new BCalReconstruction::CellType** [the_Layers];
		for (int i=1; i<the_Layers; ++i) {
			nBcCells[i]=new BCalReconstruction::CellType* [the_Rings];
				for (int j=0; j<the_Rings; ++j) {
					nBcCells[i][j]=new BCalReconstruction::CellType [the_Phis[j]];
						for (int k=0; k<the_Phis[j]; ++k){
                               				nBcCells[i][j][k]=info_detector[i][j][k];
                                                  }
                               }
                }


        std::vector<std::vector<int> > Rinput;
        Rinput.resize(22);
        for (int i = 0; i < 22; ++i) Rinput.at(i).resize(160,0);
          
         for(int l=1;l<the_Layers;l++){
                for(int r=0;r<the_Rings;r++){
                        for(int p=0;p<the_Phis[r];p++){
                                ////////////// if(sPos[0]>0 <31)
                                allInput.push_back(nBcCells[l][r][p]);
                           if(nBcCells[l][r][p].sEdepNeg > 0)
                              ++(Rinput.at(nBcCells[l][r][p].sPos[1]).at(nBcCells[l][r][p].sPos[2]) );
                        }
                }
        }
            
  //-------------------------------------------
        
        
  int count=0;
  for (unsigned int r = 0; r < Rinput.size(); ++r)
    for (unsigned int p = 0; p < Rinput[r].size(); ++p) {
      if(Rinput[r][p]>=1) count++;
    }
                                
    // cout<< "!!! count of towers after testline ->" << count << endl;   


        for(int l=1; l<the_Layers; l++){
          for(int r=0; r<the_Rings; r++){
            for(int p=0; p<the_Phis[r]; p++){
        
        
    if ( nBcCells[l][r][p].sEdepNeg > 0. && Rinput.at(nBcCells[l][r][p].sPos[1]).at(nBcCells[l][r][p].sPos[2]) >= 1 ) {
  
      int flag = 0;
    
      for (unsigned int it = 0; it < clu.size(); ++it ) {
     
        if ( nBcCells[l][r][p].sPos[1] == clu[it].sPos[1] && nBcCells[l][r][p].sPos[2] == clu[it].sPos[2]) {
          clu[it].nclu++;
          flag = 1;           
          // sum up energies after 5-th layer
          if(nBcCells[l][r][p].sPos[0]>5 && nBcCells[l][r][p].sPos[0]<20)
            clu[it].sEdepNeg += nBcCells[l][r][p].sEdepNeg;
        }
      }

      if ( flag == 0 ) {
        cluster newTower;
        newTower.sPos[0] = nBcCells[l][r][p].sPos[0];
        newTower.sPos[1] = nBcCells[l][r][p].sPos[1];
        newTower.sPos[2] = nBcCells[l][r][p].sPos[2];
        newTower.sEdepNeg   = nBcCells[l][r][p].sEdepNeg;
        newTower.sRin    = nBcCells[l][r][p].sRin;
        newTower.sRout   = nBcCells[l][r][p].sRout;
        newTower.sZstart = nBcCells[l][r][p].sZstart;
        newTower.sZend   = nBcCells[l][r][p].sZend;
        newTower.sSphi   = nBcCells[l][r][p].sSphi;
        newTower.sDphi   = nBcCells[l][r][p].sDphi;
        newTower.in      = 1;
        clu.push_back(newTower);
        //std::cout<<" cell number = "<<clu.size()<<" R ="<<newTower.sPos[1]<<" Phi = "<<newTower.sPos[2]<<" sDphi ="<<newTower.s$
       
      }
        
    }
        
            }
          }
        }
        
 // std::cout << "found " << clu.size() << " towers in this cluster" << std::endl;
        

        
 if ( clu.size() > 1 ) {
        
    //for all peaks to test for neighbours every cell should have one neighbour as minimum
    for(unsigned int nc=0; nc< clu.size(); ++nc) {
       
      double maxEn = 0.;
      unsigned int numCell = 0;
     
      //search line with maximum energy in line
      for (unsigned int np = 0; np < clu.size(); ++np ) {
        if ( clu[np].sEdepNeg > maxEn && clu[np].in==1) {
          maxEn = clu[np].sEdepNeg;
          numCell = np;
        }
      } 
      clu[numCell].in=0;
      int nnei=0; //number of neighbours
        
      //search 8 neighbours
      for(unsigned int np=0; np< clu.size(); np++) {
        //compare all cells with neighbour positions
        //for the same radius
        if((clu[np].sEdepNeg>0.)&&(clu[np].sPos[0]>0)&&(clu[np].sPos[0]<31)&&(clu[np].sPos[2]<1000)&&(clu[np].sPos[2]>0)) {
     
          if((clu[np].sPos[1]==clu[numCell].sPos[1] && clu[np].sPos[2]==clu[numCell].sPos[2]+1)
             ||(clu[np].sPos[1]==clu[numCell].sPos[1] && clu[np].sPos[2]==clu[numCell].sPos[2]-1)
             ||(clu[np].sPos[1]==clu[numCell].sPos[1]+1 && fabs(clu[numCell].sSphi-clu[np].sSphi)<=1.5*clu[numCell].sDphi) //sDp$
             ||(clu[np].sPos[1]==clu[numCell].sPos[1]-1 && fabs(clu[numCell].sSphi-clu[np].sSphi)<=1.5*clu[numCell].sDphi)) //sD$
            {
              ++nnei;
            }
      
        } // end if
        
      } // end loop np over clu
      
      if (nnei == 0) clu[numCell].sEdepNeg = 0.;
        
    } // end loop nc over clu
     
    //print out all rest cells
    for(unsigned int t=0; t < clu.size(); ++t) {
      if((clu[t].sEdepNeg>0.)&&(clu[t].sPos[0]>0)&&(clu[t].sPos[0]<31) && clu[t].sPos[2] > 0 && clu[t].sPos[2] < 1000)
        {         // std::cout<<" cell number = "<<t<<" R ="<<clu[t].sPos[1]<<" Phi = "<<clu[t].sPos[2]<<" sDphi ="<<clu[t].sDphi <<" sSphi = "<<clu[t].sSphi<< std::endl;
             
        }//if(clu[t].sEdepNeg>0.)
    }//for(int t=0; t<numclu; t++)

  } // end if clu.size() > 1   

      
  //--------Searching for clusters--------------
  //start loop over all clu and search max energetic tower
  //maxEn tower mv to allClu as a beginning of new array clu[].sEdepNeg = 0.; not to count this clu anymore
  //searching for neighbous -> loop overall allClu -> inside loop ove clu if any in neighbour ->
  //if yes clu[].sEdepNeg = 0; and add this clu to the allClu in this loop.
      
         
//  std::vector< std::vector< float > > exampleVec;
//  std::vector<float> first;
//  exampleVec.push_back(first); 
//  for () {
//    exampleVec.back()[nnnn];
  
//  }
      
  if ( clu.size() > 1 ) {
  
    double maxEn = 0.;
    int numCell = -1;
  
    //search line with maximum energy in line
    for (unsigned int np = 0; np < clu.size(); ++np ) {
      if ( clu[np].sEdepNeg > maxEn ) {
        maxEn = clu[np].sEdepNeg;
        numCell = np;
      }
    }//for np
    if( numCell != -1){
      allClu.push_back(emptyclu);//just to put some onedimentional vector
      allClu.back().push_back(clu[numCell]);
      clu[numCell].sEdepNeg = 0.;
    }
    
    
    //loop over allClu
    //for (std::vector<scoeff>::const_iterator firstLoop = allClu.begin(); firstLoop != allClu.end(); ++firstLoop) {
    //Loop over clusters
    int maxennei = 0;
    double EnMaxNei=0.;
    for (unsigned int fLoop =0; fLoop < allClu.size(); ++fLoop) {//loop over first coeff in allClu
       
      //Loop to find neighbours
      for (unsigned int np = 0; np < clu.size(); ++np ) {//loop over all clu to find neigb
        if (( clu[np].sEdepNeg > 0)&&(clu[np].sPos[0]>0)&&(clu[np].sPos[0]<31)
            &&(clu[np].sPos[2]<1000)&&(clu[np].sPos[2]>0)) {
          if ((clu[np].sPos[1]==allClu[fLoop][0].sPos[1] && clu[np].sPos[2]==allClu[fLoop][0].sPos[2]+1)
              ||(clu[np].sPos[1]==allClu[fLoop][0].sPos[1] && clu[np].sPos[2]==allClu[fLoop][0].sPos[2]-1)
              ||(clu[np].sPos[1]==allClu[fLoop][0].sPos[1]+1
                 && fabs(allClu[fLoop][0].sSphi-clu[np].sSphi)<=1.5*allClu[fLoop][0].sDphi)
              ||(clu[np].sPos[1]==allClu[fLoop][0].sPos[1]-1
                 && fabs(allClu[fLoop][0].sSphi-clu[np].sSphi)<=1.5*allClu[fLoop][0].sDphi))
            {
              allClu[fLoop].push_back(clu[np]);//cope of element from clu - is a neighbour
              //search max energetic tower and remember it
              if (clu[np].sEdepNeg > EnMaxNei){
                maxennei=allClu[fLoop].size()-1;
              }
      
              clu[np].sEdepNeg = 0.;//energy set 0
            
          
            }// if looooong stuff
        } // end if
      }//loop over np
      //std::cout << " next max " << allClu[fLoop]
                 
      if(fabs( - EnMaxNei + allClu[fLoop][0].sEdepNeg) > 0.9*allClu[fLoop][0].sEdepNeg) {
        for (unsigned int np = 0; np < clu.size(); ++np ) {//loop over all clu to find neigb
          if (( clu[np].sEdepNeg > 0)&&(clu[np].sPos[0]>0)&&(clu[np].sPos[0]<31)
              &&(clu[np].sPos[2]<1000)&&(clu[np].sPos[2]>0)) {
            if ((clu[np].sPos[1]==allClu[fLoop][maxennei].sPos[1] && clu[np].sPos[2]==allClu[fLoop][maxennei].sPos[2]+1)
                ||(clu[np].sPos[1]==allClu[fLoop][maxennei].sPos[1] && clu[np].sPos[2]==allClu[fLoop][maxennei].sPos[2]-1)
                ||(clu[np].sPos[1]==allClu[fLoop][maxennei].sPos[1]+1
                   && fabs(allClu[fLoop][maxennei].sSphi-clu[np].sSphi)<=1.5*allClu[fLoop][maxennei].sDphi)
                ||(clu[np].sPos[1]==allClu[fLoop][maxennei].sPos[1]-1
                   && fabs(allClu[fLoop][maxennei].sSphi-clu[np].sSphi)<=1.5*allClu[fLoop][maxennei].sDphi))
              {
                //      ++nnei;
                //std::cout<< "R ->" << clu[np].sPos[1] << " ; Phi ->" << clu[np].sPos[2] << std::endl;
              allClu[fLoop].push_back(clu[np]);//cope of element from clu - is a neighbour
              clu[np].sEdepNeg = 0.;//energy set 0
      
        
              }// if looooong stuff
              
          } // end if
                
        }//loop over np
      }
                
                   
      maxEn = 0.;
      numCell = -1;
                
      //search line with maximum energy in line
      for (unsigned int np = 0; np < clu.size(); ++np ) {
        if ( clu[np].sEdepNeg > maxEn ) {
          maxEn = clu[np].sEdepNeg;
          numCell = np;
        }     
      }//for np
      if( numCell != -1){
        allClu.push_back(emptyclu);//just to put some onedimentional vector
        allClu.back().push_back(clu[numCell]);
        clu[numCell].sEdepNeg = 0.;
      }            
      
      
    }//loop over firstLoop
  }//if clu.size>0
  else
    std::cerr << "no towers in calorimeter " << std::endl;
          
          
  // ------- Work with clusters ----- 2D vector -----------
  double En_cluster = 0.;
  double MaxEn_cluster = 0.;
  int clusterSize   = 0 ;
  int MaxCluSize    = 0 ;
  double CluPos        = 0. ;
  double CluPositionSearch = 0.;
  double CluPos2        = 0. ;
  double CluPositionSearch2 = 0.;
  double CluPhi        = 0. ;

  double CluPhiSearch = 0.;
  double CluPhi2        = 0. ;
  double CluPhiSearch2 = 0.;
          
  unsigned int clusterNum    = 0 ;
  //std::cout << "number of clusters" << allClu.size() << std::endl;
  //std::cout << "-------------------------------" <<  std::endl;
  
  for (unsigned int fLoop =0; fLoop < allClu.size(); ++fLoop) {//loop over first coeff in allClu
    //std::cout << "number of cluster " << fLoop << std::endl;
    En_cluster = 0.;
    clusterSize= 0;
    //std::cout << "--------Cluster Information-----------------------" <<  std::endl;
  
    for (unsigned int secLoop =0; secLoop < allClu[fLoop].size(); ++secLoop) {//loop over sec coeff
      En_cluster += allClu[fLoop][secLoop].sEdepNeg;
      //std::cout << "R=" <<allClu[fLoop][secLoop].sPos[1] <<"-> Phi=" << allClu[fLoop][secLoop].sPos[2] << "-> En=" << allClu[fL$
      clusterSize += 1;
      CluPositionSearch2 += (( allClu[fLoop][secLoop].sPos[1] + 1. ) * allClu[fLoop][secLoop].sEdepNeg);
      CluPositionSearch  += (  allClu[fLoop][secLoop].sPos[1] + 1. );
  
      CluPhiSearch2 += (( allClu[fLoop][secLoop].sPos[2] + 1. ) * allClu[fLoop][secLoop].sEdepNeg);
      CluPhiSearch  += (  allClu[fLoop][secLoop].sPos[2] + 1. );
    }//sec loop
    
    //std::cout << "Cluster energy = " << En_cluster << "; Cluster size = " << clusterSize <<std::endl;
    //    if ( MaxCluSize < clusterSize && En_cluster < 5.) {
    if ( MaxCluSize < clusterSize ) {
      MaxCluSize = clusterSize;
      clusterNum = fLoop;
      MaxEn_cluster = En_cluster;
      CluPos2 = ( (CluPositionSearch/(En_cluster)) - 1. );
      CluPos  = ( (CluPositionSearch/clusterSize)  - 1. );
      CluPhi2 = ( (CluPhiSearch/(En_cluster)) - 1. );
      CluPhi  = ( (CluPhiSearch/clusterSize)  - 1. );
  
    }
  }//first loop
    
  //std::cout << "chosen cluster ->" << clusterNum << "-> clu En =" << MaxEn_cluster << std::endl;
  // sort vertion for output in file for later analizis


  if(
      (CluPos == 0 && CluPhi >= 24) ||
      (CluPos == 1 && CluPhi >= 24) ||
      (CluPos == 2 && CluPhi >= 32) ||
      (CluPos == 3 && CluPhi >= 40) ||
      (CluPos == 4 && CluPhi >= 40) ||
      (CluPos == 5 && CluPhi >= 48) ||
      (CluPos == 6 && CluPhi >= 56) ||
      (CluPos == 7 && CluPhi >= 72) ||
      (CluPos == 8 && CluPhi >= 72) ||
      (CluPos == 9 && CluPhi >= 81) ||
      (CluPos == 10 && CluPhi >= 90) ||
      (CluPos == 11 && CluPhi >= 90) ||
      (CluPos == 12 && CluPhi >= 99) ||
      (CluPos == 13 && CluPhi >= 108) ||
      (CluPos == 14 && CluPhi >= 108) ||
      (CluPos == 15 && CluPhi >= 117) ||
      (CluPos == 16 && CluPhi >= 125))
              std::cout << "Attention!!" << " " << "Rad = "<<CluPos << "; Phi = " << CluPhi << std::endl;

  if(MaxCluSize >=2 &&
    (
      (CluPos == 0 && CluPhi < 24) ||
      (CluPos == 1 && CluPhi < 24) ||
      (CluPos == 2 && CluPhi < 32) ||
      (CluPos == 3 && CluPhi < 40) ||
      (CluPos == 4 && CluPhi < 40) ||
      (CluPos == 5 && CluPhi < 48) ||
      (CluPos == 6 && CluPhi < 56) ||
      (CluPos == 7 && CluPhi < 72) ||
      (CluPos == 8 && CluPhi < 72) ||
      (CluPos == 9 && CluPhi < 81) ||
      (CluPos == 10 && CluPhi < 90) ||
      (CluPos == 11 && CluPhi < 90) ||
      (CluPos == 12 && CluPhi < 99) ||
      (CluPos == 13 && CluPhi < 108) ||
      (CluPos == 14 && CluPhi < 108) ||
      (CluPos == 15 && CluPhi < 117) ||
      (CluPos == 16 && CluPhi < 125))
    ){
 //   std::cout <<"En="<< MaxEn_cluster<<"; Rad="<<CluPos << "; Phi=" << CluPhi << std::endl;
    reconstructed_object.RecEne[1] = GetEnergyCalib(MaxEn_cluster);
    reconstructed_object.ErrEne[1] = GetEnergyErr(CluPos,CluPhi);
    reconstructed_object.RecRad[1] = (double) CluPos;
    reconstructed_object.RecPhi[1] = (double) CluPhi;
    reconstructed_object.CoordX[1] = GetCoordRotX(CluPos,CluPhi,IP,angle);
    reconstructed_object.CoordY[1] = GetCoordY(CluPos,CluPhi);
    reconstructed_object.CoordZ[1] = GetCoordRotZ(CluPos,CluPhi,IP,angle);
    reconstructed_object.side[1] = 1;
  }  
  else {
//    std::cout << "no cluster found >=2" << std::endl;
//    std::cout <<"En=0"<<"; Rad=0"<< "; Phi=0" << std::endl;
    reconstructed_object.RecEne[1] = 0.;
    reconstructed_object.ErrEne[1] = 0.;
    reconstructed_object.RecRad[1] = 0.;
    reconstructed_object.RecPhi[1] = 0.;
    reconstructed_object.CoordX[1] = 0.;
    reconstructed_object.CoordY[1] = 0.;
    reconstructed_object.CoordZ[1] = 0.;
    reconstructed_object.side[1] = 0;
  }
 
  for (int i = 1; i < the_Layers; ++i) {
    for (int j = 0; j < the_Rings; ++j)  
      delete [] nBcCells[i][j];
                
    delete [] nBcCells[i];
  }
  delete [] nBcCells;


    return reconstructed_object;

}


double  BCalReconstruction::GetEnergyCalib(double energy){

//  return (energy*70.+4.7);  //calib BeCas
//  return (energy*71.43+10.);  //calib Marlin signal
 return (energy*78.55+0.1199);  //calib Marlin signal+bkg
//return (energy*74.91+11.45);  //calib Marlin signal+bkg 27 Feb. 2012
//return (energy*72.64 +11.69); //calib Marlin signal+bkg 1 TeV (13 May 2012)
//  return energy;
}

double  BCalReconstruction::GetCoordRotX(double ring, double pad, float IP, float angle){

  double coordX;
  double coordZ;

                //calorimeter geometry (available from GEAR file, to be read from GEAR file in the future)
                double sizepad = 7.646470588; // pad length, mm
                float inrad = 20.0; // inner radius of BeamCal, mm
                float delta = 4.0; //layer thickness, mm
                double sPhi = 3.490658504e+00; // starting Phi, rad



  float phiradial[17] = {0.232711, 0.232711, 0.174533, 0.139626, 0.139626, 0.116355, 0.0997331, 0.0872665,
                        0.0872665, 0.0775702, 0.0698132, 0.0698132, 0.0634665, 0.0581776, 0.0581776,
                        0.0537024, 0.0498666};

    if(IP > 0)
      coordZ = (double) IP + 5.*delta;
    else
      coordZ = (double) IP - 5.*delta;
    

    coordX = ((inrad+sizepad/2.) + ring*sizepad)*TMath::Cos(sPhi + (pad+1.)*phiradial[(int) ring]);

    return  (TMath::Cos(angle)*coordX - TMath::Sin(angle)*coordZ);
//    return coordX;

}

double  BCalReconstruction::GetCoordY(double ring, double pad){

                //calorimeter geometry (available in GEAR file, to be read from GEAR file in the future)
                double sizepad = 7.646470588; // mm
                float inrad = 20.0; // mm
                double sPhi = 3.490658504e+00; // rad

  float phiradial[17] = {0.232711, 0.232711, 0.174533, 0.139626, 0.139626, 0.116355, 0.0997331, 0.0872665,
                        0.0872665, 0.0775702, 0.0698132, 0.0698132, 0.0634665, 0.0581776, 0.0581776,
                        0.0537024, 0.0498666};


    
    return ((inrad+sizepad/2.) + ring*sizepad)*TMath::Sin(sPhi + (pad+1.)*phiradial[(int) ring]);

}

double  BCalReconstruction::GetCoordRotZ(double ring, double pad, float IP, float angle){

  double coordX;
  double coordZ;

                //calorimeter geometry (available in GEAR file, to be read from GEAR file in the future)
                double sizepad = 7.646470588; // mm
                float inrad = 20.0; // mm
                float delta = 4.0; //mm
                double sPhi = 3.490658504e+00; // rad


  float phiradial[17] = {0.232711, 0.232711, 0.174533, 0.139626, 0.139626, 0.116355, 0.0997331, 0.0872665,
                        0.0872665, 0.0775702, 0.0698132, 0.0698132, 0.0634665, 0.0581776, 0.0581776,
                        0.0537024, 0.0498666};
      
    if(IP > 0)
      coordZ = (double) IP + 5.*delta;
    else
      coordZ = (double) IP - 5.*delta;
      

    coordX = ((inrad+sizepad/2.) + ring*sizepad)*TMath::Cos(sPhi + (pad+1.)*phiradial[(int) ring]);
    
    return  (TMath::Sin(angle)*coordX + TMath::Cos(angle)*coordZ);
//   return coordZ;

}

double  BCalReconstruction::GetEnergyErr(double ring, double pad){


//      cout << "SUNT in GetEnergyErr" << endl;



   TFile ftmp("bg_aver_LDC_3.5T_14mrad_AntiDID_NominalBeamParam.root");

  if ( ftmp.IsZombie() ) {
      cerr << "Could not read data file. Exit." << endl;
      exit(1);
   }

//      cout << "TEST GetEnergyErr 1" << endl;

//      ftmp.ls();

   TTree *ttmp = (TTree*) ftmp.Get("tBcDensAverage");

//      cout << "TEST GetEnergyErr 2" << endl;


  if ( ttmp->IsZombie() ) {
      cerr << "Could not find data tree \"tBcDensAverage\". Exit." << endl;
      exit(1);
   }



//      cout << "TEST GetEnergyErr 3" << endl;

  Int_t Nentr = ttmp->GetEntries();
  ftmp.Close();


//      cout << "Nentr = " <<Nentr << endl;


//  scoeff bgc[Nentr];
  std::vector<scoeff>  bgc(Nentr) ;
  for (int i=0; i<Nentr; i++){
    bgc[i].sPos[0] = 0;
    bgc[i].sPos[1] = 0;
    bgc[i].sPos[2] = 0;
    bgc[i].sRin    = 0.;
    bgc[i].sRout   = 0.;
    bgc[i].sZstart = 0.;
    bgc[i].sZend   = 0.;
    bgc[i].sSphi   = 0.;
    bgc[i].sDphi   = 0.;
    bgc[i].sEdep   = 0.;
    bgc[i].sCellArea = 0.;
    bgc[i].sEnDens = 0.;
    bgc[i].sEnDensErr = 0.;
  }

  int sPos[3];
  double sRin, sRout;
  double sZstart, sZend;
  double sSphi, sDphi;
  double sEdep;
  double sCellArea,sEnDens,sEnDensErr;

    
   TFile f("bg_aver_LDC_3.5T_14mrad_AntiDID_NominalBeamParam.root");
  if ( f.IsZombie() ) {
      cerr << "Could not read data file. Exit." << endl;
      exit(1);
   }

//      f.ls();

   TTree *t = (TTree*)f.Get("tBcDensAverage");


  if(t->IsZombie()) {
      cerr<<"bgcoeff: Tree \"tBcDensAverage\" not found"<<endl;
    }


//      cout << "TEST GetEnergyErr 4" << endl;


    t->SetBranchAddress("sPos",    sPos);
    t->SetBranchAddress("sRin",   &sRin);  // cm
    t->SetBranchAddress("sRout",  &sRout); // cm
    t->SetBranchAddress("sZstart",&sZstart);
    t->SetBranchAddress("sZend",  &sZend);
    t->SetBranchAddress("sSphi",  &sSphi); // rad
    t->SetBranchAddress("sDphi",  &sDphi); // rad
    t->SetBranchAddress("sEdep",  &sEdep); // GeV
    t->SetBranchAddress("sCellArea", &sCellArea);
    t->SetBranchAddress("sEnDens", &sEnDens);
    t->SetBranchAddress("sEnDensErr", &sEnDensErr);


//      cout << "TEST GetEnergyErr 5" << endl;


    for(int i=0; i<Nentr; i++){
      t->GetEntry(i);
      bgc[i].sPos[0] = sPos[0];
      bgc[i].sPos[1] = sPos[1];
      bgc[i].sPos[2] = sPos[2];
      bgc[i].sRin = sRin;
      bgc[i].sRout = sRout;
      bgc[i].sZstart = sZstart;
      bgc[i].sZend = sZend;
      bgc[i].sSphi = sSphi;
      bgc[i].sDphi = sDphi;
      bgc[i].sEdep = sEdep;
      bgc[i].cell_area = sCellArea;  // mm2
      bgc[i].sEnDens = sEnDens; // GeV/mm2
      bgc[i].sEnDensErr = sEnDensErr;
    }

     double sqav_EnDensErr = 0.;

    for(int i=0; i<Nentr; i++)
        if(bgc[i].sPos[0] >= 0 && bgc[i].sPos[0] <= 30 && bgc[i].sPos[1] == ring && bgc[i].sPos[2] == pad){
             sqav_EnDensErr += bgc[i].sEnDensErr * bgc[i].sEnDensErr ;
//             cout << "Check rms " << " " << bgc[i].sEnDensErr << " " << bgc[i].sEnDensErr * bgc[i].sEnDensErr << " " << sqav_EnDensErr << endl;    
      }


    f.Close();

   return sqrt(sqav_EnDensErr);

}


void BCalReconstruction::Print(BCalReconstruction::RecCorr obj){

     if(obj.side[0] == 1)
        cout << "FW BeamCal " <<obj.RecEne[0]  << " * " << obj.ErrEne[0] <<" * " << obj.CoordX[0]  << " * " <<obj.CoordY[0] <<" * " << obj.CoordZ[0] <<endl;
     else if(obj.side[1] == 1)
        cout << "BW BeamCal " <<obj.RecEne[1]  << " * " << obj.ErrEne[1] <<" * " << obj.CoordX[1]  << " * " <<obj.CoordY[1] <<" * " << obj.CoordZ[1] <<endl;
}

vector<vector<int> > BCalReconstruction::getVector(int rows, int cols){

   vector<vector<int> > retValue(rows);
   vector<vector<int> >::iterator i;
   for(i = retValue.begin(); i != retValue.end(); ++i) {
	i->resize(cols);
   }
   
   return retValue;
}

void BCalReconstruction::Free2DArray(int **p2DArray){

  for (int i = 0; i < maxrings; ++i)
    delete [] p2DArray[i];
  delete [] p2DArray;

}

void BCalReconstruction::Free3DArray(BCalReconstruction::CellType ***p3DArray){ 

  for (int i = 2; i < maxlayers; ++i) {
    for (int j = 0; j < maxrings; ++j)
      delete [] p3DArray[i][j];

    delete [] p3DArray[i];
  }
  delete [] p3DArray;

}



