#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <math.h>
#include "Histogram.h"


int cindex(int *index, int* maxima, int N){
  int ind = index[N-1];
  for(int i=N-2; i>-1; i--){
    ind = ind*maxima[i]+index[i];
  }
  return ind;
}


Histogram::Histogram(std::string HistName, int Dimension){
  dim      = Dimension ;
  if(dim<1){
    std::cout << " Please let you help! How do you imagine a 0-dimensional histodram ?!" << std::endl;
    exit(-1);
  }
  histName = HistName ;
  varNames = new std::string[dim] ;
  start    = new double[dim] ;
  width    = new double[dim] ;
  bins     = new int[dim] ;

  for (int i=0; i<dim; i++){
    std::string dummy = "var" ;
    dummy += i ;
    varNames[i] = dummy ;
    start[i] = 0.;
    width[i] = 1.;
    bins[i]  = 1 ;
  }
  totNoOfBins = 1;

  content = new double[totNoOfBins] ;
  content[0] = 1.e-20;
  norm = content[0];
  vol = 1. ;
  for (int i=0; i<3; i++) density[i]=content[0];
}

Histogram::Histogram(std::string HistName, int Dimension, std::string *VarName, double *StartVal, double *BinWidth, int *NoOfBins){
  dim      = Dimension ;
  if(dim<1){
    std::cout << " Please let you help! How do you imagine a histodram of dimension < 1 ?!" << std::endl;
    exit(-1);
  }
  histName = HistName ;
  varNames = new std::string[dim] ;
  start    = new double[dim] ;
  width    = new double[dim] ;
  bins     = new int[dim] ;

  totNoOfBins = 1;
  vol = 1. ;
  for (int i=0; i<dim; i++){
    varNames[i] = VarName[i] ;
    start[i]    = StartVal[i] ;
    width[i]    = BinWidth[i] ;
    bins[i]     = NoOfBins[i] ;
    totNoOfBins *= NoOfBins[i] ;
    vol *= BinWidth[i] ;
  }

  content = new double[totNoOfBins] ;
  content[0] = 1.e-20;
  norm = totNoOfBins*content[0];
  for (int i=0; i<3; i++) density[i]=content[0];
}


Histogram::~Histogram(){
  delete [] varNames;
  delete [] start;
  delete [] width;
  delete [] bins;
  delete [] content;
}


int Histogram::SetVariable(int VarNo, std::string VarName, double StartVal, double BinWidth, int NoOfBins){
  if(VarNo<0 || VarNo>=dim){
    std::cout << " Index out of range " << std::endl;
    return -1;
  }
  for (int i=0; i<dim; i++)
    if( i!=VarNo && varNames[i]==VarName ){
      std::cout << " Name of variable already exists!" << std::endl;
      return 1;
    }

  varNames[VarNo] = VarName ;
  start[VarNo]    = StartVal ;
  width[VarNo]    = BinWidth ;
  bins[VarNo]     = NoOfBins ;
  
  totNoOfBins = 1;
  for(int i=0; i<dim; i++) totNoOfBins *= bins[i];
  
  delete [] content;
  content = new double[totNoOfBins];
  for (int i=0; i<totNoOfBins; i++) content[i] = 1.e-20;
  norm = 1.e-20;
  vol=1.;
  for (int i=0; i<dim; i++) vol *= width[i];
  for (int i=0; i<3; i++) density[i]=content[0];
  //  std::cout << "histogram new initialized - empty!" << std::endl;

  return 0;
}


int Histogram::SetVariables(std::string *VarName, double *StartVal, double *BinWidth, int *NoOfBins){
  for(int i=0; i<dim; i++)
    for(int j=0; j<dim; j++)
      if( i!=j && VarName[i]==VarName[j] ){
	std::cout << " At least one of the variable names are double!" << std::endl;
	return -1;
      }
  for (int i=0; i<dim; i++){
    varNames[i] = VarName[i] ;
    start[i]    = StartVal[i] ;
    width[i]    = BinWidth[i] ;
    bins[i]     = NoOfBins[i] ;
  }

  totNoOfBins = 1;
  for(int i=0; i<dim; i++){
    totNoOfBins *= bins[i];
    vol *= width[i];
  }

  delete [] content;
  content = new double[totNoOfBins];
  for (int i=0; i<totNoOfBins; i++) content[i] = 1.e-20;
  norm = 1.e-20;
  for (int i=0; i<3; i++) density[i]=content[0];
  //  std::cout << "histogram new initialized - empty!" << std::endl;
  
  return 0;
}


int Histogram::Fill(VObject *VO) {
  std::string *VOnames = VO->GetNames();
  bool found ;
  for(int i=0; i<dim; i++){
    found=false;
    for(int j=0; j<VO->GetNoOfVariables(); j++)
      if(VOnames[j]==varNames[i]) found = true;
    if(!found){
      std::cerr << " Histogram variable " << varNames[i] << " not contained in VObject!" << std::endl;
      return -1;
    }
  }
  
  std::cerr << "Hist DBG: start looking for bin " << dim << std::endl;
  double val;

  //FIXED:SJA  int index[dim]; //change the variable size array to dynamic memory 
  int *index ;
  index = new int[dim];
  


  for(int i=0; i<dim; i++){
    val=VO->GetValue(varNames[i]);

    std::string str="Hist DBG: val[";
    str+=varNames[i];
    str+="]";
    std::cerr << std::setw(30) << std::left << str << std::setw(3) << " = " << std::setw(10) << std::left << val << std::endl;

    if ( val<start[i] ){                   // 0-th bin:  |.x.|...|... at underflow
      index[i] = 0;
      break;
    }
    if ( val>start[i]+bins[i]*width[i] ){  // last bin:   ...|...|.x.| at overflow
      index[i] = bins[i]-1 ;
      break;
    }
    int j=0;
    while(j<=bins[i]){                      // while j< #bins
      j++;
      if( val<=start[i]+j*width[i] ){       // if val < start + #bins*width -> index = j-1
	index[i] = j-1 ;
	break;
      }
    }
  }
  
  //  std::cerr << "Hist DBG: index=" << index << std::endl;
  
  int mind = cindex(index,bins,dim);
  //  std::cerr << "Hist DBG: mind=" << mind << "  totNoOfBins=" << totNoOfBins << std::endl;
  content[mind]+=1.;
  norm+=1.;
  if(content[mind]/vol>density[0]) density[0] = content[mind]/vol;
  double mini=1.e100;
  for (int i=0; i<totNoOfBins; i++)
    if(content[i]<mini) mini=content[i];
  density[1]=mini;
  density[2]=norm/(vol*totNoOfBins);

  delete[] index; 
  return 0;  
}

int Histogram::SetBinContent(int BinInd, double Value){
  if(BinInd<0 || BinInd>=totNoOfBins){
    std::cout << " Index out of range!" << std::endl;
    return -1;
  }
  if(Value<0){
    std::cout << " Bin content smaller then 0 !" << std::endl;
    return -1;
  }
  content[BinInd]=Value;

  norm+=Value;
  if(content[BinInd]/vol>density[0]) density[0] = content[BinInd]/vol;
  double mini=1.e100;
  for (int i=0; i<totNoOfBins; i++)
    if(content[i]<mini) mini=content[i];
  density[1]=mini;
  density[2]=norm/(vol*totNoOfBins);

  return 0;
}


double Histogram::GetNorm(){
  return norm;
}

double Histogram::GetVolume(){
  return vol;
}

std::string Histogram::GetHistName(){
  return histName;
}


double* Histogram::GetDensities(){
  return density;
}


double Histogram::GetContent(VObject *VO){
  std::string *VOnames = VO->GetNames();
  bool found ;
  for(int i=0; i<dim; i++){
    found=false;
    for(int j=0; j<VO->GetNoOfVariables(); j++)
      if(VOnames[j]==varNames[i]) found = true;
    if(!found){
      std::cout << " Histogram variable " << varNames[i] << " not contained in VObject!" << std::endl;
      return -1;
    }
  }
  
  double val;

  //FIXED:SJA  int index[dim]; //change the variable size array to dynamic memory 
  int *index ;
  index = new int[dim];

  for(int i=0; i<dim; i++){
    val=VO->GetValue(varNames[i]);
    if ( val<start[i] ){                   // 0-th bin:  |.x.|...|... at underflow
      index[i] = 0;
      break;
    }
    if ( val>start[i]+bins[i]*width[i] ){  // last bin:   ...|...|.x.| at overflow
      index[i] = bins[i]-1 ;
      break;
    }
    int j=0;
    while(j<=bins[i]){                      // while j< #bins
      j++;
      if( val<=start[i]+j*width[i] ){       // if val < start + #bins*width -> index = j-1
	index[i] = j-1 ;
	break;
      }
    }
  }
  
  int mind = cindex(index,bins,dim);

  delete[] index;

  if (mind>totNoOfBins-1)
    mind = totNoOfBins-1;

  if ( mind<0 ) mind = 0;

  double result = content[mind];

  if (result <=0 ) result = 1e-20;
  
  return result;

}


double Histogram::GetNormContent(VObject *VO){
  std::string *VOnames = VO->GetNames();
  bool found ;
  for(int i=0; i<dim; i++){
    found=false;
    for(int j=0; j<VO->GetNoOfVariables(); j++)
      if(VOnames[j]==varNames[i]) found = true;
    if(!found){
      std::cout << " Histogram variable " << varNames[i] << " not contained in VObject!" << std::endl;
      return -1;
    }
  }
  
  double val;

  //FIXED:SJA  int index[dim]; //change the variable size array to dynamic memory 
  int *index ;
  index = new int[dim];

  for(int i=0; i<dim; i++){
    val=VO->GetValue(varNames[i]);
    if ( val<start[i] ){                   // 0-th bin:  |.x.|...|... at underflow
      index[i] = 0;
      break;
    }
    if ( val>start[i]+bins[i]*width[i] ){  // last bin:   ...|...|.x.| at overflow
      index[i] = bins[i]-1 ;
      break;
    }
    int j=0;
    while(j<=bins[i]){                      // while j< #bins
      j++;
      if( val<=start[i]+j*width[i] ){       // if val < start + #bins*width -> index = j-1
	index[i] = j-1 ;
	break;
      }
    }
  }
  
  int mind = cindex(index,bins,dim);

  delete[] index;

  if (norm < 1e-20) norm = 1e-20;

  if (mind>totNoOfBins-1)
    mind = totNoOfBins-1;

  double result = content[mind]/norm;
  if (result <= 0 ) result = 1e-20;
  return result;

}


double Histogram::GetBinContent(int BinInd){
  if (BinInd<0 || BinInd>=totNoOfBins){
    std::cout << " Index is out of Range!" << std::endl;
    return 0;
  }
  return content[BinInd];
} 

int Histogram::GetDimension(){
  return dim;
} 
 
int Histogram::GetTotNoOfBins(){
  return totNoOfBins;
}

std::string * Histogram::GetVarNames(){
  return varNames;
}


double * Histogram::GetVarStart(){
  return start;
}


double * Histogram::GetVarWidth(){
  return width;
}


int * Histogram::GetVarBins(){
  return bins;
}


double Histogram::CompareHist(Histogram *hist){
  if( dim!=hist->GetDimension() ){
    std::cout << " Histogramms aren't compatible - wrong dimensions" << std::endl;
    return -1;
  }
  for(int i=0; i<dim; i++){
    if(varNames[i]!=hist->GetVarNames()[i]){
      std::cout << " Histogramms aren't compatible - some variables differ (at least in order)" << std::endl;
      return -1;
    }
    if(fabs(start[i]-hist->GetVarStart()[i])>1.e-10){
      std::cout << " Histogramms aren't compatible - start values differ (at least in order)" << std::endl;
      return -1;
    }
    if(fabs(width[i]-hist->GetVarWidth()[i])>1.e-10){
      std::cout << " Histogramms aren't compatible - widthes differ (at least in order) " << width[i] << "  " << hist->GetVarWidth()[i] << std::endl;
      return -1;
    }
    if(bins[i]-hist->GetVarBins()[i]!=0){
      std::cout << " Histogramms aren't compatible - numbers of bins differ (at least in order)" << std::endl;
      return -1;
    }
  }

  double totDiff=0.;
  for(int i=0; i<totNoOfBins; i++){
    double diff = content[i] - hist->GetBinContent(i);
    totDiff += diff*diff;
  }

  return sqrt(totDiff)/(vol*totNoOfBins);
}

