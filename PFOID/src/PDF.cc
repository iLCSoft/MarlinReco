#include "PDF.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <math.h>

PDF::PDF(std::string PDFname, int NoOfCats, std::string* CatNames, int NoOfHists, int NoOfVar, std::string *VarNames){
  pdfName = PDFname;
  nCats  = NoOfCats;
  IniCount = 0;
  if(nCats < 1){
    std::cout << " PDF needs at least one category! 2 make more sense! " << std::endl;
    exit(-1);
  }
  for(int i=0; i<nCats; i++){
    for(int j=0; j<nCats; j++){
      if(i!=j && CatNames[i]==CatNames[j]){
	std::cout << " There are 2 categories with the same name! Not acceptable!!" << std::endl;
	exit(1);
      }
    }
  }
  for(int i=0; i<nCats; i++){
    Category *ct = new Category(CatNames[i], NoOfHists);
    cat.push_back(ct);
  }
  
  VO = new VObject(NoOfVar);
  if ( VO->SetNames(VarNames) !=0 ){
    std::cout << " ERROR in VObject initialization" << std::endl;
    exit(-1);
  };

  std::cout << " PDF created with categories and non-initialized histograms" << std::endl;
}


PDF::PDF(std::string Filename){
  std::ifstream in(Filename.c_str());
  if(!in){
    std::cout << " File " << Filename << " does not exist!" << std::endl;
    exit(1);
  }

  // reading
  in >> pdfName;                                             //  name of the pdf
  if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  in >> nCats;                                               //  number of categories in pdf
  if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  
  //FIXED:SJA std::string CatNames[nCats];                               //  names of the categories
  //replaced variable length array
  std::string *CatNames = new std::string[nCats];

  for (int i=0; i<nCats; i++){
    in >> CatNames[i];
    if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  }

  int NoOfVar;
  in >> NoOfVar;                                             // number of variables
  if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }

  //FIXED:SJA std::string VarNames[NoOfVar];                             // array with names of variables
  //replaced variable length array
  std::string *VarNames = new std::string[NoOfVar];

  for(int i=0; i<NoOfVar; i++){
    in >> VarNames[i];
    if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  }
  int NoOfHists;
  in >> NoOfHists;                                           //  number of histograms in each category
  if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }

  // create categories
  for(int i=0; i<nCats; i++){
    for(int j=0; j<nCats; j++){
      if(i!=j && CatNames[i]==CatNames[j]){
	std::cout << " There are 2 categories with the same name! Not acceptable!!" << std::endl;
	exit(1);
      }
    }
  }
  for(int i=0; i<nCats; i++){
    Category *ct = new Category(CatNames[i], NoOfHists);
    cat.push_back(ct);
  }

  VO = new VObject(NoOfVar);
  if ( VO->SetNames(VarNames) !=0 ){
    std::cout << " ERROR in VObject initialization" << std::endl;
    exit(-1);
  };

  std::cout << " PDF created with categories and non-initialized histograms" << std::endl;

  // -----------------

  //FIXED:SJA int dims[NoOfHists];                                       // dimensions of the histograms
  //replaced variable length array
  int *dims = new int[NoOfHists];

  for (int i=0; i<NoOfHists; i++){
    in >> dims[i];
    if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  }

  //FIXED:SJA std::string HistNames[NoOfHists];                          // Names of the histograms
  //replaced variable length array
  std::string *HistNames = new std::string[NoOfHists];
   
  for (int i=0; i<NoOfHists; i++){
    in >> HistNames[i];
    if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  }

  for (int i=0; i<NoOfHists; i++){                           // Read histogram structure

    //FIXED:SJA std::string VarName[dims[i]];
    //replaced variable length array
    std::string *VarName = new std::string[dims[i]];

    //FIXED:SJA double start[dims[i]];
    //replaced variable length array
    double *start = new double[dims[i]];
    
    //FIXED:SJA double binWidth[dims[i]];
    //replaced variable length array
    double *binWidth = new double [dims[i]];
    
    //FIXED:SJA int bins[dims[i]];
    //replaced variable length array
    int *bins = new int[dims[i]];

    for (int j=0; j<dims[i]; j++){
      in >> VarName[j];
      if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
      in >> start[j] >> binWidth[j] >> bins[j];
      if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
    }
    for(int j=0; j<nCats; j++){
      Histogram *h = new Histogram(HistNames[i],dims[i],VarName,start,binWidth,bins);
      cat[j]->hists.push_back(h);
    }
    
    delete [] VarName;
    delete [] start;
    delete [] binWidth;
    delete [] bins;

  }

  double value;
  for(int j=0; j<nCats; j++){
    for (int i=0; i<NoOfHists; i++){                           // read content
      for(int k=0; k<cat[j]->hists[i]->GetTotNoOfBins(); k++){
	in >> value;
	if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
	cat[j]->hists[i]->SetBinContent(k, value);
      }
    }
  }
  std::cout << " PDF complete" << std::endl;

  in.close();

  delete [] VarNames;
  delete [] dims;
  delete [] HistNames;

}


PDF::~PDF(){
  for(unsigned int i=0; i<cat.size(); i++) delete cat[i];
  cat.clear();
  delete VO;
}
 
int PDF::InitHistogram(std::string HistName, int Dim, std::string *VarName, double *StartVal, double *BinWidth, int *NoOfBins){
  int ret=0;

  // for all categories
  for(int i=0; i<nCats; i++){
    if(cat[i]->GetNoOfHists()==cat[i]->hists.size()){
      std::cout << " Histograms are already initialized!" << std::endl;
      return -1;
    }

    for(unsigned int j=0; j<cat[i]->hists.size(); j++)
      if(cat[i]->hists[j]->GetHistName()==HistName){
	std::cout << " There is already a histogram with name '" << HistName << "'!\n";
	return -1;
      }
    
    if(Dim>VO->GetNoOfVariables()){
      std::cout << " Dimension of histogram can be at most the # of variables in VObject " << VO->GetNoOfVariables() << std::endl;
      return -1;
    }
    
    bool found = false;
    for(int k=0; k<Dim; k++){
      found=false;
      for(int l=0; l<VO->GetNoOfVariables(); l++){
	if(VO->GetName(l)==VarName[k]) found=true;
      }
      if(!found){
	std::cout << " At least one histogram variable is not contained in VObject!\n";
	return -1;
      }
      Histogram *h = new Histogram(HistName,Dim);
      h->SetVariables(VarName, StartVal, BinWidth, NoOfBins);
      cat[i]->hists.push_back(h);
    }
    if(cat[i]->GetNoOfHists()==cat[i]->hists.size()){
      std::cout << " All histograms of category " << this->GetCatName(i) << " initialized" << std::endl;
      IniCount++;
    }
  }
  
  if(IniCount==cat.size()){
    std::cout << " PDF " << pdfName << " initialized" << std::endl;
  }
  
  return ret;
}


int PDF::FillHistograms(std::string CatName){
  int ret = -2;

  for(unsigned int i=0; i<cat.size(); i++){
    if(cat[i]->GetName()==CatName){
      if(cat[i]->GetNoOfHists()!=cat[i]->hists.size()){
	std::cout << " In category " << CatName << " are still un-initialized histograms! "; 
	std::cout << cat[i]->GetNoOfHists() << "   " << cat[i]->hists.size() << std::endl;
	return -2;
      }
      for(unsigned int j=0; j<cat[i]->hists.size(); j++){
	cat[i]->hists[j]->Fill(VO) ;
      }
      return 0;
    }
  }
  if (ret==-2) std::cout << " Error when Filling the histogram: No matching histogram found!" << std::endl;
  return ret;
}


Histogram * PDF::GetHistogram(std::string CatName, std::string HistName){
  for(unsigned int i=0; i<cat.size(); i++){
    if(cat[i]->GetName()==CatName){
      for(unsigned int j=0; j<cat[i]->hists.size(); j++){
	if(cat[i]->hists[j]->GetHistName()==HistName){
	  return cat[i]->hists[j];
	}
      }
    }
  }
  return NULL;
}


double PDF::GetLikelihood(std::string CatName){
  
  //FIXED:SJA:removed variable array:  double p [nCats][cat[0]->GetNoOfHists()];
  std::vector < std::vector < double > > p(nCats,std::vector < double >(cat[0]->GetNoOfHists()));  

  //FIXED:SJA:removed variable array:  double sum [cat[0]->GetNoOfHists()];
  double *sum = new double[cat[0]->GetNoOfHists()];

   for(unsigned int j_his=0; j_his<cat[0]->GetNoOfHists(); j_his++){
     sum[j_his]=0.;
     for(int i_cat=0; i_cat<nCats; i_cat++){
       p[i_cat][j_his] = cat[i_cat]->hists[j_his]->GetNormContent(VO); // normalized f_ij
       sum[j_his]+=p[i_cat][j_his];                                    // sum over f's
     }
   }

   //FIXED:SJA:removed variable array:  double prod [nCats];
   double *prod = new double[nCats];

   for(int i_cat=0; i_cat<nCats; i_cat++){
     prod[i_cat] = 1.;
     for(unsigned int j_his=0; j_his<cat[0]->GetNoOfHists(); j_his++){
       p[i_cat][j_his] /= sum[j_his];         // final pdfs i,j : i = category index, j = histogram index
       prod[i_cat]*=p[i_cat][j_his];          // product over p's
     }
   }

   
   //FIXED:SJA:removed variable array:  double LH [nCats];
   std::vector<double> LH(nCats, 0.0);

   double sum2=0.;

   int index=-2;
   for(int i_cat=0; i_cat<nCats; i_cat++){
     LH[i_cat] = prod[i_cat];
     if(CatName==cat[i_cat]->GetName()){
       index=i_cat;
     }
     sum2+=prod[i_cat];
   }

   if(index<0 || index>nCats-1){
     std::cout << " Error in Likelihood " << std::endl;
     delete[] sum;
     delete[] prod;
     return -1;
   }

    double sum3=0.;
    for(int i_cat=0; i_cat<nCats; i_cat++){
      sum3+=LH[i_cat];
    }  
    sum3/=sum2;
    if(fabs(sum3-1.)>1.e-10){
      std::cout << " Error in Likelihood : total LH not 1!" << std::endl;
      delete[] sum;
      delete[] prod;
      return -1;
    }

    delete[] sum;
    delete[] prod;
    return LH[index]/sum2;
}


// ------  PDF::ReadPDF  ------------
int PDF::ReadPDF(std::string Filename){
  std::ifstream in(Filename.c_str());

  // clean everything
  for(unsigned int i=0; i<cat.size(); i++) delete cat[i];
  cat.clear();
  delete VO;

  // reading
  in >> pdfName;                                             //  name of the pdf
  if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  in >> nCats;                                               //  number of categories in pdf
  if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  
  //FIXED:SJA:removed variable array:    std::string CatNames[nCats];                               //  names of the categories
  std::string *CatNames = new std::string[nCats];

  for (int i=0; i<nCats; i++){
    in >> CatNames[i];
    if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  }

  int NoOfVar;
  in >> NoOfVar;                                             // number of variables
  if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }

  //FIXED:SJA:removed variable array:  std::string VarNames[NoOfVar];                             // array with names of variables
  std::string *VarNames = new std::string[NoOfVar];                             // array with names of variables

  for(int i=0; i<NoOfVar; i++){
    in >> VarNames[i];
    if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  }
  int NoOfHists;
  in >> NoOfHists;                                           //  number of histograms in each category
  if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }

  // create categories
  for(int i=0; i<nCats; i++){
    for(int j=0; j<nCats; j++){
      if(i!=j && CatNames[i]==CatNames[j]){
	std::cout << " There are 2 categories with the same name! Not acceptable!!" << std::endl;
	exit(1);
      }
    }
  }
  for(int i=0; i<nCats; i++){
    Category *ct = new Category(CatNames[i], NoOfHists);
    cat.push_back(ct);
  }

  VO = new VObject(NoOfVar);
  if ( VO->SetNames(VarNames) !=0 ){
    std::cout << " ERROR in VObject initialization" << std::endl;
    exit(-1);
  };

  std::cout << " PDF created with categories and non-initialized histograms" << std::endl;

  // -----------------

  //FIXED:SJA:removed variable array:    int dims[NoOfHists];                                       // dimensions of the histograms
  int *dims = new int[NoOfHists];
  

  for (int i=0; i<NoOfHists; i++){
    in >> dims[i];
    if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  }

  //FIXED:SJA:removed variable array:    std::string HistNames[NoOfHists];                          // Names of the histograms
  std::string *HistNames = new std::string[NoOfHists];                          // Names of the histograms


  for (int i=0; i<NoOfHists; i++){
    in >> HistNames[i];
    if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
  }

  for (int i=0; i<NoOfHists; i++){                           // Read histogram structure

    
    //FIXED:SJA:removed variable array:    std::string VarName[dims[i]];
    std::string *VarName = new std::string[dims[i]];

    //FIXED:SJA:removed variable array:        double start[dims[i]], binWidth[dims[i]];
    double *start = new double[dims[i]];
    double *binWidth = new double[dims[i]];
    int *bins = new int[dims[i]];

    for (int j=0; j<dims[i]; j++){
      in >> VarName[j];
      if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
      in >> start[j] >> binWidth[j] >> bins[j];
      if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
    }
    for(int j=0; j<nCats; j++){
      Histogram *h = new Histogram(HistNames[i],dims[i],VarName,start,binWidth,bins);
      cat[j]->hists.push_back(h);
    }
    delete[] VarName;
    delete[] start;
    delete[] binWidth;
    delete[] bins;
  }

  double value;
  for(int j=0; j<nCats; j++){
    for (int i=0; i<NoOfHists; i++){                           // read content
      for(int k=0; k<cat[j]->hists[i]->GetTotNoOfBins(); k++){
	in >> value;
	if(!in.good()){ std::cout << " Error in reading PDF file" << std::endl;  exit(1); }
	cat[j]->hists[i]->SetBinContent(k, value);
      }
    }
  }

  in.close();
  std::cout << " PDF complete" << std::endl;

  delete[] CatNames;
  delete[] VarNames;
  delete[] dims;
  delete[] HistNames;
  return 0;
}


// ------  PDF::WritePDF  ------------
int PDF::WritePDF(std::string Filename){
  std::ofstream of(Filename.c_str());

  of << pdfName << std::endl;                                //  pdf name
  of << nCats << std::endl;                                  //  number of categories
  for (int i=0; i<nCats; i++)
    of << cat[i]->GetName() << std::endl;                    //  names of the categories

  of << VO->GetNoOfVariables() << std::endl;
  for (int i=0; i<VO->GetNoOfVariables(); i++)
    of << VO->GetName(i) << std::endl;

  of << cat[0]->GetNoOfHists() << std::endl;                 //  number of histograms
  for (unsigned int i=0; i<cat[0]->GetNoOfHists(); i++)
    of << cat[0]->hists[i]->GetDimension() << std::endl;     //  dimension of the histogram i
  for (unsigned int i=0; i<cat[0]->GetNoOfHists(); i++)
    of << cat[0]->hists[i]->GetHistName() << std::endl;      //  name of the histogram i
  
  for (unsigned int i=0; i<cat[0]->GetNoOfHists(); i++)
    for (int j=0; j<cat[0]->hists[i]->GetDimension(); j++){
      of << cat[0]->hists[i]->GetVarNames()[j] << std::endl; //  names of the variable j in hist i
      of << cat[0]->hists[i]->GetVarStart()[j] << "   " ;    //  start values of variable j
      of << cat[0]->hists[i]->GetVarWidth()[j] << "   " ;    //  bin width of variable j
      of << cat[0]->hists[i]->GetVarBins()[j] << std::endl ; //  Number of bins of variable j
    }
  
  //  double check;
  for (int k=0; k<nCats; k++){                                 //  sum over categories
    for (unsigned int i=0; i<cat[k]->GetNoOfHists(); i++){              //  sum over histograms
      //      check=0;
      for (int j=0; j<cat[k]->hists[i]->GetTotNoOfBins(); j++){//  sum over all bins
	of << cat[k]->hists[i]->GetBinContent(j) << "  ";
	//	check+=cat[k]->hists[i]->GetBinContent(j);
      }
      of << std::endl;

      //      std::cout << " Norm " << cat[k]->hists[i]->GetNorm() << " " << check << std::endl;

    }
  }

  of.close();
  
  return 0;
}


int PDF::GetNCats(){ 
  return nCats; 
}

std::string PDF::GetCatName(int CatNo){
  if(CatNo<0 || CatNo>=nCats){
    std::cout << " Index out of range! " << std::endl;
    return "";
  }
  return cat[CatNo]->GetName();
}
