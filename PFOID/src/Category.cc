#include "Category.h"
#include <iostream>

Category::Category(std::string Name, int NoOfHist){
  catName = Name ;
  nHists   = NoOfHist;
  if(nHists < 1){
    std::cout << " A Category needs at least one histogram! 2 would be better! " << std::endl;
    exit(-1);
  }
  
};


Category::~Category(){
  for(unsigned int i=0; i<hists.size(); i++) delete hists[i];
  hists.clear();
};


std::string Category::GetName(){
  return catName;
};


unsigned int Category::GetNoOfHists(){
  return nHists;
};
