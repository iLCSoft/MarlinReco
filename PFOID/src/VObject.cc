#include "VObject.h"

#include <iostream>


VObject::VObject(int NoOfVars){
  vn = NoOfVars;
  varNames = new std::string[vn];
  varValues = new double[vn];
  for(int i=0; i<vn; i++){
    varNames[i] = "";
    varValues[i] = 0.;
  }
}

VObject::~VObject(){
  delete [] varNames;
  delete [] varValues;
}

int VObject::SetName(int VarNo, std::string VarName){
  if(VarNo<0 || VarNo>=vn){
    std::cout << " Variable index out of range!" << std::endl;
    return -1;
  }
  for(int i=0; i<vn; i++)
    if( i!=VarNo && varNames[i]==VarName ){
      std::cout << " A variable with name '" << VarName << "' exists already!" << std::endl;
      return -1;
    }

  varNames[VarNo] = VarName;
  return 0;
}

int VObject::SetNames(std::string *Names){
  for(int i=0; i<vn; i++)
    for(int j=0; j<vn; j++)
      if( i!=j && Names[i]==Names[j] ){
	std::cout << " Variable name occurs twice!" << std::endl;
	return -1;
      }
  for (int i=0; i<vn; i++) varNames[i] = Names[i];

  return 0;
}

int VObject::SetValue(int VarNo, double Value){
  if(VarNo<0 || VarNo>=vn) return -1;
  varValues[VarNo] = Value;
  return 0;
}

int VObject::SetValue(std::string VarName, double Value){
  for(int i=0; i<vn; i++)
    if(varNames[i]==VarName){
      varValues[i]=Value;
      return 0;
    }

  std::cout << " No variable with that name '" << VarName << "' found!" << std::endl;
  return -1;
}

double VObject::GetValue(int VarNo) const{
  if(VarNo<0 || VarNo>=vn){
    std::cout << " Variable index is out of range!" << std::endl;
    return 0;
  }
  return varValues[VarNo];
}

double VObject::GetValue(std::string VarName) const{
  for(int i=0; i<vn; i++)
    if(varNames[i]==VarName) return varValues[i];

  std::cout << " No variable with that name '" << VarName << "' found!" << std::endl;
  return 0;
}

int VObject::GetNoOfVariables() const{
  return vn;
}

std::string VObject::GetName(int VarNo) const{
  if(VarNo<0 || VarNo>=vn){
    std::cout << " Variable index out of range!" << std::endl;
    return "";
  }
  return varNames[VarNo];
}

std::string* VObject::GetNames() const{
  return varNames;
}
