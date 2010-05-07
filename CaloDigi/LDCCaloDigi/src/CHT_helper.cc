#include "CalorimeterHitType.h"
#include "CHT_helper.h"
#include <algorithm>

/** Helper functions that should go to Marlinutil/CalorimeterHitTypes.hh */

CHT::Layout layoutFromString(const std::string& name){

  std::string str( name ) ;
  std::transform( str.begin() , str.end() , str.begin(), ::tolower ) ;

  if( str.find("endcap" ) != std::string::npos )  return CHT::endcap ;
  if( str.find("barrel" ) != std::string::npos )  return CHT::barrel ;
  if( str.find("plug" )   != std::string::npos )  return CHT::plug ;
  if( str.find("ring" )   != std::string::npos )  return CHT::ring ;
  
  std::cout << " not found :" << str << " in " << name << std::endl; 
  return CHT::any ;
}

CHT::CaloID caloIDFromString(const std::string& name){

  std::string str( name ) ;
  std::transform( str.begin() , str.end() , str.begin(), ::tolower ) ;

  if( str.find("ecal" ) != std::string::npos )  return CHT::ecal ;
  if( str.find("hcal" ) != std::string::npos )  return CHT::hcal ;
  if( str.find("yoke" ) != std::string::npos )  return CHT::yoke ;
  if( str.find("lcal" ) != std::string::npos )  return CHT::lcal ;
  if( str.find("lhcal") != std::string::npos )  return CHT::lhcal ;
  if( str.find("bcal" ) != std::string::npos )  return CHT::bcal ;

  return CHT::unknown ;
}
