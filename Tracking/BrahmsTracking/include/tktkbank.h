// A header file which defines the basics of a tpc

#ifndef Tk_Tk_Bank_h
#define Tk_Tk_Bank_h 1

#include<vector>
#include <cfortran.h>




class Tk_Tk_Bank
{

 public:
  Tk_Tk_Bank() ;
  ~Tk_Tk_Bank() ;
  void add_tk(int,int,int,int,int,int,int,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float) ; 

  void remove_tk(int) ;

  //  int tkmktkcpp(int,int,int,int,int,int,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float) ;

  int tkmktkcpp(int,int,int,int,int,int,int,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float *) ;

  void setMod_ID(int modid, int tk){tk_bank[tk].modid = modid ; } ;
  void setSubdetbits(int subdetbits, int tk){tk_bank[tk].subdetbits = subdetbits ; } ;
  void setMeasurement_code(int MesrCode, int tk){tk_bank[tk].measurement_code = MesrCode ; } ;
  void setTrackType(int tracktype, int tk){tk_bank[tk].tracktype = tracktype ; } ;
  void setNumTEs(int numtes, int tk){tk_bank[tk].numtes = numtes ; } ;
  void setCharge(int Charge, int tk){tk_bank[tk].charge = Charge ; } ;
  void setUnused(int unused, int tk){tk_bank[tk].unused = unused ; } ;
  void setNdf(int ndf, int tk){tk_bank[tk].ndf = ndf ; } ;
  void setChi2(float chi2, int tk){tk_bank[tk].chi2 = chi2 ; } ;
  void setLength(float L, int tk){tk_bank[tk].length = L ; } ;
  void setStartX(float xstart, int tk){tk_bank[tk].xstart = xstart ; } ;
  void setStartY(float ystart, int tk){tk_bank[tk].ystart = ystart ; } ;
  void setStartZ(float zstart, int tk){tk_bank[tk].zstart = zstart ; } ;
  void setEndX(float xend, int tk){tk_bank[tk].xend = xend ; } ;
  void setEndY(float yend, int tk){tk_bank[tk].yend = yend ; } ;
  void setEndZ(float zend, int tk){tk_bank[tk].zend = zend ; } ;
  void setCoord1_of_ref_point(float cord1, int tk){tk_bank[tk].coord1_of_ref_point = cord1 ; } ;
  void setCoord2_of_ref_point(float cord2, int tk){tk_bank[tk].coord2_of_ref_point = cord2 ; } ;
  void setCoord3_of_ref_point(float cord3, int tk){tk_bank[tk].coord3_of_ref_point = cord3 ; } ;
  void setTheta(float theta, int tk){tk_bank[tk].theta = theta ; } ;
  void setPhi(float phi, int tk){tk_bank[tk].phi = phi ; } ;
  void setInvp(float invp, int tk){tk_bank[tk].invp = invp ; } ;
  void setCovmatrix1(float cov1, int tk){tk_bank[tk].covmatrix1 = cov1 ; };
  void setCovmatrix2(float cov2, int tk){tk_bank[tk].covmatrix2 = cov2 ; } ;
  void setCovmatrix3(float cov3, int tk){tk_bank[tk].covmatrix3 = cov3 ; } ;
  void setCovmatrix4(float cov4, int tk){tk_bank[tk].covmatrix4 = cov4 ; } ;
  void setCovmatrix5(float cov5, int tk){tk_bank[tk].covmatrix5 = cov5 ; } ;
  void setCovmatrix6(float cov6, int tk){tk_bank[tk].covmatrix6 = cov6 ; } ;
  void setCovmatrix7(float cov7, int tk){tk_bank[tk].covmatrix7 = cov7 ; } ;
  void setCovmatrix8(float cov8, int tk){tk_bank[tk].covmatrix8 = cov8 ; } ;
  void setCovmatrix9(float cov9, int tk){tk_bank[tk].covmatrix9 = cov9 ; } ;
  void setCovmatrix10(float cov10, int tk){tk_bank[tk].covmatrix10 = cov10 ; } ;
  void setCovmatrix11(float cov11, int tk){tk_bank[tk].covmatrix11 = cov11 ; } ;
  void setCovmatrix12(float cov12, int tk){tk_bank[tk].covmatrix12 = cov12 ; } ;
  void setCovmatrix13(float cov13, int tk){tk_bank[tk].covmatrix13 = cov13 ; } ;
  void setCovmatrix14(float cov14, int tk){tk_bank[tk].covmatrix14 = cov14 ; } ;
  void setCovmatrix15(float cov15, int tk){tk_bank[tk].covmatrix15 = cov15 ; } ;

  void setPosOfFirstTEInTEList(int itkdat1, int tk){ itkdat_bank[tk].posOfFirstTEInTEList = itkdat1 ; } ;
  void setNumOfTEs(int itkdat2, int tk){ itkdat_bank[tk].numOfTEs = itkdat2 ; } ;
  void setTrackNo(int itkdat5, int tk){ itkdat_bank[tk].trackNo = itkdat5 ; } ;

  void addTE(int te, int tk){tk_bank[tk].telist.push_back(te) ; } ; 

  int   size(){return tk_bank.size() ; } ;

  int   getMod_ID(int tk){return tk_bank[tk].modid ; } ;
  int   getSubdetbits(int tk){return tk_bank[tk].subdetbits ; } ;
  int   getMeasurement_code(int tk){return tk_bank[tk].measurement_code ; } ;
  int   getTrackType(int tk){return tk_bank[tk].tracktype ; } ;
  int   getNumTEs(int tk){return tk_bank[tk].numtes ; } ;
  int   getCharge(int tk){return tk_bank[tk].charge ; } ;
  int   getUnused(int tk){return tk_bank[tk].unused ; } ;
  int   getNdf(int tk){return tk_bank[tk].ndf ; } ;
  float getChi2(int tk){return tk_bank[tk].chi2 ; } ;
  float getLength(int tk){return tk_bank[tk].length ; } ;
  float getStartX(int tk){return tk_bank[tk].xstart ; } ;
  float getStartY(int tk){return tk_bank[tk].ystart ; } ;
  float getStartZ(int tk){return tk_bank[tk].zstart ; } ;
  float getEndX(int tk){return tk_bank[tk].xend ; } ;
  float getEndY(int tk){return tk_bank[tk].yend ; } ;
  float getEndZ(int tk){return tk_bank[tk].zend ; } ;
  float getCoord1_of_ref_point(int tk){return tk_bank[tk].coord1_of_ref_point ; } ;
  float getCoord2_of_ref_point(int tk){return tk_bank[tk].coord2_of_ref_point ; } ;
  float getCoord3_of_ref_point(int tk){return tk_bank[tk].coord3_of_ref_point ; } ;
  float getTheta(int tk){return tk_bank[tk].theta ; } ;
  float getPhi(int tk){return tk_bank[tk].phi ; } ;
  float getInvp(int tk){return tk_bank[tk].invp ; } ;
  float getCovmatrix1(int tk){return tk_bank[tk].covmatrix1 ; } ;
  float getCovmatrix2(int tk){return tk_bank[tk].covmatrix2 ; } ;
  float getCovmatrix3(int tk){return tk_bank[tk].covmatrix3 ; } ;
  float getCovmatrix4(int tk){return tk_bank[tk].covmatrix4 ; } ;
  float getCovmatrix5(int tk){return tk_bank[tk].covmatrix5 ; } ;
  float getCovmatrix6(int tk){return tk_bank[tk].covmatrix6 ; } ;
  float getCovmatrix7(int tk){return tk_bank[tk].covmatrix7 ; } ;
  float getCovmatrix8(int tk){return tk_bank[tk].covmatrix8 ; } ;
  float getCovmatrix9(int tk){return tk_bank[tk].covmatrix9 ; } ;
  float getCovmatrix10(int tk){return tk_bank[tk].covmatrix10 ; } ;
  float getCovmatrix11(int tk){return tk_bank[tk].covmatrix11 ; } ;
  float getCovmatrix12(int tk){return tk_bank[tk].covmatrix12 ; } ;
  float getCovmatrix13(int tk){return tk_bank[tk].covmatrix13 ; } ;
  float getCovmatrix14(int tk){return tk_bank[tk].covmatrix14 ; } ;
  float getCovmatrix15(int tk){return tk_bank[tk].covmatrix15 ; } ;

  int getPosOfFirstTEInTEList( int tk){ return itkdat_bank[tk].posOfFirstTEInTEList ; } ;
  int getNumOfTEs( int tk){ return itkdat_bank[tk].numOfTEs ; } ;
  int getTrackNo( int tk){ return itkdat_bank[tk].trackNo ; } ;

  const std::vector<int> * getTElist(int tk){return &tk_bank[tk].telist ; } ;

 private:

struct tk_tk 
{
  int modid ;
  int subdetbits ;
  int measurement_code ;
  int tracktype ;
  int numtes ;
  int charge ;
  int unused ;
  int ndf ;
  float chi2 ;
  float length ;
  float xstart ;
  float ystart ;
  float zstart ;
  float xend ;
  float yend ; 
  float zend ;
  float coord1_of_ref_point ;
  float coord2_of_ref_point ;
  float coord3_of_ref_point ;
  float theta ;
  float phi ;
  float invp ;
  float covmatrix1 ;
  float covmatrix2 ;
  float covmatrix3 ;
  float covmatrix4 ;
  float covmatrix5 ;
  float covmatrix6 ;
  float covmatrix7 ;
  float covmatrix8 ;
  float covmatrix9 ;
  float covmatrix10 ;
  float covmatrix11 ;
  float covmatrix12 ;
  float covmatrix13 ;
  float covmatrix14 ;
  float covmatrix15 ;

  std::vector <int> telist ;

} ;

 std::vector <tk_tk> tk_bank ;

 // documentation of the ITKDAT from BRAHMS
/* * --- additional information. first index: */
/* *      1 (integer) position of first associated TE in TE list (see below) */
/* *      2 (integer) number of TEs */
/* *      5 (integer) track no. (positive if 95% of the hits belong to  */
/* *                             same track, negative if 75%, zero below.) */

struct tk_itkdat 
{
  int posOfFirstTEInTEList ; 
  int numOfTEs ;
  int trackNo ;
} ; 

 std::vector <tk_itkdat> itkdat_bank ;

} ;

// Global pointer to the Tk_Tk_Bank structure which is defined in tktkbank.cc 
extern Tk_Tk_Bank * TkTkBank ;


#endif
