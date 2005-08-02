// A header file which defines the basics of a tpc

#ifndef Tk_Te_Bank_h
#define Tk_Te_Bank_h 1

#include<vector>
#include <cfortran.h>




class Tk_Te_Bank
{

 public:
  Tk_Te_Bank() ;
  ~Tk_Te_Bank() ;
  void add_te(int,int,int,int,int,int,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float) ; 

  void remove_te(int) ;

  int tkmktecpp(int,int,int,int,int,int,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float) ;


  void   setSubdetector_ID(int subid, int te){te_bank[te].subdetector_ID = subid ; } ;
  void   setSubmodule(int submod, int te){te_bank[te].submodule = submod ; } ;
  void   setUnused(int unused, int te){te_bank[te].unused = unused ; } ;
  void   setMeasurement_code(int MesrCode, int te){te_bank[te].measurement_code = MesrCode ; } ;
  void   setPointer_to_end_of_TE(int PnteTE, int te){te_bank[te].pointer_to_end_of_TE = PnteTE ; } ;
  void   setCharge(int Charge, int te){te_bank[te].charge = Charge ; } ;
  void   setNdf(int ndf, int te){te_bank[te].ndf = ndf ; } ;
  void setChi2(float chi2, int te){te_bank[te].chi2 = chi2 ; } ;
  void setLength(float L, int te){te_bank[te].length = L ; } ;
  void setCoord1_of_ref_point(float cord1, int te){te_bank[te].coord1_of_ref_point = cord1 ; } ;
  void setCoord2_of_ref_point(float cord2, int te){te_bank[te].coord2_of_ref_point = cord2 ; } ;
  void setCoord3_of_ref_point(float cord3, int te){te_bank[te].coord3_of_ref_point = cord3 ; } ;
  void setTheta(float theta, int te){te_bank[te].theta = theta ; } ;
  void setPhi(float phi, int te){te_bank[te].phi = phi ; } ;
  void setInvp(float invp, int te){te_bank[te].invp = invp ; } ;
  void setDe_dx(float dedx, int te){te_bank[te].de_dx = dedx ; } ;
  void setCovmatrix1(float cov1, int te){te_bank[te].covmatrix1 = cov1 ; };
  void setCovmatrix2(float cov2, int te){te_bank[te].covmatrix2 = cov2 ; } ;
  void setCovmatrix3(float cov3, int te){te_bank[te].covmatrix3 = cov3 ; } ;
  void setCovmatrix4(float cov4, int te){te_bank[te].covmatrix4 = cov4 ; } ;
  void setCovmatrix5(float cov5, int te){te_bank[te].covmatrix5 = cov5 ; } ;
  void setCovmatrix6(float cov6, int te){te_bank[te].covmatrix6 = cov6 ; } ;
  void setCovmatrix7(float cov7, int te){te_bank[te].covmatrix7 = cov7 ; } ;
  void setCovmatrix8(float cov8, int te){te_bank[te].covmatrix8 = cov8 ; } ;
  void setCovmatrix9(float cov9, int te){te_bank[te].covmatrix9 = cov9 ; } ;
  void setCovmatrix10(float cov10, int te){te_bank[te].covmatrix10 = cov10 ; } ;
  void setCovmatrix11(float cov11, int te){te_bank[te].covmatrix11 = cov11 ; } ;
  void setCovmatrix12(float cov12, int te){te_bank[te].covmatrix12 = cov12 ; } ;
  void setCovmatrix13(float cov13, int te){te_bank[te].covmatrix13 = cov13 ; } ;
  void setCovmatrix14(float cov14, int te){te_bank[te].covmatrix14 = cov14 ; } ;
  void setCovmatrix15(float cov15, int te){te_bank[te].covmatrix15 = cov15 ; } ;

  void setPosOfFirstHitInHitList(int itedat1, int te){ itedat_bank[te].posOfFirstHitInHitList = itedat1 ; } ;
  void setNumOfHits(int itedat2, int te){ itedat_bank[te].numOfHits = itedat2 ; } ;
  void setPointrToFirstExclusion(int itedat3, int te){ itedat_bank[te].pointrToFirstExclusion = itedat3 ; } ;
  void setNumOfExclusions(int itedat4, int te){ itedat_bank[te].numOfExclusions = itedat4 ; } ;
  void setTrackNo(int itedat5, int te){ itedat_bank[te].trackNo = itedat5 ; } ;

  void addHit(int hit, int te){te_bank[te].hitlist.push_back(hit) ; } ; 

  int   size(){return te_bank.size() ; } ;

  int   getSubdetector_ID(int te){return te_bank[te].subdetector_ID ; } ;
  int   getSubmodule(int te){return te_bank[te].submodule ; } ;
  int   getUnused(int te){return te_bank[te].unused ; } ;
  int   getMeasurement_code(int te){return te_bank[te].measurement_code ; } ;
  int   getPointer_to_end_of_TE(int te){return te_bank[te].pointer_to_end_of_TE ; } ;
  int   getCharge(int te){return te_bank[te].charge ; } ;
  int   getNdf(int te){return te_bank[te].ndf ; } ;
  float getChi2(int te){return te_bank[te].chi2 ; } ;
  float getLength(int te){return te_bank[te].length ; } ;
  float getCoord1_of_ref_point(int te){return te_bank[te].coord1_of_ref_point ; } ;
  float getCoord2_of_ref_point(int te){return te_bank[te].coord2_of_ref_point ; } ;
  float getCoord3_of_ref_point(int te){return te_bank[te].coord3_of_ref_point ; } ;
  float getTheta(int te){return te_bank[te].theta ; } ;
  float getPhi(int te){return te_bank[te].phi ; } ;
  float getInvp(int te){return te_bank[te].invp ; } ;
  float getDe_dx(int te){return te_bank[te].de_dx ; } ;
  float getCovmatrix1(int te){return te_bank[te].covmatrix1 ; } ;
  float getCovmatrix2(int te){return te_bank[te].covmatrix2 ; } ;
  float getCovmatrix3(int te){return te_bank[te].covmatrix3 ; } ;
  float getCovmatrix4(int te){return te_bank[te].covmatrix4 ; } ;
  float getCovmatrix5(int te){return te_bank[te].covmatrix5 ; } ;
  float getCovmatrix6(int te){return te_bank[te].covmatrix6 ; } ;
  float getCovmatrix7(int te){return te_bank[te].covmatrix7 ; } ;
  float getCovmatrix8(int te){return te_bank[te].covmatrix8 ; } ;
  float getCovmatrix9(int te){return te_bank[te].covmatrix9 ; } ;
  float getCovmatrix10(int te){return te_bank[te].covmatrix10 ; } ;
  float getCovmatrix11(int te){return te_bank[te].covmatrix11 ; } ;
  float getCovmatrix12(int te){return te_bank[te].covmatrix12 ; } ;
  float getCovmatrix13(int te){return te_bank[te].covmatrix13 ; } ;
  float getCovmatrix14(int te){return te_bank[te].covmatrix14 ; } ;
  float getCovmatrix15(int te){return te_bank[te].covmatrix15 ; } ;

  int getPosOfFirstHitInHitList( int te ){ return itedat_bank[te].posOfFirstHitInHitList ; } ;
  int getNumOfHits( int te ){ return itedat_bank[te].numOfHits ; } ;
  int getPointrToFirstExclusion( int te ){ return itedat_bank[te].pointrToFirstExclusion ; } ;
  int getNumOfExclusions( int te ){ return itedat_bank[te].numOfExclusions ; } ;
  int gettrackNo( int te ){ return itedat_bank[te].trackNo ; } ;

  const std::vector<int> * getHitlist(int te){return &te_bank[te].hitlist ; } ;

 private:

struct tk_te 
{
  int subdetector_ID ;
  int submodule ;
  int unused ;
  int measurement_code ;
  int pointer_to_end_of_TE ;
  int charge ;
  int ndf ;
  float chi2 ;
  float length ;
  float coord1_of_ref_point ;
  float coord2_of_ref_point ;
  float coord3_of_ref_point ;
  float theta ;
  float phi ;
  float invp ;
  float de_dx ;
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

  std::vector <int> hitlist ;

} ;

 std::vector <tk_te> te_bank ;

 // documentation of the ITEDAT from BRAHMS
/* * --- additional information. first index: */
/* *      1 (integer) position of first associated hit in hit list (see below) */
/* *      2 (integer) number of hits */
/* *      3 (integer) pointer to first exclusion (0 if none) in exclusion list */
/* *      4 (integer) number of exclusions */
/* *      5 (integer) track no. (positive if 95% of the hits belong to  */
/* *                             same track, negative if 75%, zero below.) */
//*     remark: the exclusion list is shared with the hit bank (see there).

struct tk_itedat 
{
  int posOfFirstHitInHitList ; 
  int numOfHits ;
  int pointrToFirstExclusion ;
  int numOfExclusions ;
  int trackNo ;
} ; 

 std::vector <tk_itedat> itedat_bank ;

} ;

// Global pointer to the Tk_Te_Bank structure which is defined in tktebank.cc 
extern Tk_Te_Bank * TkTeBank ;


#endif
