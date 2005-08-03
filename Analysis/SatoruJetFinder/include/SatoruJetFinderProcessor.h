/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/*
** This file is part of the MarlinReco Project.
** Forming part of the SubPackage: BrahmsTracking.
**
** For the latest version download from Web CVS:
** www.blah.de
**
** $Id: SatoruJetFinderProcessor.h,v 1.3 2005-08-03 13:47:54 samson Exp $
**
** $Log: not supported by cvs2svn $
** Revision 1.2  2005/08/02 18:43:53  samson
** Made the processor Marlin v00-09 compliant
**
*/ 

#ifndef SatoruJetFinderProcessor_h
#define SatoruJetFinderProcessor_h 1

#include "marlin/Processor.h"

#include "lcio.h"
#include "IMPL/LCCollectionVec.h" 


using namespace lcio ;


struct SatoruPartonArray{
    int NumberOfPartons;
    float Momentum[1200];
    int PointerParticleToJets[300];
}; 
struct SatoruJetsArray{
    int NumberOfJets;
    float Momentum[80];
};

extern "C" {
    void syjkrn_( const char* GlobalModus_,
		 // Primary jet finding
		 int &NJetRequested_,float &Threshold_,
		 int &PrimaryJetFindingMode_,float* YCut_,
		 //  Reassosiation after first jet finding
		 int &MergingMode_,float &MergingThreshold_,
		 // Second Jet finding mode
		 int &SecondJetFindingMode_,
		 // input array--> array of 
		 // PPartons(DimensionOfInputArray(4-6),NumberOfPartons)
		 int &NumberOfPartons_,
		 int &DimensionOfInputArray_,float *PPartons_,
		 // Output 
		 int &MaximalNumberOfJets_,
		 int &NJetsFound_,int *PointerParicleToJet_,
		 int &DimensionOfOutputArray_,float *PJets_,
		 float &YMinus_,float &YPlus_,
		 // error flag
		 int &IError_,int GlobalModusLenght_)
	;}


namespace marlin
{


  /**
   * \verbatim
   * A                                                                                
   * universal jetfinder module devellopped by Satoru Yamashita                       
   *                                                                                  
   * steering cards:                                                                  
   * 1) Allways to specify:                                                           
   *                                                                                  
   * InputCollection  (e.g. energy flow object )of type reconstructed particle        
   * OutputCollection (e.g. jets )of type reconstructed particle                      
   * Mode                                                                             
   *                                                                                  
   * At the moment 5 modes are implemented:                                             
   *                                                                                  
   *                                                                                  
   * DurhamNjet:                                                                      
   * durham jetfinding with a fixed number of jets to be specify                      
   * with                                                                             
   * NJetRequested                                                                    
   *                                                                                  
   * DurhamYCut:                                                                      
   * durham jetfinding with a fixed ycut to be specify with                           
   * YCut                                                                             
   *                                                                                  
   * ConeBlanka                                                                       
   * cone jet finder with cone threshold energy of 0.7 GeV and                        
   * R =0.2                                                                           
   *                                                                                  
   * Saturo                                                                           
   * first durham jetfinding, afterwards reassignment of the                          
   * objects to the jet axes found in this first iteration                            
   * with the jade sceme                                                              
   * number of jets to be specify with                                                
   * NJetRequested                                                                    
   *                                                                                  
   * Manual everything to be set by hand:                                             
   *                                                                                  
   * GlobalMode  (=MD below)                                                          
   * NJetRequested  (=NJETRQ)                                                         
   * Threshold      (= THRESH)                                                        
   * PrimaryJetFindingMode (=IMODE)                                                   
   * YCut           (=YCUTS)                                                          
   * MerginMode     (=IEXAM)                                                          
   * MergingThreshold (=FRAC)                                                         
   * SecondJetFindingMode (=IEXAM4)                                                   
   *                                                                                  
   * Description                                                                      
   * Inputs:                                                                          
   *   MD      ; mode (see below)                                                     
   *   NJETRQ  : fixed number of jets ; =0 variable -- use YCUT etc...                
   *   THRESH  : threshold (GeV) for the primary jet finding                          
   *   IMODE   : Kernel Jet finding mode                                              
   *      1- 6 : YKERN                                                                
   *       IF(IMODE.EQ.1) THEN                                                        
   *          CM = 'JADE E0'                                                          
   *        ELSEIF(IMODE.EQ.2) THEN                                                   
   *          CM = 'JADE P '                                                          
   *        ELSEIF(IMODE.EQ.3) THEN                                                   
   *          CM = 'JADE P0'                                                          
   *        ELSEIF(IMODE.EQ.4) THEN                                                   
   *          CM = 'JADE E '                                                          
   *        ELSEIF(IMODE.EQ.5) THEN                                                   
   *          CM = 'DURHAM '                                                          
   *        ELSEIF(IMODE.EQ.6) THEN                                                   
   *          CM = 'GENEVA '                                                          
   *        ELSE                                                                      
   *          WRITE(6,281) IMODE                                                      
   * 281      FORMAT(/,' ### YKERN: IMODE =',I3,' INVALID; SET TO 1 ###')             
   *          IMODE = 1                                                               
   *          CM = 'JADE E0'                                                          
   *        ENDIF                                                                     
   *                                                                                  
   *                                                                                  
   *                                                                                  
   *         7 : cambridge                                                            
   *      9-11 : LUCLUS                                                               
   *        12 : cone                                                                 
   *   YCUTS   : (1) ycut or xcut or Cone Radius (2) EPS                              
   *   IEXAM   : Particle re-association method 1-6                                   
   *            EXAM(1) = THETA                                                       
   *            EXAM(2) = PART(4)*PJSTMP(4,J)-PART(1)*PJSTMP(1,J)                     
   *     +                -PART(2)*PJSTMP(2,J)-PART(3)*PJSTMP(3,J)                    
   *            EXAM(3) = DURHAM                                                      
   *            EXAM(4) = E0JADE                                                      
   *            EXAM(5) = Geneva Scheme                                               
   *            EXAM(6) = EJADE   ! it makes too bad cross talk                       
   *                                                                                  
   *   FRAC    : core rate                                                            
   *   IEXAM4  : Jet merge method 1-8                                                 
   *            EXAM(1,N) = PL(4,I)*PL(4,J)                                           
   *            EXAM(2,N) = PL(4,I)*PL(4,J)-PL(1,I)*PL(1,J)                           
   *     +           -PL(2,I)*PL(2,J)-PL(3,I)*PL(3,J)                                 
   *            EXAM(3,N) = DURHAM(I,J)                                               
   *            EXAM(4,N) = E0JADE(I,J)                                               
   *            EXAM(5,N) = EJADE(I,J)                                                
   *            EXAM(6,N) = ANGLE(I,J)                                                
   *            EXAM(7,N) = ECOS(I,J)                                                 
   *            EXAM(8,N) = PL(6,I)*PL(6,J)                                           
   *   NPAR    : Number of "particles"                                                
   *   IDIM    : first dimension of PPAR array                                        
   *   PPAR    : Particle Momenta                                                     
   *   MXJET   : Maximum Number of Jets (buffer size for PJETS = MXJETS*JDIM)         
   *   JDIM    : first dimension of PJETS array                                       
   *                                                                                  
   * Outputs:                                                                         
   *   NJET    : Number of Jets                                                       
   *   IASS    : Jet assosiation for "particle"s                                      
   *   PJETS   : Jet momenta                                                          
   *   Y34     : YHI (sometimes not filled)                                           
   *   Y45     : YLO (sometimes not filled)                                           
   *   IERR    : Error flag  0==All O.K.                                              
   *     -9999: abnormal error                                                        
   *      -999: boundary error                                                        
   *       -99: input mismatch                                                        
   *      Others : error from various internal calls                                  
   * ==========================================================================       
   *   MD                                                                             
   * There are 11 variations                                                          
   * MD   Name      : Final NJ    : 1st process  : 2nd process : 3rd process          
   * 0a. Traditional: variable    :  fixed Ycuts :     -       :   -                  
   * 0b. Traditional: fixed Njets :  fixed Njets :     -       :   -                  
   * 0c. Just Merge : fixed Njets :  fixed Ycuts : do merge    :   -                  
   * 1a. LMODE1     : variable    :  fixed Ycuts : do SYJETF1  :   -                  
   * 1b. LMODE1     : fixed final :  fixed Ycuts : do SYJETF1  : do merge             
   * 1c. LMODE1     : fixed final :  fixed Njets : do SYJETF1  :   -                  
   * 1d. LMODE1     : fixed final :  fixed Ycuts : do merge    : do SYJETF1           
   * 2a. LMODE2     : variable    :  fixed Ycuts : do SYJETF2  :   -                  
   * 2b. LMODE2     : fixed final :  fixed Ycuts : do SYJETF2  : do merge             
   * 2c. LMODE2     : fixed final :  fixed Njets : do SYJETF2  :   -                  
   * 2d. LMODE2     : fixed final :  fixed Ycuts : do merge    : do SYJETF2           
   * ===========================================================================      
   *  MD & NJETRQ &  IMODE  & YCUTS & IEXAM &  Frac   &  IEXAM4 &     ISYMOD          
   *  0a &   -    &   1-12  &  > 0  &   -   &   -     &     -   &        1            
   *  0b &  > 0   &   1-12  &   -   &   -   &   -     &     -   &        2            
   *  0c &  > 0   &   1-12  &  > 0  &   -   &   -     &    1-8  &        3            
   *  1a &   -    &   1-12  &  > 0  &  1-6  & 0.0-1.0 &     -   &       11            
   *  1b &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       12            
   *  1c &  > 0   &   1-12  &   -   &  1-6  & 0.0-1.0 &     -   &       13            
   *  1d &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       14            
   *  2a &   -    &   1-12  &  > 0  &  1-6  & 0.0-1.0 &     -   &       21            
   *  2b &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       22            
   *  2c &  > 0   &   1-12  &   -   &  1-6  & 0.0-1.0 &     -   &       23            
   *  2d &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       24            
   * ===========================================================================      
   *                                                                                  
   *                                                                                  
   *                                                                                    
   * \endverbatim
   */
  class SatoruJetFinderProcessor : public Processor {

  
  public:
    
    virtual Processor*  newProcessor() { return new SatoruJetFinderProcessor ; }
    
  
    SatoruJetFinderProcessor() ;
    
    /** Called at the begin of the job before anything is read.
     * Use to initialize the module, e.g. book histograms.
     */
    virtual void init() ;
    
    /** Called for every run.
     */
    virtual void processRunHeader( LCRunHeader* run ) ;

    /** Called for every event - the working horse.
     */
    virtual void processEvent( LCEvent * evt ) ; 
    

    /** Called after data processing for clean up.
     */
    virtual void end() ;



  private:

    std::string _inputCollection;
    std::string _outputCollection;
    std::string _jetFindingMode;
    SatoruPartonArray _partonsWorkArray;
    SatoruJetsArray _jetsWorkArray;
    //    int _nJetVal;
    std::string _globalMode;
    int _nJetRequested;
    float _threshold;
    int _primaryJetFindingMode;
    float _yCutParam;
    float _yCut[2];
    int _mergingMode;
    float _mergingThreshold;
    int _secondJetFindingMode;
    int _debug;
    void putPartons(LCEvent * evt);
    void callSatoru(LCEvent * evt);
    void goSatoru(LCEvent * evt,LCCollection *JetsCol);                     
    void getJets(LCEvent * evt,LCCollection *JetsCol);

    /* debug routines */
    void writePartons();
    void writeJets();
    void writeParameters();
  } ;



} // namespace marlin
#endif



