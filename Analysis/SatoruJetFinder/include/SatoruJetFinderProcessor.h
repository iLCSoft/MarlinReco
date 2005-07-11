#ifndef SatoruJetFinderProcessor_h
#define SatoruJetFinderProcessor_h 1

#include "marlin/Processor.h"

#include "lcio.h"
#include "IMPL/LCCollectionVec.h" 


using namespace lcio ;
using namespace marlin;

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

/**
 * A\n                                                            
 * universal jetfinder module devellopped by Satoru Yamashita\n
 *\n
 * steering cards:\n
 * 1) Allways to specify:\n
 *
 * InputCollection  (e.g. energy flow object )of type reconstructed particle\n
 * OutputCollection (e.g. jets )of type reconstructed particle\n
 * Mode\n
 * \n
 * At the moment 5 modes are implemented:                                     \n  
 *									      \n
 *									      \n
 * DurhamNjet:								      \n
 * durham jetfinding with a fixed number of jets to be specify		      \n
 * with									      \n
 * NJetRequested							      \n
 *									      \n
 * DurhamYCut:								      \n
 * durham jetfinding with a fixed ycut to be specify with		      \n
 * YCut									      \n
 *									      \n
 * ConeBlanka								      \n
 * cone jet finder with cone threshold energy of 0.7 GeV and		      \n
 * R =0.2								      \n
 * 									      \n
 * Saturo								      \n
 * first durham jetfinding, afterwards reassignment of the		      \n
 * objects to the jet axes found in this first iteration		      \n
 * with the jade sceme							      \n
 * number of jets to be specify with 					      \n
 * NJetRequested							      \n
 *									      \n
 * Manual everything to be set by hand:					      \n
 *									      \n
 * GlobalMode  (=MD below)						      \n
 * NJetRequested  (=NJETRQ)						      \n
 * Threshold      (= THRESH)						      \n
 * PrimaryJetFindingMode (=IMODE)					      \n
 * YCut           (=YCUTS)						      \n
 * MerginMode     (=IEXAM)						      \n
 * MergingThreshold (=FRAC)						      \n
 * SecondJetFindingMode (=IEXAM4)					      \n
 * 									      \n
 * Description 								      \n
 * Inputs:								      \n
 *   MD      ; mode (see below)						      \n
 *   NJETRQ  : fixed number of jets ; =0 variable -- use YCUT etc...	      \n
 *   THRESH  : threshold (GeV) for the primary jet finding		      \n
 *   IMODE   : Kernel Jet finding mode 					      \n
 *      1- 6 : YKERN							      \n
 *       IF(IMODE.EQ.1) THEN						      \n
 *          CM = 'JADE E0'						      \n
 *        ELSEIF(IMODE.EQ.2) THEN					      \n
 *          CM = 'JADE P '						      \n
 *        ELSEIF(IMODE.EQ.3) THEN					      \n
 *          CM = 'JADE P0'						      \n
 *        ELSEIF(IMODE.EQ.4) THEN					      \n
 *          CM = 'JADE E '						      \n
 *        ELSEIF(IMODE.EQ.5) THEN					      \n
 *          CM = 'DURHAM '						      \n
 *        ELSEIF(IMODE.EQ.6) THEN					      \n
 *          CM = 'GENEVA '						      \n
 *        ELSE								      \n
 *          WRITE(6,281) IMODE						      \n
 * 281      FORMAT(/,' ### YKERN: IMODE =',I3,' INVALID; SET TO 1 ###')	      \n
 *          IMODE = 1							      \n
 *          CM = 'JADE E0'						      \n
 *        ENDIF            						      \n
 *									      \n
 *									      \n
 * 									      \n
 *         7 : cambridge						      \n
 *      9-11 : LUCLUS							      \n
 *        12 : cone							      \n
 *   YCUTS   : (1) ycut or xcut or Cone Radius (2) EPS			      \n
 *   IEXAM   : Particle re-association method 1-6			      \n
 *            EXAM(1) = THETA						      \n
 *            EXAM(2) = PART(4)*PJSTMP(4,J)-PART(1)*PJSTMP(1,J)		      \n
 *     +                -PART(2)*PJSTMP(2,J)-PART(3)*PJSTMP(3,J)	      \n
 *            EXAM(3) = DURHAM						      \n
 *            EXAM(4) = E0JADE						      \n
 *            EXAM(5) = Geneva Scheme					      \n
 *            EXAM(6) = EJADE   ! it makes too bad cross talk		      \n
 *      								      \n
 *   FRAC    : core rate						      \n
 *   IEXAM4  : Jet merge method 1-8					      \n
 *            EXAM(1,N) = PL(4,I)*PL(4,J)				      \n
 *            EXAM(2,N) = PL(4,I)*PL(4,J)-PL(1,I)*PL(1,J)		      \n
 *     +           -PL(2,I)*PL(2,J)-PL(3,I)*PL(3,J)			      \n
 *            EXAM(3,N) = DURHAM(I,J)					      \n
 *            EXAM(4,N) = E0JADE(I,J)					      \n
 *            EXAM(5,N) = EJADE(I,J)					      \n
 *            EXAM(6,N) = ANGLE(I,J)					      \n
 *            EXAM(7,N) = ECOS(I,J)					      \n
 *            EXAM(8,N) = PL(6,I)*PL(6,J)				      \n
 *   NPAR    : Number of "particles"					      \n
 *   IDIM    : first dimension of PPAR array				      \n
 *   PPAR    : Particle Momenta						      \n
 *   MXJET   : Maximum Number of Jets (buffer size for PJETS = MXJETS*JDIM)   \n
 *   JDIM    : first dimension of PJETS array				      \n
 *									      \n
 * Outputs:								      \n
 *   NJET    : Number of Jets						      \n
 *   IASS    : Jet assosiation for "particle"s				      \n
 *   PJETS   : Jet momenta						      \n
 *   Y34     : YHI (sometimes not filled)				      \n
 *   Y45     : YLO (sometimes not filled)				      \n
 *   IERR    : Error flag  0==All O.K.					      \n
 *     -9999: abnormal error						      \n
 *      -999: boundary error						      \n
 *       -99: input mismatch						      \n
 *      Others : error from various internal calls			      \n
 * ========================================================================== \n
 *   MD  								      \n
 * There are 11 variations						      \n
 * MD   Name      : Final NJ    : 1st process  : 2nd process : 3rd process    \n
 * 0a. Traditional: variable    :  fixed Ycuts :     -       :   -	      \n
 * 0b. Traditional: fixed Njets :  fixed Njets :     -       :   -	      \n
 * 0c. Just Merge : fixed Njets :  fixed Ycuts : do merge    :   -	      \n
 * 1a. LMODE1     : variable    :  fixed Ycuts : do SYJETF1  :   -	      \n
 * 1b. LMODE1     : fixed final :  fixed Ycuts : do SYJETF1  : do merge	      \n
 * 1c. LMODE1     : fixed final :  fixed Njets : do SYJETF1  :   -	      \n
 * 1d. LMODE1     : fixed final :  fixed Ycuts : do merge    : do SYJETF1     \n
 * 2a. LMODE2     : variable    :  fixed Ycuts : do SYJETF2  :   -	      \n
 * 2b. LMODE2     : fixed final :  fixed Ycuts : do SYJETF2  : do merge	      \n
 * 2c. LMODE2     : fixed final :  fixed Njets : do SYJETF2  :   -	      \n
 * 2d. LMODE2     : fixed final :  fixed Ycuts : do merge    : do SYJETF2     \n
 * ===========================================================================\n
 *  MD & NJETRQ &  IMODE  & YCUTS & IEXAM &  Frac   &  IEXAM4 &     ISYMOD    \n
 *  0a &   -    &   1-12  &  > 0  &   -   &   -     &     -   &        1      \n
 *  0b &  > 0   &   1-12  &   -   &   -   &   -     &     -   &        2      \n
 *  0c &  > 0   &   1-12  &  > 0  &   -   &   -     &    1-8  &        3      \n
 *  1a &   -    &   1-12  &  > 0  &  1-6  & 0.0-1.0 &     -   &       11      \n
 *  1b &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       12      \n
 *  1c &  > 0   &   1-12  &   -   &  1-6  & 0.0-1.0 &     -   &       13      \n
 *  1d &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       14      \n
 *  2a &   -    &   1-12  &  > 0  &  1-6  & 0.0-1.0 &     -   &       21      \n
 *  2b &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       22      \n
 *  2c &  > 0   &   1-12  &   -   &  1-6  & 0.0-1.0 &     -   &       23      \n
 *  2d &  > 0   &   1-12  &  > 0  &  1-6  & 0.0-1.0 &    1-8  &       24      \n
 * ===========================================================================\n
 *									      \n
 *									      \n
 *                                                                            \n  
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


protected:
  LCWriter* _lcWrt ;
  int _nRun ;
  int _nEvt ;

private:

  std::string InputCollection;
  std::string OutputCollection;
  std::string JetFindingMode;
  SatoruPartonArray PartonsWorkArray;
  SatoruJetsArray JetsWorkArray;
  int NJetVal;
  std::string GlobalMode;
  int NJetRequested;
  float Threshold;
  int PrimaryJetFindingMode;
  float YCut[2];
  int MergingMode;
  float MergingThreshold;
  int SecondJetFindingMode;
  int Debug;
  void PutPartons(LCEvent * evt);
  void CallSatoru(LCEvent * evt);
  void GoSatoru(LCEvent * evt,LCCollection *JetsCol);                     
  void GetJets(LCEvent * evt,LCCollection *JetsCol);

  /* debug routines */
  void WritePartons();
  void WriteJets();
  void WriteParameters();
} ;
#endif



