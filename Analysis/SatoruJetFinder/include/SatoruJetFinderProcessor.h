/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/*
** This file is part of the MarlinReco Project.
** Forming part of the SubPackage: SatoruJetFinder.
**
** For the latest version download from Web CVS:
** http://www-zeuthen.desy.de/lc-cgi-bin/cvsweb.cgi/marlinreco/?cvsroot=MarlinReco
**
** $Id$
**
**
*/

#ifndef SatoruJetFinderProcessor_h
#define SatoruJetFinderProcessor_h 1

#include "marlin/Processor.h"

#include "IMPL/LCCollectionVec.h"
#include "lcio.h"

using namespace lcio;

struct SatoruPartonArray {
  int NumberOfPartons;
  float Momentum[1200];
  int PointerParticleToJets[300];
};
struct SatoruJetsArray {
  int NumberOfJets;
  float Momentum[80];
};

extern "C" {
void syjkrn_(const char* GlobalModus_,
             // Primary jet finding
             int& NJetRequested_, float& Threshold_, int& PrimaryJetFindingMode_, float* YCut_,
             //  Reassosiation after first jet finding
             int& MergingMode_, float& MergingThreshold_,
             // Second Jet finding mode
             int& SecondJetFindingMode_,
             // input array--> array of
             // PPartons(DimensionOfInputArray(4-6),NumberOfPartons)
             int& NumberOfPartons_, int& DimensionOfInputArray_, float* PPartons_,
             // Output
             int& MaximalNumberOfJets_, int& NJetsFound_, int* PointerParicleToJet_, int& DimensionOfOutputArray_,
             float* PJets_, float& YMinus_, float& YPlus_,
             // error flag
             int& IError_, int GlobalModusLenght_);
}

namespace marlin {

/**
 * A universal jetfinder module developed by Satoru Yamashita
 *
 * This processor is a wrapper to the multi algorithm jet finder
 * code written by Satoru Yamashita for the OPAL collaboratoin.
 * For further details concerning this fortran code read the
 * OPAL Technical Note TN579.
 *
 * To compile this processor, you need a FORTRAN77 compiler like g77
 * and you have to link against the CERNLIB. If you use the MARLIN package
 * mechanism you have to add the following two lines to the 'userlib.gmk'
 * file in the MARLIN directory:
 * USERLIBS += -lg2c
 * USERLIBS += -L /opt/products/cernlib/pro/lib/ -lmathlib -lkernlib
 *
 * Steering file parameters:\n
 * Allways to specify:\n
 * \b InputCollection  (e.g. energy flow object )of type reconstructed particle
 * \b OutputCollection (e.g. jets )of type reconstructed particle
 * \b Mode
 *
 * At the moment 5 modes are implemented:
 *
 * \b DurhamNjet:
 * durham jetfinding with a fixed number of jets to be specify
 * with
 * NJetRequested
 *
 * \b DurhamYCut:
 * durham jetfinding with a fixed ycut to be specify with
 * YCut
 *
 * \b ConeBlanka
 * cone jet finder with cone threshold energy of 0.7 GeV and
 * R =0.2
 *
 * \b Saturo
 * first durham jetfinding, afterwards reassignment of the
 * objects to the jet axes found in this first iteration
 * with the jade sceme
 * number of jets to be specify with
 * NJetRequested
 *
 * \b Manual everything to be set by hand:
 * \verbatim
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
 *
 * @author Satoru Yamashita (original fortran code), Thorsten Kuhl, J&ouml;rgen Samson
 * @version $Id$
 */
class SatoruJetFinderProcessor : public Processor {

public:
  virtual Processor* newProcessor() { return new SatoruJetFinderProcessor; }

  SatoruJetFinderProcessor();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the module, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

private:
  std::string _inputCollection;
  std::string _outputCollection;
  std::string _successTag;
  bool _writeTag;
  std::string _jetFindingMode;
  SatoruPartonArray _partonsWorkArray;
  SatoruJetsArray _jetsWorkArray;
  //    int _nJetVal;
  std::string _globalMode;
  int _nJetRequested;
  float _threshold;
  int _primaryJetFindingMode;
  float _yCutParam;
  float _rConeParam;
  float _epsConeParam;
  float _yCut[2];
  int _mergingMode;
  float _mergingThreshold;
  int _secondJetFindingMode;
  int _debug;
  float _YMinus;
  float _YPlus;
  void putPartons(LCEvent* evt);
  void callSatoru(LCEvent* evt);
  void goSatoru(LCEvent* evt, LCCollection* JetsCol);
  void getJets(LCEvent* evt, LCCollection* JetsCol);

  /* debug routines */
  void writePartons();
  void writeJets();
  void writeParameters();

  // true constants for mode selection
  struct MD {
    const static int JADE_E0 = 1;
    const static int JADE_P = 2;
    const static int JADE_P0 = 3;
    const static int JADE_E = 4;
    const static int DURHAM = 5;
    const static int GENEVA = 6;
    const static int CAMBRIDGE = 7;
    const static int LUCLUS_1 = 8;
    const static int LUCLUS_2 = 9;
    const static int LUCLUS_3 = 10;
    const static int LUCLUS_4 = 11;
    const static int CONE = 12;
  };
};

} // namespace marlin
#endif
