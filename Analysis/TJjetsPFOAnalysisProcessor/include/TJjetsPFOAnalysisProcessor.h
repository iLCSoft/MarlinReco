#ifndef TJjetsPFOAnalysisProcessor_h
#define TJjetsPFOAnalysisProcessor_h 1

#include "marlin/Processor.h"
#include "StandardIncludes.h"
#include "TrueJet_Parser.h"

// // TrueJet_Parser can link generator level color singlets (quarks,gluons)
// // to the jet particles they produced, thereby allowing comparison of
// // correct jet with own clustered one.
// #include "Adjusted_TrueJet_Parser.h"


class TJjetsPFOAnalysisProcessor : public Processor , public TrueJet_Parser {
  public:
    virtual Processor*  newProcessor() { return new TJjetsPFOAnalysisProcessor ; }

    TJjetsPFOAnalysisProcessor() ;

    // These two lines avoid frequent compiler warnings when using -Weffc++
    TJjetsPFOAnalysisProcessor( const TJjetsPFOAnalysisProcessor& ) = delete;
    TJjetsPFOAnalysisProcessor& operator=( const TJjetsPFOAnalysisProcessor& ) = delete;

    /** Called at the begin of the job before anything is read.
     * Use to initialize the processor, e.g. book histograms.
     */
    void init() ;

    /** Called for every run.
     */
    void processRunHeader( LCRunHeader* run ) ;

    /** Called for every event - the working horse.
     */
    void processEvent( LCEvent * event ) ;

    /** For TrueJet_Parser -> see its documentation
    */
    std::string get_recoMCTruthLink(){ return _recoMCTruthLink  ; } ;

    void check( LCEvent * event ) ;


    /** Called after data processing for clean up.
     */
    void end() ;

  protected:
    /** Definitions using typedef
    */
    typedef std::pair< MCParticleVec, ReconstructedParticleVec > JetContentPair;
    typedef std::vector< JetContentPair* > JetVec;
    typedef std::vector<const ReconstructedParticle *> ParticleVector;
    typedef std::vector<const MCParticle*> MCParticleVector;
    typedef std::set<EVENT::MCParticle*> MCParticleList;
    typedef std::set<EVENT::MCParticle*> MCParticleSet;
    typedef std::vector<std::string> StringVector;
    typedef std::vector<LCRelationNavigator*> LCRelationNavigatorVec;

    JetVec m_jets{};

    /** Input collection name.
    */
    std::string _colAllPFOs{};
    std::string _colMC{};

    // Input collections for TrueJet
    std::string _recoMCTruthLink{};


    /** Output File specifics
    */
    TFile* _otfile{};
    TTree* _datatrain{};


    std::string _rootfilename{};

    template <class Object> bool areDisjointVectors( const std::vector<Object> &v1, const std::vector<Object> &v2 ) const;
    template <class Object> std::vector<Object> turnVectorOfConstToVector( const std::vector<const Object> &v ) const;
    void makeNTuple() ;
    void findTrueJetParticles( LCEvent* event );
    bool hasSomeParentsInMCList( EVENT::MCParticle *pMCParticle, MCParticleList &mcs) const;

    /** PFOAnalysis stuff
    */

    void Clear();

    /**
    *  @brief  Extract lcio collections
    *
    *  @param  event the lc event
    */
    void ExtractCollections(JetContentPair* jet_content);

    /**
    *  @brief  Apply pfo selection rules, starting with root particles
    *
    *  @param  pMCParticle the address of a mc particle (initially call this recursive function with a top-level, root particle)
    *  @param  mcPfoCandidates to collect the list of mc pfo candidates
    */
    void ApplyPfoSelectionRules(EVENT::MCParticle *pMCParticle, MCParticleList &mcPfoCandidates) const;

    /**
    *  @brief  Make quark variables
    *
    *  @param  event the lc event
    */
    void MakeQuarkVariables(JetContentPair* jet_content);

    /**
    *  @brief  Perform pfo analysis
    */
    void PerformPfoAnalysis();

    /**
    *  @brief  Sort mc pfo targets by decreasing energy
    *
    *  @param  pLhs the address of the first mc pfo target
    *  @param  pRhs the address of the second mc pfo target
    */
    static bool SortPfoTargetsByEnergy(const EVENT::MCParticle *const pLhs, const EVENT::MCParticle *const pRhs);

    int                 m_nRun{};                                 ///<
    int                 m_nEvt{};                                 ///<
    int                 m_nJet{};                                 ///<

    int                 m_jetInitElPDG{};                         ///<
    int                 m_jetFinElPDG{};                          ///<

    int                 m_nRunSum{};                              ///<
    int                 m_nEvtSum{};                              ///<

    std::string         m_inputPfoCollection{};                   ///<
    std::string         m_mcParticleCollection{};                 ///<

    int                 m_printing{};                             ///<
    std::string         m_rootFile{};                             ///<

    int                 m_lookForQuarksWithMotherZ{};             ///<

    float               m_mcPfoSelectionRadius{};                 ///<
    float               m_mcPfoSelectionMomentum{};               ///<
    float               m_mcPfoSelectionLowEnergyNPCutOff{};      ///<

    ParticleVector      m_pfoVector{};                            ///<
    MCParticleVector    m_pfoTargetVector{};                      ///<

    int                 m_nPfosTotal{};                           ///<
    int                 m_nPfosNeutralHadrons{};                  ///<
    int                 m_nPfosPhotons{};                         ///<
    int                 m_nPfosTracks{};                          ///<
    float               m_pfoEnergyTotal{};                       ///<
    float               m_pfoEnergyNeutralHadrons{};              ///<
    float               m_pfoEnergyPhotons{};                     ///<

    float               m_pfoEnergyTracks{};                      ///<
    float               m_pfoECalToEmEnergy{};                    ///<
    float               m_pfoECalToHadEnergy{};                   ///<
    float               m_pfoHCalToEmEnergy{};                    ///<
    float               m_pfoHCalToHadEnergy{};                   ///<
    float               m_pfoMuonToEnergy{};                      ///<
    float               m_pfoOtherEnergy{};                       ///<

    float               m_pfoMassTotal{};                         ///<

    typedef std::vector<float> FloatVector;
    FloatVector         m_pfoEnergies{};                          ///<
    FloatVector         m_pfoPx{};                                ///<
    FloatVector         m_pfoPy{};                                ///<
    FloatVector         m_pfoPz{};                                ///<
    FloatVector         m_pfoCosTheta{};                          ///<

    FloatVector         m_pfoTargetEnergies{};                    ///<
    FloatVector         m_pfoTargetPx{};                          ///<
    FloatVector         m_pfoTargetPy{};                          ///<
    FloatVector         m_pfoTargetPz{};                          ///<
    FloatVector         m_pfoTargetCosTheta{};                    ///<

    typedef std::vector<int> IntVector;
    IntVector           m_pfoPdgCodes{};                          ///<
    IntVector           m_pfoTargetPdgCodes{};                    ///<

    int                 m_nPfoTargetsTotal{};                     ///<
    int                 m_nPfoTargetsNeutralHadrons{};            ///<
    int                 m_nPfoTargetsPhotons{};                   ///<
    int                 m_nPfoTargetsTracks{};                    ///<

    float               m_pfoTargetsEnergyTotal{};                ///<
    float               m_pfoTargetsEnergyNeutralHadrons{};       ///<
    float               m_pfoTargetsEnergyPhotons{};              ///<
    float               m_pfoTargetsEnergyTracks{};               ///<

    float               m_mcEnergyENu{};                          ///<
    float               m_mcEnergyFwd{};                          ///<
    float               m_eQQ{};                                  ///<
    float               m_eQ1{};                                  ///<
    float               m_eQ2{};                                  ///<
    float               m_costQQ{};                               ///<
    float               m_costQ1{};                               ///<
    float               m_costQ2{};                               ///<
    float               m_mQQ{};                                  ///<
    float               m_thrust{};                               ///<
    int                 m_qPdg{};                                 ///<

    TFile              *m_pTFile{};                               ///<
    TTree              *m_pTTree{};                               ///<
    TH1F               *m_hPfoEnergySum{};                        ///<
    TH1F               *m_hPfoEnergySumL7A{};                     ///<

} ;

// Include other headers with definition of template functions at end
#include "VectorHelper.h"
#endif
