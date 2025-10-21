// Boost Graph Library
#include <boost/graph/adjacency_list.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/graph/max_cardinality_matching.hpp>
#pragma GCC diagnostic pop

#include "EVENT/ReconstructedParticle.h"
#include "IMPL/LCCollectionVec.h"
#include "lcio.h"
#include "marlin/Processor.h"
#include <algorithm>
#include <bitset>
#include <cassert>
#include <gsl/gsl_cdf.h>
#include <iomanip>
#include <list>
#include <set>
#include <stdio.h>
#include <sys/time.h>
#include <vector>

enum { NVMAX = 100, NEMAX = 100 };
std::vector<std::vector<int>>
    comb; // Vector of vectors to store the possible edge choices for each combination of first vertices

using namespace lcio;

using namespace boost;

/** GammaGammaSolutionFinder:<br>
 *
 * @author Graham W. Wilson, University of Kansas
 */

class GammaGammaSolutionFinder : public marlin::Processor {

public:
  virtual marlin::Processor* newProcessor() { return new GammaGammaSolutionFinder; }

  GammaGammaSolutionFinder();

  /** Called at the beginning of the job before anything is read.
   *  Use to initialize the processor, e.g. book histograms.
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

  bool FindPFOs(LCEvent* evt);
  void FindGammaGammaSolutionZero(LCCollectionVec*);
  void FindGammaGammaSolutions(LCCollectionVec*);
  unsigned int CountIndependentPhotons();
  unsigned long int CombinatorialFactor(int n, int k);
  void FindCombinations(std::vector<std::vector<int>> array, unsigned int i, std::vector<int> accum);

  struct MyEdge {
    int u{}, v{};               // Edge from vertex u to vertex v with weight w
    std::bitset<NVMAX> vbits{}; // Bitset with the vertices used in this edge
    int pdgid{};                // PDG particle id (111 = pi0, 221 = eta, 331 = eta' )
    double edge_pvalue{};       // Chi-squared p-value for this edge being consistent with the mass constraint  (1 dof)
    double edge_chisq{};        // Chi-squared for this edge being consistent with the mass constraint
    double w{};                 // Edge weight
  };

  struct FirstVertex {           // struct for those vertices which correspond to min(u,v)
    int ivertex{};               // Vertex index from the initial graph definition
    int nedges{};                // Number of edges with which this vertex is the first vertex
    std::bitset<NEMAX> febits{}; // Bitset with the edge IDs available for this first vertex
    std::list<int> fvelist{};    // List with the edge IDs
    int fedge{};                 // First edge ID
  };

  struct CandidateSolution {
    int ne{};                   // Number of edges used
    int nv{};                   // Number of vertices used (should always be ne/2)
    std::bitset<NEMAX> ebits{}; // Bitset with the edges used in this solution
    double wsum{};              // Sum of the weights of the edges in the solution
    double pvalue{};            // Chi-squared p-value for wsum being consistent with ne edges (ne dof)
    int icomb{};                // Combination number
    int npi0{};                 // Number of pi0s in the solution
    int neta{};                 // Number of etas in the solution
    int netap{};                // Number of etaps in the solution
  };

private:
  std::vector<ReconstructedParticle*> _pfovec{};
  int _printing{};
  std::vector<std::string> _gammagammaCandidateCollections{};
  std::string _outputParticleCollectionName{};
  double _maxCombinationsCut{};
  int _nToRemove{}; // Number of edges less than maximal to consider
  int _algorithm{}; // Solution Finding Algorithm (1=Greedy, 2=Exhaustive)
  static bool PfoProbabilitySortFunction(ReconstructedParticle* lhs, ReconstructedParticle* rhs);

protected:
};
