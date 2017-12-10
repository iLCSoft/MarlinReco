#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "GammaGammaSolutionFinder.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

typedef CLHEP::HepLorentzVector LorentzVector ;
typedef CLHEP::Hep3Vector Vector3D ;

// Marlin stuff
#include <marlin/Global.h>
// the event display

// ROOT stuff
#include "TMath.h"
#include "TMatrixD.h"

#include <cstdlib>

using namespace lcio;

GammaGammaSolutionFinder aGammaGammaSolutionFinder;

GammaGammaSolutionFinder::GammaGammaSolutionFinder() : marlin::Processor("GammaGammaSolutionFinder") {

  registerProcessorParameter( "Printing" , 
			      "Print certain messages"  ,
			      _printing,
			       (int)1 ) ;

  std::vector<std::string> gammagammaCandidateCollections;
  gammagammaCandidateCollections.push_back(std::string("GammaGammaCandidatePi0s"));
  gammagammaCandidateCollections.push_back(std::string("GammaGammaCandidateEtas"));
  gammagammaCandidateCollections.push_back(std::string("GammaGammaCandidateEtaPrimes"));
  registerInputCollections( LCIO::RECONSTRUCTEDPARTICLE,
                            "GammaGammaCandidateCollections" ,
                            "Gamma-Gamma Candidate Collection Names" ,
                            _gammagammaCandidateCollections ,
                            gammagammaCandidateCollections);

  std::string outputParticleCollectionName = "GammaGammaSolutions";
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "OutputParticleCollectionName" , 
			     "Output Particle Collection Name "  ,
			     _outputParticleCollectionName,
			     outputParticleCollectionName);

  registerProcessorParameter( "MaxCombinationsCut" , 
                              "Maximum Number of Explorable Combinations"  ,
                              _maxCombinationsCut,
                               (double)100000000.0) ;

  registerProcessorParameter( "LessThanMaximalCut" , 
                              "Maximum Number of Fewer GammaGammaParticles"  ,
                              _nToRemove,
                               (int)0) ;            // (set to huge number - like 999 if want to explore all)

  registerProcessorParameter( "SolutionFindingAlgorithm" , 
                              "Solution Finding Algorithm"  ,
                              _algorithm,
                               (int)2) ;            

  return;

}

//===================================================================================

void GammaGammaSolutionFinder::init() { 
  if(_printing>1)printParameters();
  return;
}

//===================================================================================

void GammaGammaSolutionFinder::processRunHeader( LCRunHeader*  /*run*/) { 
  return;
}

//===================================================================================

void GammaGammaSolutionFinder::processEvent( LCEvent * evt ) { 

  // Make a new vector of particles
  LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  recparcol->setSubset(true);                      

  // Access PFO collection
  bool found = this->FindPFOs(evt);
  if(found){
    if(_printing>1)std::cout << "Find GammaGamma Solution: " << std::endl;
    unsigned int nphotons=this->CountIndependentPhotons();
    unsigned int nedges=_pfovec.size();

    if(_printing>1)std::cout << "nphotons = " << nphotons << " nedges = " << nedges << std::endl;
    unsigned long int nCr = this->CombinatorialFactor(nedges,nphotons/2);
    if(_printing>1)std::cout << "Naive combinatorial factor of " << std::setw(22) << nCr << std::endl;

    if(_algorithm == 1 ){
// Greedy algorithm
       this->FindGammaGammaSolutionZero(recparcol);        
    }
    else if(_algorithm == 2){
       if(double(nCr) < _maxCombinationsCut){
          if(_printing>1)std::cout << " Pursuing full combinatoric exploration " << std::endl; 
          this->FindGammaGammaSolutions(recparcol);
       }
       else{
          if(_printing>1)std::cout << " Too many combinations - assigning quick greedy solution " << std::endl;
          this->FindGammaGammaSolutionZero(recparcol);
       }
    }
    else{
       std::cout << "Unsupported solution finding algorithm " << _algorithm << std::endl;
    }
  }

  // Add new collection to event
  evt->addCollection( recparcol , _outputParticleCollectionName.c_str() );
  
  return;
  
}

//===================================================================================

void GammaGammaSolutionFinder::end(){ 
  return;
}

unsigned long int GammaGammaSolutionFinder::CombinatorialFactor(int n, int k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    unsigned long int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

//===================================================================================

void GammaGammaSolutionFinder::FindCombinations(std::vector<std::vector<int> > array, unsigned int i, std::vector<int> accum)
{
    if (i == array.size()) // done, no more rows
    {
        comb.push_back(accum); // assuming comb is global
    }
    else
    {
        std::vector<int> row = array[i];                     // NB   Each row of array may have a different length
        for(unsigned int j = 0; j < row.size(); ++j)
        {
            std::vector<int> tmp(accum);
            tmp.push_back(row[j]);
            FindCombinations(array,i+1,tmp);
        }
    }
}


bool GammaGammaSolutionFinder::FindPFOs( LCEvent* evt ) {

// Add all 3 GammaGammaCandidate collections to the PFO vector (pi0s , etas, eta's)

  bool tf = false;

  // clear old vector
  _pfovec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;


// All GammaGammaCandidates
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){
     for (unsigned int j=0; j < _gammagammaCandidateCollections.size(); ++j) {    
        if(*name==_gammagammaCandidateCollections[j]){
           LCCollection* col = evt->getCollection(*name);
           unsigned int nelem = col->getNumberOfElements();
           tf = true;
           for(unsigned int i=0;i<nelem;i++){
	      ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	      _pfovec.push_back(recoPart);
           }
        }
     }
  }

  if(_printing>1)std::cout << "FindPFOs : (nPFOs = " << _pfovec.size() << " )" << std::endl; 

  if(_printing>1)std::cout << "Find PFOs : " << tf << std::endl; 

  return tf;
  
}

// ===================================================================================
// Initial Greedy Algorithm
void GammaGammaSolutionFinder::FindGammaGammaSolutionZero(LCCollectionVec * recparcol) {

  if(_printing>1)std::cout << "FindGammaGammaSolutionZero : (nPFOs = " << _pfovec.size() << " )" << std::endl;

  typedef std::set<const ReconstructedParticle*> PfoSet;
  PfoSet particles_assigned;    // Set with pointers to the daughter particles (photons for now) that are already used in the solution

  // For convenience sort the GammaGammaCandidates by fit probability
  std::sort(_pfovec.begin(),_pfovec.end(),GammaGammaSolutionFinder::PfoProbabilitySortFunction);

  // Algorithm 1. Based on fit probability ordering, add GammaGammaCandidates to the output GammaGammaParticles collection 
  // if they do not contain already assigned constituent particles. For now, the constituents are always two photons  
  // - but may not hurt to keep this generic for when other decay modes (eg. Dalitz decay) are included
  for(unsigned int i = 0; i < _pfovec.size(); i++){
      const ReconstructedParticleVec particles = _pfovec[i]->getParticles();
      if(_printing>3)std::cout << "FindGammaGammaSolutionZero: (nparticles = " << particles.size() << " )" << std::endl;
      bool assignable = true;   
      for(unsigned int j = 0; j < particles.size(); j++){           // TODO - maybe double check these are indeed photons ..
          const ReconstructedParticle *particle = particles[j];
          if( particles_assigned.find(particle) != particles_assigned.end() )assignable = false; // If particle is already assigned then it is not still assignable
      }
      if(assignable){
         recparcol->addElement(_pfovec[i]);                    // Add meson to the output collection if all constituent particles are still assignable
         for(unsigned int j = 0; j < particles.size(); j++){
             const ReconstructedParticle *particle = particles[j];
             particles_assigned.insert(particle);              // Update set of particles already assigned in this solution
         }
      }         
  }
  return;
}

// ===================================================================================
// Count how many independent daughter photons are available in the GammaGammaParticle collections
unsigned int GammaGammaSolutionFinder::CountIndependentPhotons() {

  if(_printing>1)std::cout << "Count Independent Photons : (nPFOs = " << _pfovec.size() << " )" << std::endl;

  typedef std::set<const ReconstructedParticle*> PfoSet;
  PfoSet daughter_particles;    // Set with all the daughter particles

  // For convenience sort the GammaGammaCandidates by fit probability
  std::sort(_pfovec.begin(),_pfovec.end(),GammaGammaSolutionFinder::PfoProbabilitySortFunction);

  int k=-1;

  for(unsigned int i = 0; i < _pfovec.size(); i++){
      const ReconstructedParticleVec particles = _pfovec[i]->getParticles();
//      if(_printing>3)std::cout << "FindGammaGammaSolutions: (nparticles = " << particles.size() << " )" << std::endl; 
      for(unsigned int j = 0; j < particles.size(); j++){           // TODO - maybe double check these are indeed photons ..
          const ReconstructedParticle *particle = particles[j];
// GWW DEBUG
//          std::cout << "GWWWW " << i << " " << j << " " << particle << std::endl;
          if( daughter_particles.find(particle) == daughter_particles.end() ){
              daughter_particles.insert(particle);                              // Add particle to the set
              k++;                                                              // Increment photon index
              //std::cout << i << " " << k << std::endl;
          }
      }       
  }

  if(_printing>1)std::cout << "Number of independent photons available = " << daughter_particles.size() << std::endl;

  unsigned int nphotons = daughter_particles.size();

  return nphotons;
}

//===================================================================================

void GammaGammaSolutionFinder::FindGammaGammaSolutions(LCCollectionVec * recparcol) {

   struct timeval tv;
   gettimeofday(&tv,NULL);
   double t0 = tv.tv_sec+(tv.tv_usec/1000000.0);

   std::vector<MyEdge> evec;               // Vector of edges available
   evec.reserve(NEMAX*25);
   std::vector<CandidateSolution> svec;    // Vector with candidate solutions
   svec.reserve(1000);
   std::vector<FirstVertex> fvec;          // Vector with first vertices
   fvec.reserve(NVMAX*25);
   int u=-1;
   int v=-1;
   double edge_chisq;
// TODO Add some checks that the bitsets are big enough ...
   std::bitset<NVMAX> vbits;
   std::bitset<NVMAX> ored_vbits;
   std::bitset<NEMAX> ebits;
   std::bitset<NEMAX> febits;
   unsigned long int nCr;
   unsigned long int nCr2;

  if(_printing>1)std::cout << "FindGammaGammaSolutions : (nPFOs = " << _pfovec.size() << " )" << std::endl;

  typedef std::set<const ReconstructedParticle*> PfoSet;
  PfoSet daughter_particles;    // Set with all the daughter particles
  PfoSet particles_assigned;    // Set with pointers to the daughter particles (photons for now) that are already used in the solution

  // For convenience sort the GammaGammaCandidates by fit probability
  std::sort(_pfovec.begin(),_pfovec.end(),GammaGammaSolutionFinder::PfoProbabilitySortFunction);

  // First pass through - keep track of the photons involved - add each photon once to the daughter_particles set
  int k = -1;

  for(unsigned int i = 0; i < _pfovec.size(); i++){
      const ReconstructedParticleVec particles = _pfovec[i]->getParticles();
      if(_printing>3)std::cout << "FindGammaGammaSolutions: (nparticles = " << particles.size() << " )" << std::endl; 
      for(unsigned int j = 0; j < particles.size(); j++){           // TODO - maybe double check these are indeed photons ..
          const ReconstructedParticle *particle = particles[j];
// GWW DEBUG
          if(_printing>7)std::cout << "GWWWW " << i << " " << j << " " << particle << std::endl;
          if( daughter_particles.find(particle) == daughter_particles.end() ){
              daughter_particles.insert(particle);                              // Add particle to the set
              k++;                                                              // Increment photon index
              if(_printing>3)std::cout << "New photon " << particle << " " << k << std::endl;
          }
      }       
  }

  // See what we have in the set 
  // (Note that since these are pointers they are by default ordered by their address in memory).
  if(_printing>3)std::cout << "Set size " << daughter_particles.size() << std::endl;
  for(std::set<const ReconstructedParticle*>::iterator it=daughter_particles.begin(); it!=daughter_particles.end(); ++it){
      if(_printing>3)std::cout << *it << std::endl;
  }

// At this point we basically know the number of photons and GammaGammaCandidates

  if(_printing>5)std::cout << "Graph definition file information Dry-Run" << std::endl;
  if(_printing>5)std::cout << std::setw(4) << daughter_particles.size() << std::setw(4) << _pfovec.size() << std::endl;
  std::set<const ReconstructedParticle*>::iterator it;
  std::pair<std::set<const ReconstructedParticle*>::iterator,bool> retcode;

// Second pass through. This time each photon should already be in the set.
  for(unsigned int i = 0; i < _pfovec.size(); i++){
      const ReconstructedParticleVec particles = _pfovec[i]->getParticles();
      if(_printing>5)std::cout << "GammaGammaCandidate " << std::setw(4) << i << std::endl; 
      for(unsigned int j = 0; j < particles.size(); j++){           // TODO - maybe double check these are indeed photons ..
          const ReconstructedParticle *particle = particles[j];

          retcode = daughter_particles.insert(particle);    // insert should fail as all photons should already be in the set
          if (retcode.second) {
              if(_printing>5)std::cout << particle << " inserted as element ";
          }
          else {
              if(_printing>5)std::cout << particle << " already exists as element ";
          }
          if(_printing>5)std::cout << std::distance(daughter_particles.begin(),retcode.first) << std::endl;
      }       
  }
// Code here is designed to both write out a version of the problem that can be analyzed stand-alone 
// and to construct the information needed for defining the graph problem.

  if(_printing>6)std::cout << "GRAPH definition file information" << std::endl;
  int nvertices = daughter_particles.size();
  int nedges = _pfovec.size();
  typedef adjacency_list<vecS, vecS, undirectedS> my_graph;
  my_graph g(nvertices);
  int vertex_degree[NVMAX];
  for (int i = 0; i < NVMAX; ++i) {
      vertex_degree[i] = 0;
  }
  if(_printing>6)std::cout << std::setw(4) << daughter_particles.size() << std::setw(4) << _pfovec.size() << std::endl;
// Third pass through with more terse information suitable for further analysis
  for(unsigned int i = 0; i < _pfovec.size(); i++){
      const ReconstructedParticleVec particles = _pfovec[i]->getParticles();
      if(_printing>6)std::cout << std::setw(4) << i; 
      for(unsigned int j = 0; j < particles.size(); j++){           // TODO - maybe double check these are indeed photons ..
          const ReconstructedParticle *particle = particles[j];
          retcode = daughter_particles.insert(particle);    // insert should fail as all photons should already be in the set
          if(_printing>6)std::cout << "  " << std::setw(4) << std::distance(daughter_particles.begin(),retcode.first) ;
          if(j==0)u = std::distance(daughter_particles.begin(),retcode.first);
          if(j!=0)v = std::distance(daughter_particles.begin(),retcode.first);
      }
      const int particleType  = _pfovec[i]->getType();
      // GoodnessOfPID currently filled with fit probability 
      // - may at some point be more convenient to fill with the fit chi-squared ...
      const float GoodnessOfPID = _pfovec[i]->getGoodnessOfPID();
      if(_printing>6)std::cout << "  " << std::setw(4) << particleType << "  " << 
                     std::fixed << std::setw(10) << std::setprecision(6) << GoodnessOfPID << std::endl;

      evec.push_back(MyEdge());
      evec[i].u = std::min(u,v);
      evec[i].v = std::max(u,v);
      evec[i].pdgid = particleType;
      evec[i].edge_pvalue = double(GoodnessOfPID);
      edge_chisq = gsl_cdf_chisq_Qinv(double(GoodnessOfPID), 1.0);
      evec[i].edge_chisq = edge_chisq;
      evec[i].w = edge_chisq;
      vbits = 0;            // Initialize all bits in the bitset to 0
      vbits.flip(u);        // Set appropriate bit in the bitset for vertex u
      vbits.flip(v);        // Set appropriate bit in the bitset for vertex v
      evec[i].vbits = vbits;
      vertex_degree[u] = vertex_degree[u]+1;
      vertex_degree[v] = vertex_degree[v]+1;
// BGL graph
      add_edge(std::min(u,v),std::max(u,v),g);
  }

// Lots of print-out specifying the problem

   if(_printing>6)std::cout << " " << std::endl;
   if(_printing>6)std::cout << "Edge summary " << std::endl;
   for (int i = 0; i < nedges; ++i) {
       if(_printing>6)std::cout << std::setw(3) << i << " " << evec[i].vbits << " " 
                                << std::setw(2) << evec[i].u << " " << std::setw(2) << evec[i].v << std::endl;
   }
   if(_printing>6)std::cout << " " << std::endl;
   if(_printing>6)std::cout << "Vertex degrees " << std::endl;
   int vdsum = 0;
   for (int i = 0; i < nvertices; ++i){
       if(_printing>6)std::cout << std::setw(3) << i << " " << std::setw(3) << vertex_degree[i] << std::endl;
       vdsum += vertex_degree[i];
   }
   if(_printing>6)std::cout << "vdsum = " << vdsum << std::endl;
 
   if(_printing>6)std::cout << "Vertex - Vertex partnerships" << std::endl;
   for(int iv =0; iv < nvertices; ++iv){
       if(_printing>6)std::cout << "Vertex " << std::setw(2) << iv << ": " ;
       for (int je = 0; je < nedges; ++je) {
           if(evec[je].u == iv){
              if(_printing>6)std::cout << std::setw(2) << evec[je].v << " ";
           }
       }
       if(_printing>6)std::cout << std::endl;
   }

   if(_printing>6)std::cout << "Vertex - Edge partnerships" << std::endl;
   for(int iv =0; iv < nvertices; ++iv){
       if(_printing>6)std::cout << "Vertex " << std::setw(2) << iv << ": " ;
       febits = 0;
       int nfvedges = 0;
       for (int je = 0; je < nedges; ++je) {
           if(evec[je].u == iv){
              if(_printing>6)std::cout << std::setw(2) << je << " ";
              nfvedges ++;
              febits.flip(je);        // Set appropriate bit in the bitset for this first vertex u
           }              
       }
       if(nfvedges>0){
          fvec.push_back(FirstVertex());
          fvec.back().ivertex = iv;
          fvec.back().nedges = nfvedges;
          fvec.back().febits = febits;
          std::list<int> fvelist;
          for (int je = 0; je < nedges; ++je) {
              int nfound = 0;
              if(febits.test(je)){
                 fvelist.push_back(je);
                 nfound++;
                 if(nfound==1)fvec.back().fedge = je;
              }
          }
          fvec.back().fvelist = fvelist;          
       }
       if(_printing>6)std::cout << std::endl;
   }
   if(_printing>6)std::cout << " " << std::endl;
   if(_printing>6)std::cout << "Check FirstVertex structure" << std::endl;
   if(_printing>6)std::cout << std::setw(2) << fvec.size() << " FirstVertices out of " 
                            << std::setw(2) << nvertices << " vertices in total" << std::endl;
   if(_printing>6)std::cout << " i iv ne           febits                 Edges" << std::endl;

   int nfirst_vertices = 0;
   for(unsigned int i = 0; i < fvec.size(); ++i) {
       nfirst_vertices++;     
       if(_printing>6)std::cout << std::setw(2) << i << " " << std::setw(2) << fvec[i].ivertex << " " 
                                << std::setw(2) << fvec[i].nedges << " " << fvec[i].febits << "  { "; 
       for (std::list<int>::iterator it=fvec[i].fvelist.begin(); it!=fvec[i].fvelist.end(); ++it){
           if(_printing>6)std::cout << std::setw(2) << *it << " ";
       } 
       if(_printing>6)std::cout << "}" << std::endl;
   }
   if(_printing>6)std::cout << "nfirst_vertices = "<<nfirst_vertices << std::endl;

// Graph analysis by boost
   std::vector<graph_traits<my_graph>::vertex_descriptor> mate(nvertices);
  // find the maximum cardinality matching. we'll use a checked version
  // of the algorithm, which takes a little longer than the unchecked
  // version, but has the advantage that it will return "false" if the
  // matching returned is not actually a maximum cardinality matching
  // in the graph.
   bool success = checked_edmonds_maximum_cardinality_matching(g, &mate[0]);
   assert(success);
   if(_printing>6)std::cout << std::endl << "Found a matching of size " << matching_size(g, &mate[0]) << std::endl;
   if(_printing>6)std::cout << "The matching is:" << std::endl;
   graph_traits<my_graph>::vertex_iterator vi, vi_end;
   for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi)
      if (mate[*vi] != graph_traits<my_graph>::null_vertex() && *vi < mate[*vi])
         if(_printing>6)std::cout << "{" << *vi << ", " << mate[*vi] << "}" << std::endl;
   if(_printing>6)std::cout << std::endl;
// Probably wouldn't hurt to check that this is a matching ......

   int n = nfirst_vertices;             // n: Number of edges to choose from.
   int nmatch = 0;
   const unsigned int rmax = matching_size(g, &mate[0]);
   if(_printing>2)std::cout << "Setting rmax to " << rmax << std::endl;

   int rminA = 1;
   int rminB = rmax - _nToRemove;
   int rminC;
   rminC = rminA;
   if(rminB>rminA)rminC = rminB;
   const unsigned int rmin = rminC;

//   const int rmin = std::max(1,rmax - _nToRemove);
   if(_printing>2)std::cout << "Setting rmin to " << rmin << std::endl; 

   gettimeofday(&tv,NULL);
   double t1 = tv.tv_sec+(tv.tv_usec/1000000.0);
   double t2 = 0.0;
   if(_printing>3)std::cout << "Initialization took " << std::fixed << std::setw(12) << std::setprecision(6) << t1-t0 << " (s) " << std::endl;

// Now the initial combination finding is with respect to first-vertex edges.
   for (unsigned int r=rmax; r>=rmin; --r){                // Try various size matchings 
 
   if(_printing>2)std::cout << "subsidiary matching problem: find " << r << " edges from the set of size " << n << " with no overlaps " << std::endl;
   nCr =  this->CombinatorialFactor(nedges,r);     // Choose r edges from nedges
   nCr2 = this->CombinatorialFactor(n,r);          // Choose r first-vertices
   if(_printing>2)std::cout << "combination counting ... nfv = " << n << " r = " << std::setw(2) << r 
                            << "  nCr = " << std::setw(22) << nCr << "  nCr (fv) = " << std::setw(22) << nCr2 << std::endl;
   if(_printing>2)std::cout << std::endl;
   int nmatch_withr = 0; 

   unsigned long int icombi=0;

// Check each combination with r edges. Use nextpermutation method from StackOverflow for generating combinations
   std::vector<bool> vb(n);   
   std::vector<int> d;                               // Vector to store integers denoting the first vertices in the potential match
   std::fill(vb.begin() + r, vb.end(), true);        // fill elements vb[r] .... vb[n-1] with true -- other elements should be false
   do {
       for (int i = 0; i < n; ++i) {
           if (!vb[i]) {
              d.push_back(i);                        // Keep track of the indices in vector d where the bool is false
           }
       }

// Probably best to create a vector of vectors with all the combinations of edges associated with each first vertex, and loop through those.
// http://stackoverflow.com/questions/8620030/generate-all-combination-of-elements-in-2d-vector may be helpful

       std::vector<int> c;                               // Vector to store integers denoting the edges in the potential match
       c.reserve(NEMAX);

// Form array with edge possibilities when nedges > 1
       std::vector< std::vector<int> > array;
       int nfilled = 0;       
       for(unsigned int i=0; i<r; ++i) { 
//           std::cout << " First vertex " << i << ": ";
           std::vector<int> tmp;
           if(fvec[d[i]].nedges>1){    // Only do this for those first vertices with more than one edge
// Iterate over the elements of the list associated with first vertex d[i]
              for (std::list<int>::iterator it=fvec[d[i]].fvelist.begin(); it!=fvec[d[i]].fvelist.end(); ++it){
//                std::cout << std::setw(2) << *it << " ";
                  tmp.push_back(*it);
              }
//            std::cout << std::endl;
              array.push_back(tmp);
           }
           else{
              c.push_back(fvec[d[i]].fedge);               // Should be OK at front of vector ??
              nfilled++;
           }
       }

  // Generate all combinations consisting of 1 element from each row
  std::vector<int> accum;             // empty on input
  if(!comb.empty())comb.clear();      // Need to also clear the output array - otherwise it just keeps accumulating.
  this->FindCombinations(array,0,accum);        // recursive algorithm - start with 0
  for ( std::vector<std::vector<int> >::size_type i = 0; i < comb.size(); ++i )
  {
       icombi++;
//       if(icombi%1000000==0)std::cout << "C " << std::setw(16) << icombi << std::endl;
       if(c.size()==r && comb[i].size()>0){
          c.erase(c.begin()+nfilled,c.end());                // Remove the entries at the end.
       }
       for ( std::vector<int>::size_type j = 0; j < comb[i].size(); ++j )
       {
           c.push_back(comb[i][j]);                                       // May need to be more careful about potential reallocation
       }

// Check to see how many vertices are used in this edge set
       ored_vbits = evec[c[0]].vbits;
       for (unsigned int i=1; i<r; ++i){
            ored_vbits = (ored_vbits | evec[c[i]].vbits);
            if(ored_vbits.count()!=2*(i+1))break;
       }
       size_t nused = ored_vbits.count();


       if(nused==2*r){             // Number of vertices used equals twice the number of requested edges ( = MATCH )

          ebits = 0;          // Initialize all bits in the bitset to 0

          double wsum = 0.0;
          int npi0 = 0;
          int neta = 0;
          int netap = 0;
          for (unsigned int i=0; i<r; ++i){
/*              std::cout << "Assigned edge details" << std::endl;
              std::cout << std::setw(2) << c[i] << " " << evec[c[i]].vbits << " " << std::setw(2) << evec[c[i]].u << " " 
                        << std::setw(2) << evec[c[i]].v << " " << evec[c[i]].pdgid
                        << std::fixed << std::setw(10) << std::setprecision(5) << evec[c[i]].w << " " 
                        << std::fixed << std::setw(10) << std::setprecision(6) << evec[c[i]].edge_pvalue << std::endl;   */
              ebits.flip(c[i]);                 // Set appropriate bit in the bitset for edge c[i]
              wsum += evec[c[i]].w;
              if(evec[c[i]].pdgid==111)npi0++;
              if(evec[c[i]].pdgid==221)neta++;
              if(evec[c[i]].pdgid==331)netap++;
          }
// Store the CandidateSolution for later perusal
          svec.push_back(CandidateSolution());
          svec.back().ne = r;
          svec.back().nv = nused;
          svec.back().wsum = wsum;
          svec.back().ebits = ebits;
          svec.back().npi0 = npi0;
          svec.back().neta = neta;
          svec.back().netap = netap;
// Not really needed to calculate the pvalue - but this was measured to be fairly fast (0.1s degradation for 47000 matches).
          double dof = double(r);
          double chisq_value = wsum;
          double pvalue = gsl_cdf_chisq_Q(chisq_value,dof);
          svec.back().pvalue = pvalue;
          svec.back().icomb = icombi;
/*          std::cout << "Candidate solution" << std::endl;
          std::cout << std::setw(2) << r << " " << ored_vbits << " " 
                    << std::setw(2) << npi0 << " " << std::setw(2) << neta << " " << std::setw(2) << netap << " " 
                    << std::fixed << std::setw(10) << std::setprecision(5) << wsum << " " 
                    << std::fixed << std::setw(10) << std::setprecision(6) << pvalue << std::endl; */
          nmatch++;
          nmatch_withr++;
       }
  }

       if(!d.empty())d.clear();                              // Clear vector d
       if(!array.empty())array.clear();                      // Clear vector of vectors array

       icombi++;                                             // Increment combination index
   } while ( std::next_permutation(vb.begin(), vb.end() ) );  

   if(_printing>2)std::cout << "Found " << std::setw(7) << nmatch_withr << " matches with " 
                            << std::setw(3) << r << " edges from " << std::setw(22) << nCr << " combinations " << std::endl;

   gettimeofday(&tv,NULL);
   t2 = tv.tv_sec+(tv.tv_usec/1000000.0);
   if(_printing>3)std::cout << "Match search for n = " << std::setw(2) << n << " r = " << std::setw(2) << r << " " 
                            << std::fixed << std::setw(12) << std::setprecision(6) << t2-t1 << " (s) " 
                            << std::fixed << std::setw(12) << std::setprecision(6) << 1.0e6*(t2-t1)/double(nCr) 
                            << " micro-seconds per combination " << std::endl;
   t1 = t2;

 }   // r loop

   if(_printing>2)std::cout << "Found " << std::setw(7) << nmatch << " matches total " << std::endl;

   if(_printing>5)std::cout << "Two particle hypothesis results (pi0, eta) with rmax = " << rmax << std::endl;
// hard-coded min and max searches (could probably do more elegantly ...)
   int ibest = -1;
   int imin;
   int imax;
   for (unsigned int r=rmax; r>0; --r){
       double min_wsum =  999999.0;
       double max_wsum = -999999.0;
       imin = -1;
       imax = -1;
//       std::cout << "Test r = " << r << std::endl;
       for(unsigned int i = 0; i < svec.size(); ++i) {
           if(svec[i].ne == int(r) && svec[i].netap ==0){
//              std::cout << "Testing i = " << i << " " << svec[i].wsum << std::endl;
              if(svec[i].wsum > max_wsum){
                 max_wsum = svec[i].wsum;
                 imax = i;
              }

              if(svec[i].wsum < min_wsum){
                 min_wsum = svec[i].wsum;
                 imin = i;
              }
//              std::cout << "Sums etc " << imin << " " << min_wsum << " " << imax << " " << max_wsum << std::endl;
           }
       }
//       std::cout << " Checking solutions " << r << " " << imin << " " << imax << std::endl;
       if(r==rmax && imin >= 0)ibest = imin;              // Flag this as the "best solution".
       if(imin>=0 && imax >=0)
       if(_printing>5)std::cout << "wsum extrema for r = " << std::setw(2) << r << " " << std::setw(10) << min_wsum << " " << std::setw(10) << max_wsum 
                                << " pvalues " << std::setw(12) << svec[imin].pvalue << "  " << std::setw(12) << svec[imax].pvalue 
                                << "  ( " << imin << " " << imax << " ) " << " Combinations " << svec[imin].icomb << " " << svec[imax].icomb << std::endl;
   }
// Examine best solution
   if(ibest>=0){
      double pvalue = gsl_cdf_chisq_Q(svec[ibest].wsum,svec[ibest].ne);
      if(_printing>4)std::cout << "Best candidate solution " << ibest << std::endl;
      if(_printing>4)std::cout << "npi0, neta, netap " << svec[ibest].npi0 << " " << svec[ibest].neta << " " << svec[ibest].netap << std::endl;
      if(_printing>4)std::cout << std::setw(4) << ibest << " " << std::setw(2) << svec[ibest].ne << " " 
                               << std::setw(2) << svec[ibest].nv << "   " << svec[ibest].ebits << " " 
                               << std::fixed << std::setw(10) << std::setprecision(5) << svec[ibest].wsum << " " 
                               << std::fixed << std::setw(10) << std::setprecision(6) << pvalue << std::endl;
// Now loop over the edges
      for (unsigned int i=0; i<NEMAX; ++i){
          if( (svec[ibest].ebits).test(i) ){
              if(_printing>4)std:: cout << "Edge " << std::setw(2) << i << " " << std::setw(2) << evec[i].u << " " 
                                        << std::setw(2) << evec[i].v << " " << evec[i].pdgid << " " 
                                        << std::fixed << std::setw(10) << std::setprecision(6) << evec[i].edge_pvalue << std::endl;
          }
      }
   }

   gettimeofday(&tv,NULL);
   double t3 = tv.tv_sec+(tv.tv_usec/1000000.0);
   if(_printing>2)std::cout << "Solution analysis time = " << std::fixed << std::setw(12) << std::setprecision(6) << t3-t2 << " (s) " << std::endl;
   if(_printing>2)std::cout << "Total elapsed time     = " << std::fixed << std::setw(12) << std::setprecision(6) << t3-t0 << " (s) " << std::endl;

  // Algorithm 1. Based on fit probability ordering, add GammaGammaCandidates to the output GammaGammaParticles collection 
  // if they do not contain already assigned constituent particles. For now, the constituents are always two photons  
  // - but may not hurt to keep this generic for when other decay modes (eg. Dalitz decay) are included


  // Algorithm 2. Maximal solution containing pi0s and/or etas (but no etaprimes) with lowest chi-squared.

  if(_printing>3)std::cout << "ALGORITHM2 RESULTS" << std::endl;
  for(unsigned int i = 0; i < _pfovec.size(); ++i){
      if( (svec[ibest].ebits).test(i) ){
         const ReconstructedParticleVec particles = _pfovec[i]->getParticles();
         bool assignable = true;   
         for(unsigned int j = 0; j < particles.size(); ++j){           // TODO - maybe double check these are indeed photons ..
             const ReconstructedParticle *particle = particles[j];
             if( particles_assigned.find(particle) != particles_assigned.end() )assignable = false; // If particle is already assigned then it is not still assignable
         }
         if(assignable){
            recparcol->addElement(_pfovec[i]);                    // Add meson to the output collection if all constituent particles are still assignable
            if(_printing>3)std::cout << std::setw(4) << i;
            for(unsigned int j = 0; j < particles.size(); ++j){
                const ReconstructedParticle *particle = particles[j];
                particles_assigned.insert(particle);              // Update set of particles already assigned in this solution
// also keep track of photon indices as above
                retcode = daughter_particles.insert(particle);    // insert should fail as all photons should already be in the set
                if(_printing>3)std::cout << "  " << std::setw(4) << std::distance(daughter_particles.begin(),retcode.first) ;
            }
            const int particleType  = _pfovec[i]->getType();
         // GoodnessOfPID currently filled with fit probability - may at some point be more convenient to fill with the fit chi-squared ...
            const float GoodnessOfPID = _pfovec[i]->getGoodnessOfPID();
            if(_printing>3)std::cout << "  " << std::setw(4) << particleType << "  " 
                                     << std::fixed << std::setw( 10 ) << std::setprecision( 6 ) << GoodnessOfPID << std::endl;         
         }
      }         
  }
  return;
}

bool GammaGammaSolutionFinder::PfoProbabilitySortFunction(ReconstructedParticle* lhs,ReconstructedParticle* rhs){

  // Sort gamma gamma candidates by type and by decreasing fit probability

  //  true if lhs goes before

/*
   const double lhs_energy  = lhs->getEnergy();
   const double rhs_energy  = rhs->getEnergy();
*/
   const int lhs_particleType  = lhs->getType();
   const int rhs_particleType  = rhs->getType();

  // GoodnessOfPID currently filled with fit probability - may at some point be more convenient to fill with the fit chi-squared ...
   const float lhs_GoodnessOfPID = lhs->getGoodnessOfPID();
   const float rhs_GoodnessOfPID = rhs->getGoodnessOfPID();

   if(lhs_particleType==rhs_particleType)return (lhs_GoodnessOfPID>rhs_GoodnessOfPID);  // This rank makes sense when GoodnessOfPID has the fit probability
   return (lhs_particleType<rhs_particleType);                                          // Favor lower valued types (pi0s cf etas cf eta's)

}

