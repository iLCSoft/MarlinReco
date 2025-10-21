#ifndef TrackZVertexGrouping_h
#define TrackZVertexGrouping_h 1

#include "lcio.h"
#include "marlin/Processor.h"
#include <gsl/gsl_rng.h>
#include <string>
#include <vector>

using namespace lcio;
using namespace marlin;

class TH1F;

/** Group Tracks into clusters with consistent z-positions of their vertex, based on the Z0 significance.
 *  Algorithm developed by S.Sasikumar, DESY.
 *
 *
 * @author F.Gaede, DESY, September 2018
 */

class TrackZVertexGrouping : public Processor {

public:
  virtual Processor* newProcessor() { return new TrackZVertexGrouping; }

  TrackZVertexGrouping(const TrackZVertexGrouping&) = delete;
  TrackZVertexGrouping& operator=(const TrackZVertexGrouping&) = delete;

  TrackZVertexGrouping();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

protected:
  /** Input collection name with Tracks
   */
  std::string _colNameTracks{};
  std::string _colNameTrkGroupPFOs{};
  std::string _colNameTrkGroupVertices{};

  float _z0SignificanceCut{};

  int _nRun{};
  int _nEvt{};

  std::vector<TH1F*> _h{};
};

#endif
