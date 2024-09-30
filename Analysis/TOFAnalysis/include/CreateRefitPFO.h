#ifndef CreateRefitPFO_h
#define CreateRefitPFO_h 1

#include <vector>
#include "marlin/Processor.h"
#include "UTIL/LCRelationNavigator.h"
#include "EVENT/Track.h"
#include "TLorentzVector.h"

class CreateRefitPFO : public marlin::Processor{
	public:
		marlin::Processor*  newProcessor(){ return new CreateRefitPFO; }
		CreateRefitPFO();
		~CreateRefitPFO() = default;

		void init();
		void processEvent( EVENT::LCEvent *event );
        void end(){};

        int getTrackPDG(EVENT::Track* track, UTIL::LCRelationNavigator& nav);
		int getTrackIndex(EVENT::LCCollection* trackCollection, EVENT::Track* selectedTrack);
		TLorentzVector getTrackFourMomentum(EVENT::Track* track , double mass);
		std::vector<float> updateChargedPFOCovMat(EVENT::Track* track , double mass);
    private:
        int _nEvt;
        double _bField;        
};

#endif
