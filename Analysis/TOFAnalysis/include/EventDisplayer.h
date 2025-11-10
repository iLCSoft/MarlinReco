#ifndef EventDisplayer_h
#define EventDisplayer_h 1

#include "DD4hep/Detector.h"
#include "marlin/Processor.h"
#include "EVENT/LCEvent.h"
#include "marlinutil/DDMarlinCED.h"

/**
 * Mixin class for attaching an event display to any inheriting processor.
 * The inheriting processor has to make this class a friend as well, and make
 * sure to call the EventDisplayer constructor passing itself for registering
 * an additional processor parameter, steering whether to actually run an event
 * display or not.
 */
class EventDisplayer {
    public:
        template<typename T>
        EventDisplayer(T* proc);

        /// Initialize the display. Call this in init of your processor
        void initDisplay(marlin::Processor* proc);

        /**
         * Draw something to the event display. Where something is determined by
         * by what happens inside Func. Func is an arbitrary function taking
         * an arbitrary number of arguments that also have to be passed along
         * to this function.
         */
        template<typename Func, typename... Args>
        void drawDisplay(marlin::Processor* proc, EVENT::LCEvent* evt, Func&& func, Args&&... args);

    private:
        dd4hep::Detector& _detector = dd4hep::Detector::getInstance();
        bool _eventDisplay{false};
        bool _isInit{false};

};

void EventDisplayer::initDisplay(marlin::Processor* proc){
    if (_eventDisplay && !_isInit) {
        DDMarlinCED::init(proc);
        _isInit = true;
    }
}

template<typename Func, typename... Args>
void EventDisplayer::drawDisplay(marlin::Processor* proc, EVENT::LCEvent* evt, Func&& func, Args&&... args) {
    if (!_eventDisplay) {
        return;
    }
    if (!_isInit) {
        initDisplay(proc);
    }

    DDMarlinCED::newEvent(proc, evt);
    DDMarlinCED::drawDD4hepDetector(_detector, false, std::vector<std::string>{""});
    DDCEDPickingHandler& pHandler= DDCEDPickingHandler::getInstance();
    pHandler.update(evt);

    func(std::forward<Args>(args)...);

    DDMarlinCED::draw(proc);
}


template<typename T>
EventDisplayer::EventDisplayer(T* proc) {
    proc->registerProcessorParameter("eventDisplay",
                                     "eventDisplay",
                                    _eventDisplay,
                                    bool(false));
}


#endif
