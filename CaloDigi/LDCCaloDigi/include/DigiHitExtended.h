#ifndef DigiHitExtended_H
#define DigiHitExtended_H 1

#include "lcio.h"
#include "EVENT/LCIO.h"
#include "EVENT/SimCalorimeterHit.h"

using namespace lcio;

/**
\addtogroup CaloDigi CaloDigi
@{

*/

typedef std::vector<SimCalorimeterHit * > SimCalorimeterHitVec;

class DigiHitExtended {

    public:

    DigiHitExtended();
    ~DigiHitExtended();

    void setAmpl(float ampl);
    void setCellID(int cellid);
    void setPosition(float * pos);
    void addSimHit(SimCalorimeterHit * simhit);

    float getAmpl();
    int getCellID();
    float * getPosition();
    SimCalorimeterHitVec & getSimHitVector();


    private:

    SimCalorimeterHitVec _simHitVec{};
    
    int _cellid{};
    float _pos[3]{};
    float _ampl{};

};

/** @}*/

#endif
