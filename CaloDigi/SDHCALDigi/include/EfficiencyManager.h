#ifndef EfficiencyManager_h
#define EfficiencyManager_h

#include <map>
#include <string>

struct AsicKey;
class SimDigitalGeomCellId;

class EfficiencyManager {
public:
  EfficiencyManager();
  virtual ~EfficiencyManager();

  virtual float getEfficiency(SimDigitalGeomCellId* cellID) = 0;
};

class UniformEfficiency : public EfficiencyManager {
public:
  UniformEfficiency(float val = 1.0);
  virtual ~UniformEfficiency();

  virtual float getEfficiency(SimDigitalGeomCellId* cellID);

protected:
  float value;
};

class AsicEfficiency : public UniformEfficiency {
public:
  AsicEfficiency(std::string fileName, float globalVal = 1.0);
  virtual ~AsicEfficiency();

  virtual float getEfficiency(SimDigitalGeomCellId* cellID);

protected:
  void readFile(std::string fileName);

  std::map<AsicKey, float> effMap;
};

#endif // EfficiencyManager_h
