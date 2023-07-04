#ifndef ModelMgr_h
#define ModelMgr_h 1

#include <TrainingModel.h>

#include <map>
#include <set>

namespace cpid
{
  class ModelMgr
  {
    friend class TrainingModel;

    public:

      static ModelMgr* instance();

      virtual ~ModelMgr(){};

      // create new algorithm to be used in processor
      // the caller takes ownership here - don't forget to delete!
      TrainingModel* createModel(const std::string& type);

      // get archetype algorithm from map
      TrainingModel* getModel(const std::string& type);

      void printAvailableModelTypes();

    protected:

      void registerModel(TrainingModel* model);

      std::vector<std::string> getAvailableModelTypes();

      ModelMgr();

    private:

      static ModelMgr* _me;
      std::map<const std::string, TrainingModel*> _map{};
  };
}

#endif
