#include <ModelMgr.h>

namespace cpid
{
  ModelMgr* ModelMgr::_me = 0;

  ModelMgr::ModelMgr(){}

  ModelMgr* ModelMgr::instance()
  {
    if(_me == 0) _me = new ModelMgr;
    return _me;
  }

  void ModelMgr::registerModel(TrainingModel* model)
  {
    const std::string& name = model->type();

        if(_map.find(name) != _map.end()) return;
        else _map[name] = model;
  }

  std::set<std::string> ModelMgr::getAvailableModelTypes()
  {
    std::set<std::string> mtypes;
    for(std::map<const std::string, TrainingModel*>::iterator i=_map.begin() ; i!= _map.end() ; i++) mtypes.insert(i->first);
    return mtypes;
  }

  TrainingModel* ModelMgr::createModel(const std::string& type)
  {
    TrainingModel* model = getModel(type);
    TrainingModel* newModel = model->newModel();
    return newModel;
  }

  TrainingModel* ModelMgr::getModel(const std::string& type)
  {
    return _map[type];
  }

  void ModelMgr::printAvailableModelTypes()
  {
    std::cout << "Available training models: " << std::endl;
    for (const auto& [key, value] : _map) {sloM << " " << key << std::endl;}
    sloM << "-------------------------" << std::endl;
  }

}
