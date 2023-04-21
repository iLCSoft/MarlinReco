#include <AlgorithmMgr.h>

namespace cpid
{
  AlgorithmMgr* AlgorithmMgr::_me = 0;

  AlgorithmMgr::AlgorithmMgr(){}

  AlgorithmMgr* AlgorithmMgr::instance()
  {
    if(_me == 0) _me = new AlgorithmMgr;
    return _me;
  }

  void AlgorithmMgr::registerAlgorithm(InputAlgorithm* algorithm)
  {
    const std::string& name = algorithm->type();

        if(_map.find(name) != _map.end()) return;
        else _map[name] = algorithm;
  }

  std::set<std::string> AlgorithmMgr::getAvailableAlgorithmTypes()
  {
    std::set<std::string> atypes;
    for(std::map<const std::string, InputAlgorithm*>::iterator i=_map.begin() ; i!= _map.end() ; i++) atypes.insert(i->first);
    return atypes;
  }

  InputAlgorithm* AlgorithmMgr::createAlgorithm(const std::string& type)
  {
    InputAlgorithm* algorithm = getAlgorithm(type);
    InputAlgorithm* newAlgorithm = algorithm->newAlgorithm();
    return newAlgorithm;
  }

  InputAlgorithm* AlgorithmMgr::getAlgorithm(const std::string& type)
  {
    return _map[type];
  }

  void AlgorithmMgr::printAvailableAlgorithmTypes()
  {
    std::cout << "available algorithms: " << std::endl;
    for (const auto& [key, value] : _map) {std::cout << key << std::endl;}
  }

}
