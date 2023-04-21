#ifndef AlgorithmMgr_h
#define AlgorithmMgr_h 1

#include <InputAlgorithm.h>

#include <map>
#include <set>

namespace cpid
{
  class AlgorithmMgr
  {
    friend class InputAlgorithm;

    public:

      static AlgorithmMgr* instance();

      virtual ~AlgorithmMgr(){};

      // create new algorithm to be used in processor
      InputAlgorithm* createAlgorithm(const std::string& type);

      // get archetype algorithm from map
      InputAlgorithm* getAlgorithm(const std::string& type);

      void printAvailableAlgorithmTypes();

    protected:

      void registerAlgorithm(InputAlgorithm* algorithm);

      std::set< std::string > getAvailableAlgorithmTypes();

      AlgorithmMgr();

    private:

      static AlgorithmMgr* _me;
      std::map<const std::string, InputAlgorithm*> _map{};
      std::list<InputAlgorithm*> _list{};
  };
}

#endif
