#ifndef VECTORHELPER_H
#define VECTORHELPER_H 1
#include "TJjetsPFOAnalysisProcessor.h"

template <class Object> bool TJjetsPFOAnalysisProcessor::areDisjointVectors( const std::vector<Object> &v1, const std::vector<Object> &v2 ) const {
  bool have_overlap = false;
  for (unsigned int i=0; i<v1.size(); i++) {
    Object v1_element = v1[i];
    if(std::find(v2.begin(), v2.end(), v1_element) != v2.end()) {
      have_overlap = true;
      break;
    }
  }
  return !have_overlap;
}
#endif
