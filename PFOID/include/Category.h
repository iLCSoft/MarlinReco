#ifndef CATEGORY_
#define CATEGORY_ 1

#include "Histogram.h"
#include <string>
#include <vector>

typedef std::vector<Histogram*> HistVec;

class Category{            // Container for Histograms
private:
  std::string catName{};    // name of the category, e.g. electrons, muons, pions
  int nHists{};             // number of histograms

public:
  HistVec hists{};          // histograms (public to be access from outside)

  // constructor: creates arrays
  Category(std::string Name, int NoOfHist);
  ~Category();                              // destructor: deletes arrays

  std::string GetName();   // returns category name 'catName'
  unsigned int GetNoOfHists();      // returns number of histograms contained in this category

};


#endif
