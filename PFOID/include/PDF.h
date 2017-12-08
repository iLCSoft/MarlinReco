#ifndef PDFclass_
#define PDFclass_ 1

#include "Category.h"
#include "Histogram.h"
#include "VObject.h"
#include <string>

typedef std::vector<Category*> CatVec;

class PDF{
private:
  std::string pdfName{};         // name of the PDF object, e.g. charged particles
  int nCats{};                   // number of categories
  CatVec cat{};                   // array of categories in the pdf
  unsigned int IniCount{};

public:
  VObject *VO{};
  
  // constructor: creates arrays
  PDF(std::string PDFname, int NoOfCats, std::string* CatNames, int NoOfHists, int NoOfVar, std::string *VarNames);
  PDF(std::string Filename);                // overloaded constructor: creates arrays, fills from file
  ~PDF();                                   // destructor:  deletes arrays

  // for filling the histograms in the categories manually
  int InitHistogram(std::string HistName, int Dim, std::string *VarName, double *StartVal, double *BinWidth, int *NoOfBins);
  int FillHistograms(std::string CatName);
  Histogram * GetHistogram(std::string CatName, std::string HistName);

  // returns likelihood of belonging to category 'CatName'
  double GetLikelihood(std::string CatName);
  
  int ReadPDF(std::string Filename);        // reads a pdf object from file
  int WritePDF(std::string Filename);       // writes a pdf object to file
  int GetNCats();                           // returns nCats
  std::string GetCatName(int CatNo);        // returns name of cat 'CatNo'

};


#endif
