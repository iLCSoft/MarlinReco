#ifndef HISTOGRAM_
#define HISTOGRAM_ 1

#include "VObject.h"
#include <string>

class Histogram{
private:
  std::string histName{};     // name of the histogram
  int dim{};                  // dimension of the histogram = number of contained variables
  std::string *varNames{};    // array of names of the contained variables [dim] Order!!
  double *start{};            // array with start values of variables      [dim] Order!!
  double *width{};            // array with bin width of variables         [dim] Order!!
  int *bins{};                // array with number of bins per variable    [dim] Order!!
  int totNoOfBins{};          // total number of bins = bins[0]*...bins[dim-1]

  double *content{};          // array with histogram content [totNoOfBins]

  double density[3]{};        // [0] = maximum d., [1] = minimum d., [2] = average d. of histogram
  double norm{};              // norm of the histogramm (number of entries)
  double vol{};               // volume of one bin  // vol*norm = correct norm

public:
  Histogram(std::string HistName, int Dimension);          // creates arrays
  Histogram(std::string HistName, int Dimension, std::string *VarName, double *StartVal, double *BinWidth, int *NoOfBins);          // creates arrays
  ~Histogram();                                            // deletes arrays
  
  // set methods return 0 on succes and -1 on error (+ error message)
  // set variable VarNo's name, start value, bin width and number of bins
  int SetVariable(int VarNo, std::string VarName, double StartVal, double BinWidth, int NoOfBins);
  // set variable's  name, start value, bin width and number of bins all at once (not recommended to use)
  int SetVariables(std::string *VarName, double *StartVal, double *BinWidth, int *NoOfBins);
  
  int Fill(VObject *VO) ;     // fills histgramm with variable values in VO, returns 0 on success , -1 on error + message
  int SetBinContent(int BinInd, double Value) ; // fills bin with value, 0 on success, -1 on error + message
 
  double GetNorm() ;          // returns the Norm
  double GetVolume() ;        // returns the Bin volume (vol)
  std::string GetHistName() ; // returns the histogram name

  double *GetDensities();              // returns the density array
  double GetContent(VObject *VO);      // returns content of bin where variables in VO lie in, or 0 on error + message
  double GetNormContent(VObject *VO);  // returns to 'norm' normalized content of bin where VO is in, or 0 on error
  double GetBinContent(int BinInd);    // returns content of bin 'BinInd', or 0 on error + message
  
  int GetDimension();                  // returns dim
  int GetTotNoOfBins();                // returns totNoOfBins
  std::string *GetVarNames();          // returns varNames
  double *GetVarStart();               // returns start
  double *GetVarWidth();               // returns width
  int *GetVarBins();                   // returns bins

  double CompareHist(Histogram *hist); // returns mean squared distance of bins divided by number totNoOfBins
                                       // of -1 on error and a message
  
};

#endif
