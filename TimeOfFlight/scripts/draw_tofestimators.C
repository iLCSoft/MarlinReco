#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"

#include <string>

/*

.x ../TimeOfFlight/scripts/draw_tofestimators.C("./tof_estimators.root", "TOFEstimators50ps" , "hbetaFirstHitsChrg" )
.x ../TimeOfFlight/scripts/draw_tofestimators.C("./tof_estimators.root", "TOFEstimators50ps" , "hbetaCloseHitsChrg" )
.x ../TimeOfFlight/scripts/draw_tofestimators.C("./tof_estimators.root", "TOFEstimators50ps" , "hbetaFirstHitsNeut" )
.x ../TimeOfFlight/scripts/draw_tofestimators.C("./tof_estimators.root", "TOFEstimators50ps" , "hbetaCloseHitsNeut" )
.x ../TimeOfFlight/scripts/draw_tofestimators.C("./tof_estimators.root", "TOFEstimators50ps" , "hbetaLastTrkHit" )

 */


void draw_tofestimators(const std::string& fileName,
			std::string treeName,
			const std::string& histName){

  std::string outFile = treeName + "_" + histName + ".png" ;


  treeName += "/" ;
  treeName += histName ;

  

  TFile file( fileName.c_str() );
  
  TH2*  h = (TH2*) file.Get( treeName.c_str() ) ;

  std::cout << "  got histo   " << treeName << " : h = " << h << std::endl ;

  if( h != 0 ) 
    h->Draw() ;

 
  TF1 *fe  = new TF1("fe","x/sqrt(x*x+0.0005109989461*0.0005109989461)",0.1,10.);
  TF1 *fmu = new TF1("fmu","x/sqrt(x*x+0.105658*0.105658)",0.1,10.);
  TF1 *fpi = new TF1("fpi","x/sqrt(x*x+0.13957018*0.13957018)",0.1,10.);
  TF1 *fk  = new TF1("fk","x/sqrt(x*x+0.493677*0.493677)",0.1,10.);
  TF1 *fp  = new TF1("fp","x/sqrt(x*x+0.9382720813*0.9382720813)",0.1,10.);

  fe->SetLineColor( kBlack ) ;
  fmu->SetLineColor( kBlack ) ;
  fpi->SetLineColor( kBlack ) ;
  fk->SetLineColor( kBlack ) ;
  fp->SetLineColor( kBlack ) ;

  fe-> SetLineWidth( 1 ) ;
  fmu->SetLineWidth( 1 ) ;
  fpi->SetLineWidth( 1 ) ;
  fk-> SetLineWidth( 1 ) ;
  fp-> SetLineWidth( 1 ) ;

  fe->Draw("same") ;
  fmu->Draw("same") ;
  fpi->Draw("same") ;
  fk->Draw("same") ;
  fp->Draw("same") ;

  gPad->Print( outFile.c_str() ) ;

}
