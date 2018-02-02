#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"

#include <string>

void draw_cluster_time(const std::string& histName){

  std::string treeName("MyTOFPlots/") ;
  treeName += histName ;


  TFile fka("aida_file_kaon_random_REC.root");
  TH1* tka = (TH1*) fka.Get( treeName.c_str() ) ;

  TFile fpi("aida_file_pion_random_REC.root") ;
  TH1* tpi = (TH1*) fpi.Get( treeName.c_str() ) ;

  TFile fmu("aida_file_muon_random_REC.root") ;
  TH1* tmu = (TH1*) fmu.Get( treeName.c_str() ) ;

  TFile fpr("aida_file_proton_random_REC.root") ;
  TH1* tpr = (TH1*) fpr.Get( treeName.c_str() ) ;

  TFile fel("aida_file_electron_random_REC.root") ;
  TH1* tel = (TH1*) fel.Get( treeName.c_str() ) ;


 
  tka->SetLineColor( kRed ) ;
  tpi->SetLineColor( kGreen ) ;
  tmu->SetLineColor( kBlue ) ;
  tpr->SetLineColor( kYellow+1 ) ;
  tel->SetLineColor( kOrange ) ;

  tka->SetMarkerColor( kRed ) ;
  tpi->SetMarkerColor( kGreen ) ;
  tmu->SetMarkerColor( kBlue ) ;
  tpr->SetMarkerColor( kYellow+1 ) ;
  tel->SetMarkerColor( kOrange ) ;


  TLegend leg(.1,.75,.2,.9,"");
  leg.SetFillColor(0);
  tmu->SetFillColor(0);
  tpi->SetFillColor(0);
  tka->SetFillColor(0);
  tpr->SetFillColor(0);
  tel->SetFillColor(0);

  tmu->Draw() ;
  tka->Draw("same") ;
  tpi->Draw("same") ;
  tel->Draw("same") ;
  tpr->Draw("same") ;

  TLegendEntry * l = leg.AddEntry(&fka,"kaons","l") ;
  l->SetLineColor( kRed ) ;
  l->SetMarkerColor( kRed ) ;

  l = leg.AddEntry(&fel,"electrons","l");
  l->SetLineColor( kOrange );
  l->SetMarkerColor( kOrange );

  l = leg.AddEntry(&fpi,"pions","l") ;
  l->SetLineColor( kGreen ) ;
  l->SetMarkerColor( kGreen ) ;

  l = leg.AddEntry(&fmu,"muons","l") ;
  l->SetLineColor( kBlue ) ;
  l->SetMarkerColor( kBlue ) ;

  l = leg.AddEntry(&fpr,"protons","l") ;
  l->SetLineColor( kYellow+1 ) ;
  l->SetMarkerColor( kYellow+1 ) ;

  // l = leg.AddEntry(&fph,"photons","l") ;
  // l->SetLineColor( kViolet ) ;

  leg.DrawClone("Same");

  std::string filename( histName ) ;
  filename += ".png" ;

  gPad->Print( filename.c_str() ) ;

  //  gPad->WaitPrimitive() ;

//  fka.Close() ;
//  fpi.Close() ;
//  fmu.Close() ;
//  fpr.Close() ;
//  fel.Close() ;
//  fph.Close() ;
}
