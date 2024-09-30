#include "BohdanDrawing.h"
#include "TOF.h"
#include "EVENT/TrackState.h"
#include "UTIL/Operators.h"
#include "UTIL/TrackTools.h"
#include "marlinutil/DDMarlinCED.h"
#include "marlinutil/CalorimeterHitType.h"

#include "TF1.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TColor.h"
#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <numeric>
#include <cstdlib>
	
using CLHEP::RandGauss;
using dd4hep::rec::Vector3D;

static int picture_counter = 1;

TStyle* getMyStyle(){
    //Mostly CMS style with some ILD style comments
    TStyle* myStyle = new TStyle("myStyle", "My Style");

    // For the canvas:
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetCanvasColor(kWhite);  // ILD 10
    myStyle->SetCanvasDefH(600); //Height of canvas
    myStyle->SetCanvasDefW(600); //Width of canvas
    myStyle->SetCanvasDefX(0);   //Position on screen
    myStyle->SetCanvasDefY(0);

    // For the Pad:
    myStyle->SetPadBorderMode(0);
    myStyle->SetPadGridX(false);
    myStyle->SetPadGridY(false);
    myStyle->SetPadColor(kWhite); // ILD 10
    myStyle->SetGridColor(0);
    myStyle->SetGridStyle(3);
    myStyle->SetGridWidth(1);

    // For the frame:
    myStyle->SetFrameBorderMode(0);
    myStyle->SetFrameFillColor(0); // ILD 10
    myStyle->SetFrameLineWidth(1); // ILD 2
    myStyle->SetFrameBorderSize(1);
    myStyle->SetFrameFillStyle(0);
    myStyle->SetFrameLineColor(1);
    myStyle->SetFrameLineStyle(1);

    // For the histo:
    // myStyle->SetHistFillColor(1);
    // myStyle->SetHistFillStyle(0);
    myStyle->SetHistLineColor(1);
    myStyle->SetHistLineStyle(0);
    myStyle->SetHistLineWidth(1); // ILD 2
    // myStyle->SetLegoInnerR(Float_t rad = 0.5);
    myStyle->SetNumberContours(256); //default 20

    myStyle->SetEndErrorSize(2);
    myStyle->SetErrorX(0.);

    myStyle->SetMarkerStyle(20);

    //For the fit/function:
    myStyle->SetOptFit(1); // ILD 0
    myStyle->SetFitFormat("5.4g");
    myStyle->SetFuncColor(2);
    myStyle->SetFuncStyle(1);
    myStyle->SetFuncWidth(1); // ILD 2

    //For the date:
    myStyle->SetOptDate(0);
    // myStyle->SetDateX(Float_t x = 0.01);
    // myStyle->SetDateY(Float_t y = 0.01);

    // For the statistics box:
    myStyle->SetOptFile(0);
    myStyle->SetOptStat(0); // mean and RMS: SetOptStat("mr");
    myStyle->SetStatColor(kWhite); // ILD 10
    myStyle->SetStatFont(42);
    myStyle->SetStatFontSize(0.025); // ILD 0.07
    myStyle->SetStatTextColor(1);
    myStyle->SetStatFormat("6.4g");
    myStyle->SetStatBorderSize(1); // ILD 0
    myStyle->SetStatH(0.1);
    myStyle->SetStatW(0.15);
    // myStyle->SetStatStyle(Style_t style = 1001);
    // myStyle->SetStatX(Float_t x = 0);
    // myStyle->SetStatY(Float_t y = 0);

    // Margins:
    auto margin = 0.18;
    auto leftMarginFraction = 0.89;
    auto bottomMarginFraction = 0.72;
    myStyle->SetPadTopMargin((1-bottomMarginFraction)*margin);
    myStyle->SetPadBottomMargin(bottomMarginFraction*margin);
    myStyle->SetPadLeftMargin(leftMarginFraction*margin);
    myStyle->SetPadRightMargin((1-leftMarginFraction)*margin);

    // For the Global title:
    myStyle->SetOptTitle(0);    // 0=No Title
    myStyle->SetTitleFont(42);
    myStyle->SetTitleColor(1);
    myStyle->SetTitleTextColor(1);
    myStyle->SetTitleFillColor(10); // ILD 0
    myStyle->SetTitleFontSize(0.05);
    // myStyle->SetTitleH(0); // Set the height of the title box
    // myStyle->SetTitleW(0); // Set the width of the title box
    // myStyle->SetTitleX(0); // Set the position of the title box
    // myStyle->SetTitleY(0.985); // Set the position of the title box
    // myStyle->SetTitleStyle(Style_t style = 1001);
    // myStyle->SetTitleBorderSize(2); // ILD 0

    // For the axis titles:
    myStyle->SetTitleColor(1, "XYZ");
    myStyle->SetTitleFont(42, "XYZ");
    myStyle->SetTitleSize(0.06, "XYZ"); // ILD 0.07
    // myStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
    // myStyle->SetTitleYSize(Float_t size = 0.02);
    myStyle->SetTitleXOffset(0.9); // ILD 1
    myStyle->SetTitleYOffset(1.25); // ILD 1.1
    // myStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

    // For the axis labels:
    myStyle->SetLabelColor(1, "XYZ");
    myStyle->SetLabelFont(42, "XYZ");
    myStyle->SetLabelOffset(0.007, "XYZ"); // ILD 0.015
    myStyle->SetLabelSize(0.05, "XYZ"); // ILD 0.06

    // For the axis:
    myStyle->SetAxisColor(1, "XYZ");
    myStyle->SetStripDecimals(kTRUE);
    myStyle->SetTickLength(0.03, "XYZ");
    myStyle->SetNdivisions(510, "XYZ"); // ILD 506
    myStyle->SetPadTickX(0);  // 0=Text labels (and tics) only on bottom, 1=Text labels on top and bottom. ILD 1
    myStyle->SetPadTickY(1);

    // Change for log plots:
    myStyle->SetOptLogx(0);
    myStyle->SetOptLogy(0);
    myStyle->SetOptLogz(0);

    // Postscript options:
    myStyle->SetPaperSize(20.,20.);
    // myStyle->SetLineScalePS(Float_t scale = 3);
    // myStyle->SetLineStyleString(Int_t i, const char* text);
    // myStyle->SetHeaderPS(const char* header);
    // myStyle->SetTitlePS(const char* pstitle);

    // myStyle->SetBarOffset(Float_t baroff = 0.5);
    // myStyle->SetBarWidth(Float_t barwidth = 0.5);
    // myStyle->SetPaintTextFormat(const char* format = "g");
    // myStyle->SetTimeOffset(float_t toffset);
    // myStyle->SetHistMinimumZero(kTRUE);

    // myStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    myStyle->SetPalette(1,0); // from ILD

    return myStyle;
}

void displayPFO(EVENT::ReconstructedParticle* pfo){
    if (pfo == nullptr) return;
    std::vector<EVENT::Track*> tracks;
    if ( not pfo->getTracks().empty() )
        tracks = getSubTracks(pfo->getTracks()[0]);

    //plot settings for tracker hits
    int type = 0; // point
    int size = 4; // larger point
    unsigned long color = 0x000000;
    for(auto* track: tracks){
        auto hits = track->getTrackerHits();
        for (auto* hit : hits ){
            auto pos = hit->getPosition();
    
            ced_hit(pos[0], pos[1], pos[2], type, size, color);
            // hitIdx++;
        }
    }
    //plot ecal state for fun
    if (not tracks.empty()){
        const TrackState* tsCalo = tracks.back()->getTrackState( TrackState::AtCalorimeter );
        auto posCalo = tsCalo->getReferencePoint();
        ced_hit(posCalo[0], posCalo[1], posCalo[2] + tsCalo->getZ0(), 0, 4, 0x000000); // track state at calorimeterer
    }


    std::vector<EVENT::Cluster*> clusters = pfo->getClusters();
    type = 0; // point
    size = 4; // larger point
    color = 0x000000;
    for(auto* cluster: clusters){
        auto hits = cluster->getCalorimeterHits();
        for (auto* hit: hits){
            auto pos = hit->getPosition();
            ced_hit(pos[0], pos[1], pos[2], type, size, color);
        }
    }

}


void displayFTDSimHits(EVENT::LCEvent* evt){
    EVENT::LCCollection* hits = evt->getCollection("FTDCollection");

    for (int i=0; i < hits->getNumberOfElements(); i++) {
        SimTrackerHit* hit = static_cast<EVENT::SimTrackerHit*>(hits->getElementAt(i));
        auto pos = hit->getPosition();
        int type = 0; // point
        int size = 4; // larger point
        unsigned long color = 0x000000;
        ced_hit(pos[0], pos[1], pos[2], type, size, color);
    }
}

void plotECALTimes(EVENT::Cluster* cluster, Vector3D posAtEcal, Vector3D momAtEcal, MCParticle* mc){
    std::cout<<"PDG: "<<mc->getPDG()<<std::endl;
    Vector3D mom(mc->getMomentum());
    std::cout<<"Momentum: "<<mom.r()<<" ( "<<mom.rho()<<", "<<mom.z()<<")"<<std::endl;

    auto hits = cluster->getCalorimeterHits();
    auto frankHits = selectFrankEcalHits(cluster, posAtEcal, momAtEcal, 10);

    //fill maps
    std::map <CalorimeterHit*, float> hit2dToTrack;
    std::map <CalorimeterHit*, float> hit2time;
    std::map <CalorimeterHit*, float> hit2timeSmeared;
    for(auto hit : hits){
        bool isECALHit = ( CHT( hit->getType() ).caloID() == CHT::ecal );
        if ( (!isECALHit) ) continue;
        hit2dToTrack[hit] = (Vector3D(hit->getPosition()) - posAtEcal).r();
        hit2time[hit] = hit->getTime();
        hit2timeSmeared[hit] = RandGauss::shoot(hit->getTime(), 0.1);
    }
    if ( hit2time.empty() ) return;

    TStyle* myStyle = getMyStyle();
    myStyle->cd();
    gROOT->ForceStyle();
    TCanvas c = TCanvas("c", "All ECAL hits");

    std::vector<float> x_all, y_all, y_all_smeared;
    for(auto hit : hits){
        bool isECALHit = ( CHT( hit->getType() ).caloID() == CHT::ecal );
        if ( (!isECALHit) ) continue;
        x_all.push_back( hit2dToTrack[hit] );
        y_all.push_back( hit2time[hit] );
        y_all_smeared.push_back( hit2timeSmeared[hit] );
    }
    std::vector<float> y_err(x_all.size(), 0.1);
    TGraphErrors gr( x_all.size(), x_all.data(), y_all.data(), nullptr, y_err.data() );
    gr.GetXaxis()->SetTitle("d to track (mm)");
    gr.GetYaxis()->SetTitle("Hit time (ns)");
    gr.Draw("AP");
    TLatex t;
    t.SetNDC();
    t.DrawLatex(0.2,  0.87, Form("#splitline{PDG: %d}{p: %.1f GeV}", mc->getPDG(), Vector3D(mc->getMomentum()).r() ));
    TLegend* legend = new TLegend(0.57, 0.84, 0.91, 0.93);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(62);

    legend->AddEntry(&gr, "ECAL hits (true time)", "pe");
    legend->Draw();
    c.Modified(); c.Update();
    c.SaveAs( Form("./plots/steps/%d_step1.png", picture_counter) );
    // while ( !gSystem->ProcessEvents() ){
    //     gSystem->Sleep(5);
    // }

    gr.SetMarkerColor(15);
    gr.SetLineColor(15);
    TGraph gr2( x_all.size(), x_all.data(), y_all_smeared.data());
    gr2.SetTitle("100 ps smeared hit time;d to track (mm); Hit time (ns)");
    gr2.Draw("Psame");

    legend->Clear();
    legend->AddEntry(&gr, "Previous step", "pe");
    legend->AddEntry(&gr2, "Smear time 100 ps", "p");
    legend->Draw();

    c.Modified(); c.Update();
    c.SaveAs( Form("./plots/steps/%d_step2.png", picture_counter) );
    // while ( !gSystem->ProcessEvents() ){
    //     gSystem->Sleep(5);
    // }

    std::vector<float> x_frank, y_frank_smeared;
    for(auto hit : frankHits){
        bool isECALHit = ( CHT( hit->getType() ).caloID() == CHT::ecal );
        if ( (!isECALHit) ) continue;
        x_frank.push_back( hit2dToTrack[hit] );
        y_frank_smeared.push_back( hit2timeSmeared[hit] );
    }
 
    TGraph gr3( x_frank.size(), x_frank.data(), y_frank_smeared.data() );
    TF1 f("f", Form("%f + (x-%f)/%f", y_all_smeared.front(), x_all.front(), CLHEP::c_light), 0., *std::max_element(x_all.begin(), x_all.end()) );
    gr3.SetTitle("Selected hits;d to track (mm); Hit time (ns)");
    gr2.SetMarkerColor(15);
    gr2.SetLineColor(15);
    gr2.Draw("P");
    gr3.Draw("Psame");

    f.Draw("same");

    legend->Clear();
    legend->AddEntry(&gr2, "Previous step", "p");
    legend->AddEntry(&gr3, "Selected hits", "p");
    legend->AddEntry(&f, "Speed of light", "l");
    legend->Draw();

    c.Modified(); c.Update();
    c.SaveAs( Form("./plots/steps/%d_step3.png", picture_counter) );
    // while ( !gSystem->ProcessEvents() ){
    //     gSystem->Sleep(5);
    // }



    std::vector<float> y_frank_smeared_corr;
    for(size_t i=0; i < x_frank.size(); i++){
        y_frank_smeared_corr.push_back( y_frank_smeared[i] - x_frank[i]/CLHEP::c_light );
    }

    TGraph gr4( x_frank.size(), x_frank.data(), y_frank_smeared_corr.data() );
    float average = std::reduce(y_frank_smeared_corr.begin(), y_frank_smeared_corr.end()) / float(y_frank_smeared_corr.size());
    TF1 f2("f", Form("%f", average), 0., *std::max_element(x_all.begin(), x_all.end()) );

    gr4.SetTitle("Corrected hit time;d to track (mm); Corrected hit time (ns)");
    gr.Draw("AP");
    gr.SetLineStyle(0);
    gr.SetLineColor(0);
    gr.SetMarkerStyle(0);
    gr.SetLineColor(0);
    gr4.Draw("Psame");

    f2.Draw("same");

    legend->Clear();
    legend->AddEntry(&gr4, "Time corrected", "p");
    legend->AddEntry(&f2, Form("Average (TOF) = %.1f ns", average), "l");
    legend->Draw();

    c.Modified(); c.Update();
    c.SaveAs( Form("./plots/steps/%d_step4.png", picture_counter) );
    // while ( !gSystem->ProcessEvents() ){
    //     gSystem->Sleep(5);
    // }
    picture_counter++;

}

void plotTrackParams(const std::vector<HitState>& trackHitStates, EVENT::MCParticle* mc, float bField){
    // DEBUG output
    int pdg = std::abs( mc->getPDG() );
    std::cout<<"PDG: "<<pdg<<std::endl;

    Vector3D mcMom(mc->getMomentum());
    std::cout<<"Momentum: "<<mcMom.r()<<" ( "<<mcMom.rho()<<", "<<mcMom.z()<<")"<<std::endl;

    auto omega2pt = [bField](float omega){
        return 0.000299792458 * bField / omega;
    };
    auto pt2omega = [bField](float pt){
        return 0.000299792458 * bField / pt;
    };

    auto tanL2theta = [bField](float tanL){
        return 90. - 180.*std::atan(tanL)/M_PI;
    };

    TStyle* myStyle = getMyStyle();

    struct Margin{
        float left, right, top, bottom;
        Margin(float l, float r, float t, float b){
            left = l;
            right = r;
            top = t;
            bottom = b;
        }
    };
    Margin margin(0.175, 0.175, 0.08, 0.15);
    myStyle->SetPadTopMargin(margin.top); // ILD 0.08
    myStyle->SetPadBottomMargin(margin.bottom); // ILD 0.18
    myStyle->SetPadLeftMargin(margin.left); // ILD 0.17
    myStyle->SetPadRightMargin(margin.right); // ILD 0.08
    myStyle->SetPadTickY(0);

    myStyle->cd();
    gROOT->ForceStyle();
    TCanvas c = TCanvas("c", "Track parameters evolution through the track", 600/(1. - margin.left - margin.right), 600/(1. - margin.top - margin.bottom));
    c.SetGridx(true);
    std::vector<float> z, omega, tanL, trueZ, trueOmega, trueTanL;
    for(auto hitState : trackHitStates){
        auto simHit = hitState.simHit;
        auto ts = hitState.ts;
        if (simHit == nullptr) continue;
        z.push_back( std::abs(ts.getReferencePoint()[2] + ts.getZ0()) );
        omega.push_back( std::abs(ts.getOmega()) );
        tanL.push_back( std::abs(ts.getTanLambda()) );

        trueZ.push_back( std::abs(simHit->getPosition()[2]) );
        Vector3D mom(simHit->getMomentum());
        trueOmega.push_back( pt2omega(mom.rho()) );
        trueTanL.push_back( std::abs(mom.z())/mom.rho() );
    }
    TGraph grOmega( z.size(), z.data(), omega.data());
    TGraph grTrueOmega( trueZ.size(), trueZ.data(), trueOmega.data());
    TGraph grTanL( z.size(), z.data(), tanL.data());
    TGraph grTrueTanL( trueZ.size(), trueZ.data(), trueTanL.data());



    grOmega.SetLineColor(TColor::GetColor("#377eb8"));
    grOmega.SetMarkerColor(TColor::GetColor("#377eb8"));
    grOmega.SetLineWidth(3);
    grTrueOmega.SetLineColor(TColor::GetColor("#e41a1c"));
    grTrueOmega.SetMarkerColor(TColor::GetColor("#e41a1c"));
    grTrueOmega.SetLineWidth(3);
    TMultiGraph mg;
    mg.Add(&grOmega);
    mg.Add(&grTrueOmega);
    mg.GetXaxis()->SetNoExponent(true);
    mg.GetXaxis()->SetTitle("z (mm)");
    mg.GetXaxis()->SetNdivisions(508);
    mg.GetYaxis()->SetTitle("#Omega (1/mm)");
    mg.GetYaxis()->SetTitleOffset(1.5);
    mg.GetXaxis()->CenterTitle(true);
    mg.GetYaxis()->CenterTitle(true);
    mg.Draw("APL");
    c.Update();
    float range = c.GetUymax() - c.GetUymin();
    mg.GetYaxis()->SetRangeUser(c.GetUymin() - 0.5*range, c.GetUymax() + 0.5*range);
    c.Modified(); c.Update();

    // draw an axis on the right side
    auto axis = TGaxis(c.GetUxmax(),c.GetUymax(), c.GetUxmax(), c.GetUymin(), omega2pt(c.GetUymax()), omega2pt(c.GetUymin()), 510, "L-");
    axis.SetTitle("p_{T} (GeV/c)");
    axis.SetTitleOffset(1.5);
    axis.CenterTitle(true);
    axis.SetTitleColor(1);
    axis.SetTitleFont(42);
    axis.SetTitleSize(0.06); // ILD 0.07
    // For the axis labels:
    axis.SetLabelColor(1);
    axis.SetLabelFont(42);
    axis.SetLabelOffset(0.007); // ILD 0.015
    axis.SetLabelSize(0.05); // ILD 0.06

    // For the axis:
    axis.SetTickLength(0.03);
    axis.Draw();

    TLegend legend = TLegend(0.52, 0.755, 0.765, 0.895);
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    legend.SetTextFont(62);

    legend.SetHeader(Form("PDG: %d", pdg));
    legend.AddEntry(&grOmega, "Reconstructed", "l");
    legend.AddEntry(&grTrueOmega, "True", "l");
    legend.Draw();
    c.Modified(); c.Update();
    c.SaveAs( Form("./plots/track_params/pdg_%d_pt_%2.0f_pz_%2.0f_omega.png", pdg, mcMom.rho()*1000, mcMom.z()*1000) );
    // while ( !gSystem->ProcessEvents() ){
    //     gSystem->Sleep(5);
    // }


    // COPYPASTE for TanL
    TCanvas c2 = TCanvas("c2", "Track parameters evolution through the track", 600/(1. - margin.left - margin.right), 600/(1. - margin.top - margin.bottom));
    c2.SetGridx(true);

    grTanL.SetLineColor(TColor::GetColor("#377eb8"));
    grTanL.SetMarkerColor(TColor::GetColor("#377eb8"));
    grTanL.SetLineWidth(3);
    grTrueTanL.SetMarkerColor(TColor::GetColor("#e41a1c"));
    grTrueTanL.SetLineColor(TColor::GetColor("#e41a1c"));
    grTrueTanL.SetLineWidth(3);
    TMultiGraph mg2;
    mg2.Add(&grTanL, "PL");
    mg2.Add(&grTrueTanL, "PL");
    mg2.GetXaxis()->SetTitle("z (mm)");
    mg2.GetXaxis()->SetNdivisions(508);
    mg2.GetXaxis()->SetNoExponent(true);
    mg2.GetYaxis()->SetTitle("tan #lambda");
    mg2.GetYaxis()->SetTitleOffset(1.5);
    mg2.GetXaxis()->CenterTitle(true);
    mg2.GetYaxis()->CenterTitle(true);
    mg2.Draw("APL");
    c2.Update();
    float range2 = c2.GetUymax() - c2.GetUymin();
    mg2.GetYaxis()->SetRangeUser(c2.GetUymin() - 0.5*range2, c2.GetUymax() + 0.5*range2);
    c2.Modified(); c2.Update();

    // draw an axis on the right side
    auto axis2 = TGaxis(c2.GetUxmax(),c2.GetUymax(), c2.GetUxmax(), c2.GetUymin(), tanL2theta(c2.GetUymax()), tanL2theta(c2.GetUymin()), 510, "L-");
    axis2.SetMaxDigits(2);
    axis2.SetTitle("#theta (#circ)");
    axis2.SetTitleOffset(1.5);
    axis2.CenterTitle(true);
    axis2.SetTitleColor(1);
    axis2.SetTitleFont(42);
    axis2.SetTitleSize(0.06); // ILD 0.07
    // For the axis labels:
    axis2.SetLabelColor(1);
    axis2.SetLabelFont(42);
    axis2.SetLabelOffset(0.007); // ILD 0.015
    axis2.SetLabelSize(0.05); // ILD 0.06

    // For the axis:
    axis2.SetTickLength(0.03);
    axis2.Draw();

    TLegend legend2 = TLegend(0.52, 0.755, 0.765, 0.895);
    legend2.SetFillStyle(0);
    legend2.SetBorderSize(0);
    legend2.SetTextFont(62);

    legend2.SetHeader(Form("PDG: %d", pdg));
    legend2.AddEntry(&grTanL, "Reconstructed", "l");
    legend2.AddEntry(&grTrueTanL, "True", "l");
    legend2.Draw();
    c2.Modified(); c2.Update();
    c2.SaveAs( Form("./plots/track_params/pdg_%d_pt_%2.0f_pz_%2.0f_tanl.png", pdg, mcMom.rho()*1000, mcMom.z()*1000) );
    while ( !gSystem->ProcessEvents() ){
        gSystem->Sleep(5);
    }
    picture_counter++;
}

void displayTOFExplanation(std::vector<EVENT::CalorimeterHit*> allHits, std::vector<EVENT::CalorimeterHit*> selectedHits, double x, double y, double z, double px, double py, double pz){

    int type = 0; // point
    int size = 4; // larger point
    unsigned long color = 0x000000;
    ced_line(x, y, z, x+px*100, y+py*100, z+pz*100 , type , 1., color);

    for (auto hit:allHits){
        auto pos = hit->getPosition();
        bool isSelectedHit = std::find(selectedHits.begin(), selectedHits.end(), hit) != selectedHits.end();
        ced_hit(pos[0], pos[1], pos[2], type, size*2, isSelectedHit ? 0xfc0303 : 0x000000);
    }
}

