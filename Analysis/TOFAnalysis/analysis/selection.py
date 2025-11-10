import ROOT
import numpy as np
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetPaintTextFormat("3.0f")
ROOT.EnableImplicitMT()

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")
df = df.Filter("tofClosest0 > 6.")
df = df.Define("beta", "trackLengthToEcal_IKF_zedLambda/(tofClosest0*299.792458)")
df = df.Define("mass2", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda*( 1./(beta*beta) - 1.)")

#DONT FORGET TO REMOVE THIS
# df = df.Filter("cleanClosestHit == 1 && cleanTrack == 1 && layerClosest == 0")

h1 = df.Histo1D(("h1", ";clean track; %", 2, 0, 2), "cleanTrack")
h2 = df.Histo1D(("h2", ";clean closest hit; %", 2, 0, 2), "cleanClosestHit")
h3 = df.Histo1D(("h3", ";layer closest; %", 15,0,15), "layerClosest")
h4 = df.Histo1D(("h4", ";calo type; %", 15,0,15), "typeClosest")
h5 = df.Histo1D(("h5", ";calo id; %", 15,0,15), "caloIDClosest")
h6 = df.Histo1D(("h6", ";calo layout; %", 15,0,15), "layoutClosest")
# histos = [h1, h2, h3, h4, h5, h6]
# canvas = ROOT.TCanvas()
# for h in histos:
#     h.Scale(100./h.GetEntries())
#     h.SetLineWidth(3)
#     h.GetYaxis().SetRangeUser(0., 100.)
#     h.Draw("histo text45")
#     canvas.Update()
#     input("wait")

# ROOT.Experimental.AddProgressBar(df)
n_mom_bins, mom_min, mom_max = 2000, 0, 20
n_mass_bins, mass_min, mass_max = 2000, -3, 3


cuts = ["true","cleanClosestHit == 1", "cleanClosestHit == 1 && cleanTrack == 1", "cleanClosestHit == 1 && cleanTrack == 1 && layerClosest == 0"]
# cuts = ["true","cleanClosestHit != 1", "cleanTrack != 1", "layerClosest != 0"]
cuts = ["true","cleanClosestHit == 1 && cleanTrack == 1 && layerClosest == 0 && layoutClosest == 4"]

ROOT.gStyle.SetPadLeftMargin(0.17)
ROOT.gStyle.SetPadRightMargin(0.14)
ROOT.gStyle.SetPadTopMargin(0.04)
ROOT.gStyle.SetPadBottomMargin(0.14)
canvas2 = ROOT.TCanvas("c_2d_m_vs_p_total","",600,600)
histos ={}
n_events ={}
for cut in cuts:
    df_sel = df.Filter(cut)
    histos[cut] = df_sel.Histo2D((f"h2_{cut}", "; Momentum (GeV/c); Mass^{2} (GeV^{2}/c^{4})", n_mom_bins, mom_min, mom_max, n_mass_bins, mass_min, mass_max), "harmonicMomToEcal_IKF_zedLambda","mass2")
    n_events[cut] = df_sel.Count()

print("Total events: ", n_events["true"].GetValue())
for cut in cuts:
    print("Passing ", cut, 100.*n_events[cut].GetValue()/n_events["true"].GetValue())
    # DRAW FANCY 2D plot
    histos[cut].Draw("colz")
    histos[cut].GetYaxis().SetTitleOffset(1.1)
    print(histos[cut].GetMaximum())
    histos[cut].SetMinimum(1)
    histos[cut].SetMaximum(10000)
    canvas2.SetLogz()
    canvas2.SetGridx(0)
    canvas2.SetGridy(0)
    canvas2.Update()
    palette = histos[cut].GetListOfFunctions().FindObject("palette")
    palette.SetX1(mom_min + 1.01*(mom_max-mom_min))
    palette.SetX2(mom_min + 1.05*(mom_max-mom_min))
    palette.SetLabelOffset(0.001)
    canvas2.Modified()
    canvas2.Update()
    input("wait")
