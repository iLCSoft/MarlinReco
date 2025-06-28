import ROOT
import numpy as np
from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.EnableImplicitMT()


colors = ['#0444b3', '#0068cc', '#008add', '#1caaea', '#61caf4', '#98e8ff']
# colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99']
colors = [ ROOT.TColor.GetColor(c) for c in colors]

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root").Filter("abs(pdg) == 321")
df = df.Filter("tofClosest0 > 6.")
# ROOT.Experimental.AddProgressBar(df)
n_mass_bins, mass_min, mass_max = 2000, -0.5, 1

resolutions = [0, 10, 30, 100, 300]
tof_columns = [f"tofClosest{i}" for i in resolutions]
histos = []
for tof_column in tof_columns:
    df_beta = df.Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")
    df_mass2 = df_beta.Define("mass2", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda*( 1./(beta*beta) - 1.)")

    h = df_mass2.Histo1D((f"h_{tof_column}", "; Mass^{2} (GeV^{2}/c^{4}); N entries", n_mass_bins, mass_min, mass_max), "mass2")
    # df_mass = df_mass2.Filter("mass2>0").Define("mass", "sqrt(mass2)")
    # h = df_mass.Histo1D((f"h_{tof_column}", "; Mass (GeV/c^{2}); N entries", n_mass_bins, mass_min, mass_max), "mass")
    histos.append(h)

ROOT.gStyle.SetPadLeftMargin(0.17)
ROOT.gStyle.SetPadRightMargin(0.14)
ROOT.gStyle.SetPadTopMargin(0.08)
ROOT.gStyle.SetPadBottomMargin(0.17)
canvas = ROOT.TCanvas("m_vs_tofres_1d","",600,600)

legend = ROOT.TLegend(0.61, 0.63, 0.86, 0.9)
legend.SetFillStyle(0)
legend.SetBorderSize(0)


for i, h in enumerate(histos):
    h.Draw("" if i == 0 else "same")
    h.SetLineColor(colors[i])
    h.SetLineWidth(2)
    legend.AddEntry(h.GetPtr(), f"{resolutions[i]} ps","l")
    if i == 0:
        h.GetYaxis().SetTitleOffset(1.3)
        h.GetXaxis().SetMaxDigits(3)
        h.GetYaxis().SetMaxDigits(3)
        h.GetYaxis().SetRangeUser(0, 110000)

legend.Draw()

# line = ROOT.TLine(0.13957039*0.13957039, 0., 0.13957039*0.13957039, 110000.)
line = ROOT.TLine(kaon.mass2, 0., kaon.mass2, 110000.)
line.Draw()
line.SetLineWidth(3)
line.SetLineStyle(9)
line.SetLineColor(ROOT.kViolet-3)

canvas.Update()
input("wait")
