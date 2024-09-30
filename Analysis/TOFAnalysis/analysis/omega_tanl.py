import ROOT
import numpy as np
from utils import *
ROOT.EnableImplicitMT()

MODE_NAME, VAR_NAME = "change", "tanL"
if MODE_NAME == "residual" and VAR_NAME == "omega":
    x_title, y_title = "#Omega_{reco} - #Omega_{true} (1/mm)", "Normalised N entries"
    n_bins, x_min, x_max = 100, -20e-6, 20e-6
    x_text_pos, y_text_pos = 2.58e-6, 13.34e-3
elif MODE_NAME == "residual" and VAR_NAME == "tanL":
    x_title, y_title = "tan(#lambda)_{reco} - tan(#lambda)_{true}", "Normalised N entries"
    n_bins, x_min, x_max = 100, -10e-3, 10e-3
    x_text_pos, y_text_pos = 3e-3, 10e-3
elif MODE_NAME == "change" and VAR_NAME == "omega":
    x_title, y_title = "#Omega_{ECAL} - #Omega_{IP} (1/mm)", "Normalised N entries"
    n_bins, x_min, x_max = 200, -0.5e-6, 1e-6
elif MODE_NAME == "change" and VAR_NAME == "tanL":
    x_title, y_title = "tan(#lambda)_{ECAL} - tan(#lambda)_{IP}", "Normalised N entries"
    n_bins, x_min, x_max = 200, -5e-3, 5e-3



df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
df = df.Filter("tofClosest0 > 6.").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")
df = df.Define("ptTrue", "sqrt(mcPx*mcPx + mcPy*mcPy)")
df = df.Define("ptIp", "sqrt(recoIpPx*recoIpPx + recoIpPy*recoIpPy)")
df = df.Define("ptCalo", "sqrt(recoCaloPx*recoCaloPx + recoCaloPy*recoCaloPy)")
df = df.Define("pTrue", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
df = df.Define("pIp", "sqrt(recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz)")
df = df.Define("pCalo", "sqrt(recoCaloPx*recoCaloPx + recoCaloPy*recoCaloPy + recoCaloPz*recoCaloPz)")
df = df.Define("omegaTrue", "1e-6 * 299.792458 * 3.5 / ptTrue")
df = df.Define("omegaIp", "1e-6 * 299.792458 * 3.5 / ptIp")
df = df.Define("omegaCalo", "1e-6 * 299.792458 * 3.5 / ptCalo")
df = df.Define("tanLTrue", "mcPz/ptTrue")
df = df.Define("tanLIp", " recoIpPz/ptIp")
df = df.Define("tanLCalo", "recoCaloPz/ptCalo")
df = df.Define("residual_omega", "omegaIp - omegaTrue")
df = df.Define("residual_tanL", "tanLIp - tanLTrue")
df = df.Define("change_omega", "omegaCalo - omegaIp")
df = df.Define("change_tanL", "tanLCalo - tanLIp")
if MODE_NAME == "residual":
    df = df.Filter("pTrue < 1.")

h_pion_default = df.Filter("abs(pdg) == 211").Histo1D(("h_pion_default", pion.legend, n_bins, x_min, x_max), f"{MODE_NAME}_{VAR_NAME}")
h_kaon_default = df.Filter("abs(pdg) == 321").Histo1D(("h_kaon_default", kaon.legend, n_bins, x_min, x_max), f"{MODE_NAME}_{VAR_NAME}")
h_proton_default = df.Filter("abs(pdg) == 2212").Histo1D(("h_proton_default", proton.legend, n_bins, x_min, x_max), f"{MODE_NAME}_{VAR_NAME}")

canvas = create_canvas(0.3, 0.55, 0.7)
legend =create_legend(0.2, 0.65, 0.62, 0.87)
histos = [h_pion_default, h_kaon_default, h_proton_default]

for h in histos:
    h.Scale(1./h.GetEntries())

max_y = 1.05*max([hist.GetMaximum() for hist in histos])
for p, h in zip(particles, histos):
    h.SetLineColor(p.color)
    h.SetLineWidth(3)
    legend.AddEntry(p.legend_graph, p.legend, "f")
    if p.name == "pion":
        h.Draw("histo")
        h.GetXaxis().SetTitle(x_title)
        h.GetXaxis().SetTitleOffset(1.1)
        h.GetXaxis().SetMaxDigits(3)
        h.GetXaxis().SetNdivisions(506)
        h.GetYaxis().SetTitle(y_title)
        h.GetYaxis().SetMaxDigits(3)
        h.SetMinimum(0.)
        h.SetMaximum(max_y)
    else:
        h.Draw("histo same")

if MODE_NAME == "residual":
    latex.DrawLatex(x_text_pos, y_text_pos, "p < 1 GeV/c")

legend.Draw()
canvas.Update()
input("wait")