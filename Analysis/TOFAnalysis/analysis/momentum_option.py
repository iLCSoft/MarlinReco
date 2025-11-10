import ROOT
import numpy as np
from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(0)
ROOT.EnableImplicitMT()

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
        .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212").Filter("tofClosest0 > 6.")\
        .Define("beta", "trackLengthToEcal_IKF_zedLambda/(tofClosest0*299.792458)")\
        .Define("mom2_ip", "recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz")\
        .Define("mom_ip", "sqrt(mom2_ip)")\
        .Define("mom2_ecal", "recoCaloPx*recoCaloPx + recoCaloPy*recoCaloPy + recoCaloPz*recoCaloPz")\
        .Define("mom_ecal", "sqrt(mom2_ecal)")\
        .Define("mom2_hm", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda")\
        .Define("mom_hm", "harmonicMomToEcal_IKF_zedLambda")\
        .Define("mass2_ip", "mom2_ip*(1./(beta*beta) - 1.)")\
        .Define("mass2_ecal", "mom2_ecal*(1./(beta*beta) - 1.)")\
        .Define("mass2_hm", "mom2_hm*(1./(beta*beta) - 1.)")

options = ["ip", "ecal", "hm"]

histos_2d = []
histos_1d_pi = []
histos_1d_k = []
histos_1d_p = []
for opt in options:
    h2d = df.Histo2D( (get_rand_string(), "", 1000, 0, 10, 1000, -0.3, 1.2), f"mom_{opt}", f"mass2_{opt}" )
    h1d_pi = df.Histo1D( (get_rand_string(), "", 500, 16e-3, 22e-3), f"mass2_{opt}" )
    h1d_k = df.Histo1D( (get_rand_string(), "", 250, 0.22, 0.26), f"mass2_{opt}" )
    h1d_p = df.Histo1D( (get_rand_string(), "", 250, 0.84, .92), f"mass2_{opt}" )

    histos_2d.append(h2d)
    histos_1d_pi.append(h1d_pi)
    histos_1d_k.append(h1d_k)
    histos_1d_p.append(h1d_p)

canvases = []
legends = []
lls = []
# for h2d in histos_2d:
#     h2d.GetXaxis().SetTitle("Momentum (GeV/c)")
#     h2d.GetYaxis().SetTitle("Mass^{2} (GeV^{2}/c^{4})")
#     canvases.append( draw_2d_plot(h2d) )


colors = [ROOT.TColor.GetColor('#000000'), ROOT.TColor.GetColor('#688E26'), ROOT.TColor.GetColor('#FAA613')]
for p, h_to_plot in zip(particles, [histos_1d_pi, histos_1d_k, histos_1d_p]):
    canvas = create_canvas(margin=0.33, left_margin_fraction=0.55, bottom_margin_fraction=0.65)
    max_y = 1.18*max([hist.GetMaximum() for hist in h_to_plot])
    for h, color in zip(h_to_plot, colors):
        h.SetTitle(";Mass^{2} (GeV^{2}/c^{4}); N entries")
        h.GetXaxis().SetNdivisions(506)
        h.GetXaxis().SetMaxDigits(3)
        h.GetXaxis().SetTitleOffset(1.1)
        h.GetYaxis().SetTitleOffset(1.4)
        h.GetYaxis().SetMaxDigits(3)
        h.SetLineColor(color)
        h.SetLineWidth(3)
        h.SetMinimum( 0. )
        h.SetMaximum( max_y )
        h.Draw("L") if h == h_to_plot[0] else h.Draw("L same")

    legend = create_legend(0.205, 0.58, 0.47, 0.88)

    legend.SetTextFont(42)

    leg1 = ROOT.TGraph()
    leg1.SetFillColor(colors[0])
    legend.AddEntry(leg1, "p_{ip}", "f")

    leg3 = ROOT.TGraph()
    leg3.SetFillColor(colors[1])
    legend.AddEntry(leg3, "p_{ECAL}", "f")

    leg5 = ROOT.TGraph()
    leg5.SetFillColor(colors[2])
    legend.AddEntry(leg5, "p_{HM}", "f")

    lines = draw_vertical_mass_lines(max_y)
    legend.AddEntry(lines[p], f"true {p.legend} mass", "l")

    legend.Draw()
    lls.append(lines)
    # latex.SetTextSize(0.04)
    # latex.DrawLatexNDC(0.2, 0.8, text)
    canvas.Update()
    canvases.append(canvas)
    legends.append(legend)

input("wait")












# colors = [ROOT.TColor.GetColor('#ff7f00') ,ROOT.TColor.GetColor('#984ea3') ,ROOT.TColor.GetColor('#4daf4a') ,ROOT.TColor.GetColor('#377eb8') ,ROOT.TColor.GetColor('#e41a1c')]
colors = [ROOT.TColor.GetColor('#1b9e77'), ROOT.TColor.GetColor('#d95f02'), ROOT.TColor.GetColor('#7570b3'), ROOT.TColor.GetColor('#e7298a')]

def draw_lines():
    lines = {}
    pdgs = [211, 321, 2212]
    m_pdg = {211 : 0.13957039*1000, 321 : 0.493677*1000, 2212 : 0.938272088*1000}
    for pdg in pdgs:
        lines[pdg] = ROOT.TLine(0., m_pdg[pdg], 15., m_pdg[pdg])
        lines[pdg].SetLineColor(2)
        lines[pdg].SetLineWidth(1)
        lines[pdg].SetLineStyle(9)
        lines[pdg].Draw()
    return lines

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
df = df.Filter("tofClosest0 > 6.").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")

#PLOT 1D
histos = []
df = df.Define("beta", "trackLengthToEcal_IKF_zedLambda/(tofClosest0*299.792458)").Filter("beta >= 0 && beta <= 1").\
        Define("mom_ip", "sqrt(recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz)").\
        Define("mom_ecal", "sqrt(recoCaloPx*recoCaloPx + recoCaloPy*recoCaloPy + recoCaloPz*recoCaloPz)").\
        Define("mom_hm", "harmonicMomToEcal_IKF_zedLambda")

h_ip = df.Define("mass", "mom_ip*sqrt( 1./(beta*beta) - 1.)*1000").Histo1D(("h_ip", "; mass [MeV]; N entries", 500, 130, 150), "mass")
h_ecal = df.Define("mass", "mom_ecal*sqrt( 1./(beta*beta) - 1.)*1000").Histo1D(("h_ecal", "; mass [MeV]; N entries", 500, 130, 150), "mass")
h_hm = df.Define("mass", "mom_hm*sqrt( 1./(beta*beta) - 1.)*1000").Histo1D(("h_hm", "; mass [MeV]; N entries", 500, 130, 150), "mass")
histos.append(h_ip)
histos.append(h_ecal)
histos.append(h_hm)

ROOT.gStyle.SetPadLeftMargin(0.18)
ROOT.gStyle.SetPadRightMargin(0.08)
canvas = ROOT.TCanvas("mass_1d_mom_options",
                        "",
                        int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
                        int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
canvas.SetGridx()
canvas.SetGridy()

legend = ROOT.TLegend()

names = ["p_{IP}", "p_{ECAL}", "p_{HM}"]
for i, (name, h) in enumerate(zip(names, histos)):
    h.Draw("L" if i == 0 else "Lsame")
    h.SetLineColor(colors[i])
    h.SetLineWidth(3)
    legend.AddEntry(h.GetPtr(),name,"l")

legend.Draw()
histos[0].GetYaxis().SetTitleOffset(1.45)
histos[0].SetMinimum(0)
canvas.Modified()
canvas.Update()
input("wait")