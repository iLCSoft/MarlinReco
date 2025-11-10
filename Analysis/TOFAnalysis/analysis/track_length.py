import ROOT
import numpy as np
from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.EnableImplicitMT()

def compare_IDR_vs_latest():
    histos_2d = []
    histos_1d_pi = []
    histos_1d_k = []
    histos_1d_p = []

    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
            .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212").Filter("tofClosest0 > 6.")\

    df_idr = df.Define("beta", "trackLength_IDR/(tofClosest0*299.792458)")\
            .Define("mom2", "recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz")\
            .Define("mom", "sqrt(mom2)")\
            .Define("mass2", "mom2*(1./(beta*beta) - 1.)")
    h2d = df_idr.Histo2D( (get_rand_string(), "", 1000, 0, 10, 1000, -0.3, 1.2), "mom","mass2" )
    h1d_pi = df_idr.Histo1D( (get_rand_string(), "", 500, -10e-3, 50e-3), "mass2" )
    h1d_k = df_idr.Histo1D( (get_rand_string(), "", 250, 0.19, 0.31), "mass2" )
    h1d_p = df_idr.Histo1D( (get_rand_string(), "", 250, 0.78, 1.), "mass2" )

    histos_2d.append(h2d)
    histos_1d_pi.append(h1d_pi)
    histos_1d_k.append(h1d_k)
    histos_1d_p.append(h1d_p)

    df_new = df.Define("beta", "trackLengthToEcal_IKF_zedLambda/(tofClosest0*299.792458)")\
            .Define("mom2", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda")\
            .Define("mom", "sqrt(mom2)")\
            .Define("mass2", "mom2*(1./(beta*beta) - 1.)")

    h2d = df_new.Histo2D( (get_rand_string(), "", 1000, 0, 10, 1000, -0.3, 1.2), "mom","mass2" )
    h1d_pi = df_new.Histo1D( (get_rand_string(), "", 500, -10e-3, 50e-3), "mass2" )
    h1d_k = df_new.Histo1D( (get_rand_string(), "", 250, 0.19, 0.31), "mass2" )
    h1d_p = df_new.Histo1D( (get_rand_string(), "", 250, 0.78, 1.), "mass2" )

    histos_2d.append(h2d)
    histos_1d_pi.append(h1d_pi)
    histos_1d_k.append(h1d_k)
    histos_1d_p.append(h1d_p)


    canvases = []
    legends = []
    lls = []

    for h2d in histos_2d:
        h2d.GetXaxis().SetTitle("Momentum (GeV/c)")
        h2d.GetYaxis().SetTitle("Mass^{2} (GeV^{2}/c^{4})")
        canvases.append( draw_2d_plot(h2d) )


    colors = [ROOT.TColor.GetColor('#000000'), ROOT.TColor.GetColor('#FAA613')]
    for p, h_to_plot in zip(particles, [histos_1d_pi, histos_1d_k, histos_1d_p]):
        canvas = create_canvas(margin=0.33, left_margin_fraction=0.55, bottom_margin_fraction=0.65)
        max_y = 1.18*max([hist.GetMaximum() for hist in h_to_plot])
        for h, color in zip(h_to_plot, colors):
            h.SetTitle(";Mass^{2} (GeV^{2}/c^{4}); N entries")
            h.GetXaxis().SetMaxDigits(3)
            h.GetXaxis().SetTitleOffset(1.1)
            h.GetYaxis().SetTitleOffset(1.4)
            h.GetYaxis().SetMaxDigits(3)
            h.SetLineColor(color)
            h.SetLineWidth(3)
            h.SetMinimum( 0. )
            h.SetMaximum( max_y )
            h.Draw("L") if h == h_to_plot[0] else h.Draw("L same")

        legend = create_legend(0.555, 0.435, 0.82, 0.88)

        legend.SetTextFont(42)

        leg1 = ROOT.TGraph()
        leg1.SetFillColor(colors[0])
        legend.AddEntry(leg1, "ILD IDR (2020)", "f")

        leg3 = ROOT.TGraph()
        leg3.SetFillColor(colors[1])
        legend.AddEntry(leg3, "The best in this study", "f")

        lines = draw_vertical_mass_lines(max_y)
        legend.AddEntry(lines[p], f"true {p.legend} mass", "l")

        legend.Draw()
        lls.append(lines)
        canvas.Update()
        canvases.append(canvas)
        legends.append(legend)

    input("wait")

compare_IDR_vs_latest()

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
        .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212").Filter("tofClosest0 > 6.")\
        .Define("mom2HM", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda")
        # .Define("mom2", "recoIpPx*recoIpPx + recoIpPy*recoIpPy + recoIpPz*recoIpPz")\
        # .Define("mom", "sqrt(mom2)")

trk_len_methods = ["trackLengthToEcal_IKF_zedLambda"]
text = ""

# trk_len_methods = ["trackLengthToEcal_IKF_phiLambda", "trackLengthToEcal_IKF_phiZed", "trackLengthToEcal_IKF_zedLambda"]
# text = "L=#sum_{i}L_{i}(#Omega_{i}, tan#lambda_{i})"

# trk_len_methods = ["trackLengthToEcal_SHA_phiLambda_ECAL", "trackLengthToEcal_SHA_phiZed_ECAL", "trackLengthToEcal_SHA_zedLambda_ECAL"]
# text = "#Omega_{ECAL}, tan#lambda_{ECAL}"

histos_2d = []
histos_1d_pi = []
histos_1d_k = []
histos_1d_p = []
for trk_len in trk_len_methods:
    df_new = df.Define("beta", f"{trk_len}/(tofClosest0*299.792458)")\
               .Define("mass2", "mom2HM*(1./(beta*beta) - 1.)")
            #    .Define("mass2", "mom2HM*2.*(1./beta - 1.)")
    h2d = df_new.Histo2D( (get_rand_string(), "", 1000, 0, 10, 1000, -0.3, 1.2), "harmonicMomToEcal_IKF_zedLambda","mass2" )
    h1d_pi = df_new.Histo1D( (get_rand_string(), "", 500, -10e-3, 50e-3), "mass2" )
    h1d_k = df_new.Histo1D( (get_rand_string(), "", 250, 0.19, 0.31), "mass2" )
    h1d_p = df_new.Histo1D( (get_rand_string(), "", 250, 0.78, 1.), "mass2" )

    histos_2d.append(h2d)
    histos_1d_pi.append(h1d_pi)
    histos_1d_k.append(h1d_k)
    histos_1d_p.append(h1d_p)

canvases = []
legends = []
lls = []
for h2d in histos_2d:
    h2d.GetXaxis().SetTitle("Momentum (GeV/c)")
    h2d.GetYaxis().SetTitle("Mass^{2} (GeV^{2}/c^{4})")
    canvases.append( draw_2d_plot(h2d) )


colors = [ROOT.TColor.GetColor('#000000'), ROOT.TColor.GetColor('#688E26'), ROOT.TColor.GetColor('#FAA613')]
for p, h_to_plot in zip(particles, [histos_1d_pi, histos_1d_k, histos_1d_p]):
    canvas = create_canvas(margin=0.33, left_margin_fraction=0.55, bottom_margin_fraction=0.65)
    max_y = 1.18*max([hist.GetMaximum() for hist in h_to_plot])
    for h, color in zip(h_to_plot, colors):
        h.SetTitle(";Mass^{2} (GeV^{2}/c^{4}); N entries")
        h.GetXaxis().SetMaxDigits(3)
        h.GetXaxis().SetTitleOffset(1.1)
        h.GetYaxis().SetTitleOffset(1.4)
        h.GetYaxis().SetMaxDigits(3)
        h.SetLineColor(color)
        h.SetLineWidth(3)
        h.SetMinimum( 0. )
        h.SetMaximum( max_y )
        h.Draw("L") if h == h_to_plot[0] else h.Draw("L same")

    legend = create_legend(0.555, 0.435, 0.82, 0.88)

    legend.SetTextFont(42)

    leg1 = ROOT.TGraph()
    leg1.SetFillColor(colors[0])
    legend.AddEntry(leg1, "#left|#frac{#Delta#varphi}{#Omega}#right|#sqrt{1+tan^{2}#lambda}", "f")

    leg3 = ROOT.TGraph()
    leg3.SetFillColor(colors[1])
    legend.AddEntry(leg3, "#sqrt{#left(#frac{#Delta#varphi}{#Omega}#right)^{2} + #Deltaz^{2}}", "f")

    leg5 = ROOT.TGraph()
    leg5.SetFillColor(colors[2])
    legend.AddEntry(leg5, "#left|#frac{#Deltaz}{tan#lambda}#right|#sqrt{1+tan^{2}#lambda}", "f")

    lines = draw_vertical_mass_lines(max_y)
    legend.AddEntry(lines[p], f"true {p.legend} mass", "l")

    legend.Draw()
    lls.append(lines)
    latex.SetTextSize(0.04)
    latex.DrawLatexNDC(0.2, 0.8, text)
    canvas.Update()
    canvases.append(canvas)
    legends.append(legend)

input("wait")

