import numpy as np
import ROOT
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from utils import *

VAR_NAME = 'm2'
PLOT_NAME = "dtdl" # approx
if VAR_NAME == 'm':
    y_title = "Mass (GeV/c^{2})"
    y_min, y_max = 0, 3.
    lin_uncert_legend = " m#[]{#frac{#Deltap}{p} \u2295 #gamma^{2}#(){#frac{#DeltaL}{L} \u2295 #frac{#DeltaT}{T}}}"
    total_uncert_legend = " m(p#pm#Deltap, L#mp#DeltaL, T#pm#DeltaT)"
    x_text_pos, y0_text_pos, dy = 0.44, 1.77, 0.2
    legend = ROOT.TLegend(0.2, 0.71, 0.84, 0.92)
    legend.SetColumnSeparation(0.2)
elif VAR_NAME == 'm2':
    y_title = "Mass^{2} (GeV^{2}/c^{4})"
    y_min, y_max = -1., 3.
    lin_uncert_legend = " 2m^{2}#[]{#frac{#Deltap}{p} \u2295 #gamma^{2}#(){#frac{#DeltaL}{L} \u2295 #frac{#DeltaT}{T}}}"
    total_uncert_legend = "m^{2}(p#pm#Deltap, L#mp#DeltaL, T#pm#DeltaT)"
    x_text_pos, y0_text_pos, dy = 0.27, 1.82, 0.22
    legend = ROOT.TLegend(0.2, 0.75, 0.76, 0.91)
    legend.SetColumnSeparation(0.05)

if PLOT_NAME == 'dtdl':
    pass
elif PLOT_NAME == 'approx':
    pass

ROOT.gStyle.SetCanvasPreferGL(True)
plt.rcParams.update({'font.size':18, 'text.usetex':True})

margin = 0.22
ROOT.gStyle.SetPadLeftMargin(0.8*margin)
ROOT.gStyle.SetPadRightMargin(0.2*margin)
ROOT.gStyle.SetPadTopMargin(0.3*margin)
ROOT.gStyle.SetPadBottomMargin(0.7*margin)
canvas = ROOT.TCanvas(get_rand_string(),"", 600, 600)
legend.SetNColumns(2)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetMargin(0.2)
legend.SetTextFont(62)

fill_alpha = .3

momentum, track_length = np.arange(0.1, 8, 0.01), 2200.
# # STAR DATA from the paper https://arxiv.org/abs/nucl-ex/0308022
dp = 0.013 # relative
dl = 0.002 # relative
# dt = 0.100 # ns
dt = 0.020 # ns

particles = create_list_of_particles(momentum, track_length)
(pion, kaon, proton) = particles


graphs = {}
if PLOT_NAME == "approx":
    for p in particles:
        m_down_approx, m_up_approx = get_linear_uncertainty(p, dp, dl, dt, func=VAR_NAME)
        m_down, m_up = get_uncertainty(p, dp, dl, dt, func=VAR_NAME)
        graphs[p, 'approx'] = get_filled_graph(momentum, m_down_approx, m_up_approx)
        graphs[p, 'total'] = get_filled_graph(momentum, m_down, m_up)
        graphs[p, 'mass'] = ROOT.TGraph(len(momentum), momentum, np.full_like(momentum, p.mass))

        if p.name == 'pion':
            graphs[p, 'mass'].Draw("AL")
            graphs[p, 'mass'].GetXaxis().SetTitle("Momentum (GeV/c)")
            graphs[p, 'mass'].GetYaxis().SetTitle(y_title)
            graphs[p, 'mass'].GetYaxis().SetTitleOffset(1.2)
            graphs[p, 'mass'].GetXaxis().SetRangeUser(0, 8)
            graphs[p, 'mass'].GetYaxis().SetRangeUser(y_min, y_max)

        else:
            graphs[p, 'mass'].Draw("L same")

        graphs[p, 'mass'].SetLineColor( p.color )
        graphs[p, 'approx'].SetLineWidth(3)
        graphs[p, 'approx'].SetLineStyle(9)
        graphs[p, 'approx'].SetLineColor(p.color)
        graphs[p, 'total'].SetFillColorAlpha(p.color, fill_alpha)
        graphs[p, 'approx'].Draw("Lsame")
        graphs[p, 'total'].Draw("fsame")



    #row1 col1
    gr_legend_dl = ROOT.TGraph()
    gr_legend_dl.SetLineStyle(9)
    gr_legend_dl.SetLineWidth(3)
    legend.AddEntry(gr_legend_dl, lin_uncert_legend, "l")

    #row1 col2
    gr_legend_pi = ROOT.TGraph()
    gr_legend_pi.SetFillColor(pion.color)
    legend.AddEntry(gr_legend_pi, pion.legend, "f")

    #row2 col1
    gr_legend_true = ROOT.TGraph()
    gr_legend_true.SetLineColor(ROOT.kBlack)
    legend.AddEntry("", "", "")

    #row2 col2
    gr_legend_k = ROOT.TGraph()
    gr_legend_k.SetFillColor(kaon.color)
    legend.AddEntry(gr_legend_k, kaon.legend, "f")

    #row3 col1
    gr_legend_dldt = ROOT.TGraph()
    gr_legend_dldt.SetFillColor(ROOT.kBlack)
    gr_legend_dldt.SetFillColorAlpha(ROOT.kBlack, fill_alpha)
    legend.AddEntry(gr_legend_dldt, total_uncert_legend, "f")

    #row3 col2
    gr_legend_p = ROOT.TGraph()
    gr_legend_p.SetFillColor(proton.color)
    legend.AddEntry(gr_legend_p, proton.legend, "f")


    legend.Draw()

elif PLOT_NAME == "dtdl":
    for p in particles:
        m_down_dt, m_up_dt = get_uncertainty(p, 0., 0., dt, func=VAR_NAME)
        m_down_dl, m_up_dl = get_uncertainty(p, 0., dl, 0., func=VAR_NAME)
        graphs[p, 'dt'] = get_filled_graph(momentum, m_down_dt, m_up_dt)
        graphs[p, 'dl'] = get_filled_graph(momentum, m_down_dl, m_up_dl)
        graphs[p, 'mass'] = ROOT.TGraph(len(momentum), momentum, np.full_like(momentum, p.mass))

        if p.name == 'pion':
            graphs[p, 'mass'].Draw("AL")
            graphs[p, 'mass'].GetXaxis().SetTitle("Momentum (GeV/c)")
            graphs[p, 'mass'].GetYaxis().SetTitle(y_title)
            graphs[p, 'mass'].GetYaxis().SetTitleOffset(1.2)
            graphs[p, 'mass'].GetXaxis().SetRangeUser(0, 8)
            graphs[p, 'mass'].GetYaxis().SetRangeUser(y_min, y_max)
        else:
            graphs[p, 'mass'].Draw("L same")

        graphs[p, 'mass'].SetLineColor( p.color )
        graphs[p, 'dl'].SetLineWidth(3)
        graphs[p, 'dl'].SetLineStyle(9)
        graphs[p, 'dl'].SetLineColor(p.color)
        graphs[p, 'dt'].SetFillColorAlpha(p.color, fill_alpha)
        graphs[p, 'dl'].Draw("Lsame")
        graphs[p, 'dt'].Draw("fsame")



    #row1 col1
    gr_legend_dl = ROOT.TGraph()
    gr_legend_dl.SetLineStyle(9)
    gr_legend_dl.SetLineWidth(3)
    legend.AddEntry(gr_legend_dl, "#Delta L", "l")

    #row1 col2
    gr_legend_pi = ROOT.TGraph()
    gr_legend_pi.SetFillColor(pion.color)
    legend.AddEntry(gr_legend_pi, pion.legend, "f")

    #row2 col1
    gr_legend_true = ROOT.TGraph()
    gr_legend_true.SetLineColor(ROOT.kBlack)
    legend.AddEntry("", "", "")

    #row2 col2
    gr_legend_k = ROOT.TGraph()
    gr_legend_k.SetFillColor(kaon.color)
    legend.AddEntry(gr_legend_k, kaon.legend, "f")

    #row3 col1
    gr_legend_dldt = ROOT.TGraph()
    gr_legend_dldt.SetFillColor(ROOT.kBlack)
    gr_legend_dldt.SetFillColorAlpha(ROOT.kBlack, fill_alpha)
    legend.AddEntry(gr_legend_dldt, "#Delta T", "f")

    #row3 col2
    gr_legend_p = ROOT.TGraph()
    gr_legend_p.SetFillColor(proton.color)
    legend.AddEntry(gr_legend_p, proton.legend, "f")


    legend.Draw()

latex.DrawLatex(x_text_pos, y0_text_pos, f"L = {track_length:.0f}" + " mm")
if PLOT_NAME == "approx":
    latex.DrawLatex(x_text_pos, y0_text_pos - dy, "#Deltap/p = " + f"{dp:.3f}")
    latex.DrawLatex(x_text_pos, y0_text_pos - 2*dy, "#DeltaL/L = " + f"{dl:.3f}")
    latex.DrawLatex(x_text_pos, y0_text_pos - 3*dy, "#DeltaT = " + f"{1000*dt:.0f}" + " ps")
elif PLOT_NAME == "dtdl":
    latex.DrawLatex(x_text_pos, y0_text_pos - dy, "#DeltaL/L = " + f"{dl:.3f}")
    latex.DrawLatex(x_text_pos, y0_text_pos - 2*dy, "#DeltaT = " + f"{1000*dt:.0f}" + " ps")

canvas.Update()
input("wait")
