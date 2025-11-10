import ROOT
import numpy as np
from random import choice
from string import ascii_letters

ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)

DEDX_COLOR = ROOT.TColor.GetColor("#c52700")


SPEED_OF_LIGHT = 299.792458 # mm / ns
PION_MASS = 0.13957039 # GeV/ c^{2}
KAON_MASS = 0.493677 # GeV/ c^{2}
PROTON_MASS = 0.93827208816 # GeV/ c^{2}

latex = ROOT.TLatex()
latex.SetTextFont(52)
latex.SetTextSize(0.04)

class Particle:
    def __init__(self, name, short_name, legend, color, mass, momentum, track_length, pdg):
        self.name = name
        self.short_name = short_name
        self.legend = legend
        self.color = color
        self.mass = mass
        self.mass2 = mass*mass
        self.momentum = momentum
        self.track_length = track_length
        self.pdg = pdg
        self.tof = track_length/SPEED_OF_LIGHT * np.sqrt( 1 + (mass/momentum)**2 )
        self.beta = track_length/(SPEED_OF_LIGHT * self.tof)
        self.gamma = 1/np.sqrt(1 - self.beta**2)
        self.legend_graph = ROOT.TGraph()
        self.legend_graph.SetFillColor(self.color)

def create_list_of_particles(momentum, track_length):
    pion = Particle('pion', "pi", ' #pi^{#pm}', ROOT.TColor.GetColor('#1b9e77'), PION_MASS, momentum, track_length, 211)
    kaon = Particle('kaon', "k", ' K^{#pm}', ROOT.TColor.GetColor('#d95f02'), KAON_MASS, momentum, track_length, 321)
    proton = Particle('proton', "p", ' p', ROOT.TColor.GetColor('#7570b3'), PROTON_MASS, momentum, track_length, 2212)
    return [pion, kaon, proton]

particles = create_list_of_particles(1., 1) # just for colours! overwrite when studying mass_uncertainty!
(pion, kaon, proton) = particles

def get_rand_string():
    return ''.join(choice(ascii_letters) for i in range(16))


def get_filled_graph(x, y_min, y_max):
    '''
    Return a TGraph filled between y_min and y_max.
    '''
    n_points = len(x)
    graph = ROOT.TGraph(2*n_points)
    for i in range(n_points):
        graph.SetPoint(i, x[i], y_max[i])
        graph.SetPoint(n_points + i, x[n_points-i-1], y_min[n_points-i-1])
    return graph


def get_linear_uncertainty(particle, dp, dl, dt, func="m"):
    '''
    dp, dl - relative. dt - absolute
    '''
    common_uncrt = np.sqrt(dp**2 + (dl**2 + (dt/particle.tof)**2)*particle.gamma**4)
    f0, df = 0, 0
    if func == "m":
        f0 = particle.mass
        df = f0 * common_uncrt
    elif func == "m2":
        f0 = particle.mass*particle.mass
        df = 2 * f0 * common_uncrt
    return np.nan_to_num(f0 - df, nan=0), np.nan_to_num(f0 + df, nan=0)


def get_uncertainty(particle, dp, dl, dt, func="m"):
    '''
    dp, dl - relative. dt - absolute
    Return mass of the particle in momentum bins for the plot
    '''
    mom_down, mom_up = particle.momentum*(1. - dp), particle.momentum*(1. + dp)
    track_length_down, track_length_up = particle.track_length*(1. - dl), particle.track_length*(1. + dl)
    tof_down, tof_up = particle.tof - dt, particle.tof + dt

    beta_down, beta_up = track_length_down/(tof_up*SPEED_OF_LIGHT), track_length_up/(tof_down*SPEED_OF_LIGHT)
    m2_down, m2_up = mom_down*mom_down*(1./(beta_up*beta_up) - 1.), mom_up*mom_up*(1./(beta_down*beta_down) - 1.)
    if func == "m":
        return np.nan_to_num(np.sqrt(m2_down), nan=0), np.nan_to_num(np.sqrt(m2_up), nan=0)
    return m2_down, m2_up



######################## DRAWING #################

def draw_2d_plot(h, margin=0.33, left_margin_fraction=0.58, bottom_margin_fraction=0.65):
    '''Draw a 2D histogram in the appropriate styling and pause.'''
    canvas = create_canvas(margin, left_margin_fraction, bottom_margin_fraction)

    h.Draw("colz")
    h.GetXaxis().SetTitleOffset(1.1)
    if h.GetYaxis().GetTitleOffset() < 1.4:
         h.GetYaxis().SetTitleOffset(1.4)

    print(f"The highest bin is: {h.GetMaximum()}", end=" ")

    canvas.SetLogz()
    canvas.SetGridx(0)
    canvas.SetGridy(0)
    canvas.Update()

    palette = h.GetListOfFunctions().FindObject("palette")
    x_min = h.GetXaxis().GetXmin()
    x_max = h.GetXaxis().GetXmax()
    palette.SetX1(x_min + 1.01*(x_max-x_min))
    palette.SetX2(x_min + 1.05*(x_max-x_min))
    palette.SetMaxDigits(3)
    palette.SetLabelOffset(0.006)

    canvas.Modified()
    canvas.Update()
    return canvas

def create_legend(x1=0.2, y1=0.75, x2=0.76, y2=0.91):
    legend = ROOT.TLegend(x1, y1, x2, y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(62)
    return legend

def create_canvas(margin=0.22, left_margin_fraction=0.8, bottom_margin_fraction=0.7):
    '''
    Create a square canvas 1200x1200 with equal horizontal and vertical margins to ensure also a square plot inside the axes!
    '''
    ROOT.gStyle.SetPadLeftMargin( left_margin_fraction * margin )
    ROOT.gStyle.SetPadRightMargin( (1 - left_margin_fraction) * margin )
    ROOT.gStyle.SetPadBottomMargin(bottom_margin_fraction * margin)
    ROOT.gStyle.SetPadTopMargin( (1 - bottom_margin_fraction) * margin)
    canvas = ROOT.TCanvas(get_rand_string(), "", 1200, 1200)
    return canvas

def draw_vertical_mass_lines(maxy):
    lines = {}
    for p in particles:
        lines[p] = ROOT.TLine(p.mass2, 0., p.mass2, maxy)
        lines[p].SetLineColor(15)
        lines[p].SetLineWidth(2)
        lines[p].SetLineStyle(1)
        lines[p].Draw()
    return lines





##################### RMS90 ###################

def interval_quantile_(x, quant=0.9):
    '''Calculate the shortest interval that contains the desired quantile'''
    # the number of possible starting points
    n_low = int(len(x) * (1. - quant))
    # the number of events contained in the quantil
    n_quant = len(x) - n_low

    # Calculate all distances in one go
    distances = x[-n_low:] - x[:n_low]
    i_start = np.argmin(distances)

    return i_start, i_start + n_quant

def fit90(x):
    x = np.sort(x)
    n10percent = int(round(len(x)*0.1))
    n90percent = len(x) - n10percent

    start, end = interval_quantile_(x, quant=0.9)

    rms90 = np.std(x[start:end])
    mean90 = np.mean(x[start:end])
    mean90_err = rms90/np.sqrt(n90percent)
    rms90_err = rms90/np.sqrt(2*n90percent)   # estimator in root
    return mean90, rms90, mean90_err, rms90_err



### SEP POWER UTILS

def convert_p_value_to_sep_power(p_value):
    '''Return separation power equivalent to the given p-value (1 - efficiency)'''
    return 2*ROOT.Math.gaussian_quantile_c(p_value, 1)

def convert_sep_power_to_eff(sep_power):
    '''Return efficiency from the given sep power'''
    return ROOT.Math.gaussian_cdf(0.5*sep_power)

def find_optimal_cut(h1, h2, debug=False):
    '''Calcualte the x value between the two histograms where eff=1-misid'''
    cut_mass = []
    efficiencies = []
    mis_ids = []
    diff = []
    n_signal = h1.GetEntries()
    n_background = h2.GetEntries()

    for bin in range(1, h1.GetXaxis().GetNbins() + 1):
        bin_x = h1.GetXaxis().GetBinUpEdge(bin)
        # integral includes first/last bins
        eff = h1.Integral(0, bin)/n_signal if n_signal !=0 else 0
        mis_id = h2.Integral(0, bin)/n_background if n_background !=0 else 0
        if mis_id > eff:
            # h1 is the right histogram and h2 is the left! Flip the effieincy!
            eff, mis_id = 1-eff, 1-mis_id

        cut_mass.append( bin_x )
        efficiencies.append(eff)
        mis_ids.append(mis_id)
        diff.append( abs(1-eff - mis_id) )

    optimal_idx = np.argmin( diff )

    optimal_cut = cut_mass[optimal_idx]
    optimal_eff = efficiencies[optimal_idx]
    optimal_mis_id = mis_ids[optimal_idx]

    # print(f"Cut: {round(optimal_cut, 2)} / eff: {round(optimal_eff, 2)} / mis id: {round(optimal_mis_id, 2)}")
    if debug:
        draw_optimal_cut(h1, h2, optimal_cut)

    return optimal_cut, optimal_eff, optimal_mis_id

def draw_optimal_cut(h1, h2, cut):
    '''Draw two histograms and the calculated cut wehere eff=1-misid. Used for debugging'''
    # my default pi/k/p colours
    h1, h2 = h1.Clone(), h2.Clone()
    h1.Scale(1./h1.GetEntries())
    h2.Scale(1./h2.GetEntries())

    h1_fill = h1.Clone()
    h2_fill = h2.Clone()
    for bin in range(1, h1_fill.GetXaxis().GetNbins() + 1):
        if h1_fill.GetBinCenter(bin) < cut:
            h1_fill.SetBinContent(bin, 0.)
        else:
            h2_fill.SetBinContent(bin, 0.)

    margin = 0.22
    ROOT.gStyle.SetPadLeftMargin(0.9*margin)
    ROOT.gStyle.SetPadRightMargin(0.1*margin)
    ROOT.gStyle.SetPadTopMargin(0.3*margin)
    ROOT.gStyle.SetPadBottomMargin(0.7*margin)
    canvas = ROOT.TCanvas(get_rand_string(),"", 600, 600)

    h1.Draw("hist")
    h1.GetXaxis().SetTitle("Mass^{2} (GeV^{2}/c^{4})")
    h1.GetYaxis().SetTitle("Normalised N entries")
    # h1.GetYaxis().SetTitleOffset(1.2)
    h1.GetYaxis().SetMaxDigits(3)
    h2.Draw("hist same")
    h1_fill.Draw("hist f same")
    h2_fill.Draw("hist f same")

    ymax = 1.05*max(h1.GetMaximum(), h2.GetMaximum())
    h1.SetMaximum(ymax)

    h1.SetLineWidth(4)
    h1.SetLineColor(pion.color)
    h1_fill.SetFillStyle(3001)
    h1_fill.SetFillColor(pion.color)
    h1_fill.SetLineColor(pion.color)

    h2.SetLineWidth(4)
    h2.SetLineColor(kaon.color)
    h2_fill.SetFillStyle(3001)
    h2_fill.SetFillColor(kaon.color)
    h2_fill.SetLineColor(kaon.color)

    line = ROOT.TLine(cut, 0., cut, ymax)
    line.SetLineColor(15)
    line.SetLineWidth(2)
    line.SetLineStyle(9)
    line.Draw()

    canvas.Update()
    input("wait")

def calculate_pid_graphs(h1, h2):
    '''Return graphs of the separation power, efficiency, mis-id, cut values versus momentum based on two 2D histograms'''
    gr_cut = ROOT.TGraph()
    gr_eff = ROOT.TGraph()
    gr_misid = ROOT.TGraph()
    gr_sep_power = ROOT.TGraph()

    for i in range(1, h1.GetXaxis().GetNbins() + 1):
        x_low, x, x_up = h1.GetXaxis().GetBinLowEdge(i), h1.GetXaxis().GetBinCenter(i), h1.GetXaxis().GetBinUpEdge(i)
        h1_proj = h1.ProjectionY("h1_proj", i, i)
        h2_proj = h2.ProjectionY("h2_proj", i, i)

        # print(f"Momentum: {x_low:.2f} -- {x_up:.2f}")
        cut, eff, misid = find_optimal_cut(h1_proj, h2_proj)
        sep_power = convert_p_value_to_sep_power(1-eff)

        gr_cut.SetPoint(i-1, x, cut)
        gr_eff.SetPoint(i-1, x, eff)
        gr_misid.SetPoint(i-1, x, misid)
        gr_sep_power.SetPoint(i-1, x, sep_power)
    return gr_cut, gr_eff, gr_misid, gr_sep_power

def combine_two_graphs(gr1, gr2):
    '''Return TGraph as a sum in quadratures of gr1 and gr2 point by point'''
    n_points = gr1.GetN()
    if n_points is not gr2.GetN():
        raise ValueError(f'Two graphs must have equal number of points! First has {gr1.GetN()} and second has {gr2.GetN()}')

    gr = ROOT.TGraph()
    for i in range( n_points ):
        x = gr1.GetPointX(i)
        y1 = gr1.GetPointY(i)
        y2 = gr2.GetPointY(i)
        y = np.sqrt(y1*y1 + y2*y2)
        gr.SetPoint(i, x, y)
    return gr

def convert_graph_sp_to_eff(gr_sp):
    gr_eff = ROOT.TGraph()
    for i in range( gr_sp.GetN() ):
        x = gr_sp.GetPointX(i)
        sp = gr_sp.GetPointY(i)
        eff = convert_sep_power_to_eff(sp)
        gr_eff.SetPoint(i, x, eff)
    return gr_eff

def draw_resolution_sep_powers(graphs):
    '''Draw the separation power graphs for different TOF resolutions'''
    canvas = create_canvas(0.38, 0.4, 0.65)
    canvas.SetTicky(False)

    legend = create_legend(0.4, 0.63, 0.98, 0.94)
    for i, (res, gr) in enumerate(graphs.items()):
        if i == 0:
            # canvas.DrawFrame(0., MIN_SEP_POWER, 19., MAX_SEP_POWER)
            x_last = MAX_MOMENTUM

            gr.Draw("ALP")
            gr.GetYaxis().SetTitleOffset(0.9)
            gr.GetXaxis().SetRangeUser(0, x_last)
            gr.GetXaxis().SetNdivisions(512)
            gr.GetYaxis().SetRangeUser(MIN_SEP_POWER, MAX_SEP_POWER)
            canvas.Modified()
            canvas.Update()
            # draw an axis on the right side
            x_pos = canvas.GetUxmax()
            axis_eff = ROOT.TGaxis(x_pos, canvas.GetUymin(), x_pos-0.001, canvas.GetUymax(), MIN_SEP_POWER, MAX_SEP_POWER)
            axis_eff.SetTitleColor( ROOT.gStyle.GetTitleColor("Y") )
            axis_eff.SetTitleFont( ROOT.gStyle.GetTitleFont("Y") )
            axis_eff.SetTitleSize( ROOT.gStyle.GetTitleSize("Y") )
            axis_eff.CenterTitle(True)
            axis_eff.SetTitleOffset(2.2)
            axis_eff.SetTitle("Efficiency (%)")
            axis_eff.SetLabelColor( ROOT.gStyle.GetLabelColor("Y") )
            axis_eff.SetLabelFont( ROOT.gStyle.GetLabelFont("Y") )
            axis_eff.SetLabelOffset(-0.14)
            axis_eff.SetLabelSize( ROOT.gStyle.GetLabelSize("Y") )
            axis_eff.SetTickLength(0.03)
            for j in range(int(MIN_SEP_POWER), int(MAX_SEP_POWER)+1):
                axis_eff.ChangeLabel(j+1, -1, -1, -1, ROOT.kBlack, -1, f"{convert_sep_power_to_eff(j)*100:.2f}") # only in the new ROOT versions
            axis_eff.DrawClone()
        else:
            gr.Draw("LPsame")

        gr.SetLineColor(COLORS_RESOLUTION[i])
        gr.SetMarkerColor(COLORS_RESOLUTION[i])
        gr.SetLineWidth(4)
        gr.SetMarkerStyle(20)
        legend.AddEntry(gr, f"{res} ps","lp")

    legend.DrawClone()

    # latex = ROOT.TLatex()
    # latex.SetTextFont(52)
    # latex.DrawLatex(12, 6.06, "ILD preliminary")

    canvas.Update()
    return canvas

def draw_dedx_sep_powers(gr_tof, gr_dedx):
    '''Draw the separation power graphs for TOF and dEdx'''
    COLORS_DEDX = [ ROOT.TColor.GetColor(c) for c in ["#00aaff", "#cd54b9", "#c52700"] ]
    MIN_SEP_POWER, MAX_SEP_POWER = -0.2, 7.5

    canvas = create_canvas(0.38, 0.4, 0.65)
    canvas.SetLogx()
    canvas.SetTicky(False)

    legend = create_legend(0.35, 0.68, 0.74, 0.84)
    gr_tof.Draw("ALP")
    gr_tof.GetXaxis().SetRangeUser(0.5, 20)
    gr_tof.GetXaxis().SetMoreLogLabels()
    gr_tof.GetXaxis().SetNoExponent()
    gr_tof.GetXaxis().SetNdivisions(506)
    gr_tof.GetYaxis().SetTitleOffset(0.9)
    gr_tof.GetYaxis().SetRangeUser(MIN_SEP_POWER, MAX_SEP_POWER)
    gr_tof.SetLineColor(COLORS_DEDX[0])
    gr_tof.SetMarkerColor(COLORS_DEDX[0])
    gr_tof.SetLineWidth(4)
    gr_tof.SetMarkerStyle(20)
    legend.AddEntry(gr_tof, gr_tof.GetTitle(),"lp")
    canvas.Modified()
    canvas.Update()
    # draw an axis on the right side
    x_pos = 10**canvas.GetUxmax() # NOTE: FREAKING ROOT https://root-forum.cern.ch/t/getumin-getumax-show-wrong-results-for-the-canvases-with-the-log-scale/58867
    axis_eff = ROOT.TGaxis(x_pos, canvas.GetUymin(), x_pos-0.001, canvas.GetUymax(), MIN_SEP_POWER, MAX_SEP_POWER)
    axis_eff.SetTitleColor( ROOT.gStyle.GetTitleColor("Y") )
    axis_eff.SetTitleFont( ROOT.gStyle.GetTitleFont("Y") )
    axis_eff.SetTitleSize( ROOT.gStyle.GetTitleSize("Y") )
    axis_eff.CenterTitle(True)
    axis_eff.SetTitleOffset(2.2)
    axis_eff.SetTitle("Efficiency (%)")
    axis_eff.SetLabelColor( ROOT.gStyle.GetLabelColor("Y") )
    axis_eff.SetLabelFont( ROOT.gStyle.GetLabelFont("Y") )
    axis_eff.SetLabelOffset(-0.14)
    axis_eff.SetLabelSize( ROOT.gStyle.GetLabelSize("Y") )
    axis_eff.SetTickLength(0.03)
    for j in range(int(MIN_SEP_POWER), int(MAX_SEP_POWER)+1):
        axis_eff.ChangeLabel(j+1, -1, -1, -1, ROOT.kBlack, -1, f"{convert_sep_power_to_eff(j)*100:.2f}")
    axis_eff.DrawClone()
    gr_dedx.Draw("LPsame")
    gr_dedx.SetLineColor(COLORS_DEDX[2])
    gr_dedx.SetMarkerColor(COLORS_DEDX[2])
    gr_dedx.SetLineWidth(4)
    gr_dedx.SetMarkerStyle(20)
    legend.AddEntry(gr_dedx, gr_dedx.GetTitle(),"lp")

    gr_combined = ROOT.TGraph()
    gr_combined.SetTitle("dE/dx #oplus TOF")
    for i in range( gr_tof.GetN() ):
        x = gr_dedx.GetPointX(i)
        sp_dedx = gr_dedx.GetPointY(i)
        sp_tof = gr_tof.GetPointY(i)
        gr_combined.SetPoint(i, x, np.sqrt(sp_dedx*sp_dedx + sp_tof*sp_tof))

    gr_combined.SetLineColor(COLORS_DEDX[1])
    gr_combined.SetMarkerColor(COLORS_DEDX[1])
    gr_combined.SetLineWidth(4)
    gr_combined.SetMarkerStyle(20)
    legend.AddEntry(gr_combined, gr_combined.GetTitle(),"lp")
    gr_combined.DrawClone("LPsame")

    legend.DrawClone()

    # latex = ROOT.TLatex()
    # latex.SetTextFont(52)
    # latex.DrawLatex(12, 6.06, "ILD preliminary")

    canvas.Update()
    return canvas
