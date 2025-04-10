#!/usr/bin/env python3

"""Produces momentum distribution plots for chapter 8"""

import numpy as np
import ROOT
from utils import *

ROOT.EnableImplicitMT()

dark24 = ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038']
dark24 = [ROOT.TColor.GetColor(c) for c in dark24]

#Linear momentum bins
N_MOMENTUM_BINS, MIN_MOMENTUM, MAX_MOMENTUM = 70, 0, 20 # GeV/c
#Log momentum bins for dE/dx
N_LOG_MOMENTUM_BINS, MIN_LOG_MOMENTUM, MAX_LOG_MOMENTUM = 30, -0.3, 1.3 # (0.5 to ~20) GeV/c
# LOG_MOMENTUM_BINS = np.array([ 10**(MIN_LOG_MOMENTUM + (MAX_LOG_MOMENTUM-MIN_LOG_MOMENTUM)*i/N_LOG_MOMENTUM_BINS ) for i in range(N_LOG_MOMENTUM_BINS+1) ])
LOG_MOMENTUM_BINS = np.logspace(MIN_LOG_MOMENTUM, MAX_LOG_MOMENTUM, N_LOG_MOMENTUM_BINS+1)
N_DEDX_BINS, MIN_DEDX, MAX_DEDX = 3000, 0, 1e-6
DEDX_BINS = np.linspace(MIN_DEDX, MAX_DEDX, N_DEDX_BINS+1)
N_MASS2_BINS, MIN_MASS2, MAX_MASS2 = 3000, -3, 3  # GeV^2/c^4
# MASS2_BINS = np.array([ MIN_MASS2 + (MAX_MASS2 - MIN_MASS2)*i/N_MASS2_BINS for i in range(N_MASS2_BINS+1) ])
MASS2_BINS = np.linspace(MIN_MASS2, MAX_MASS2, N_MASS2_BINS+1)
MOMENTUM_COLUMN = "harmonicMomToEcal_IKF_zedLambda" 
TRACK_LENGTH_COLUMN = "trackLengthToEcal_IKF_zedLambda"
RESOLUTIONS = [0, 1, 5, 10, 30, 50, 100, 300] # ps
COLORS_RESOLUTION = [ ROOT.TColor.GetColor(c) for c in ["#00aaff", "#0091ea", "#0079d3", "#0061bd", "#004aa5", "#00348d", "#001d75", "#00045c"] ]
COLORS_DEDX = [ ROOT.TColor.GetColor(c) for c in ["#00aaff", "#cd54b9", "#c52700"] ]
MIN_SEP_POWER, MAX_SEP_POWER = 0, 6






def test():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root")
    df = df.Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")\
            .Define("pt", "sqrt(mcPx*mcPx + mcPy*mcPy)")\
            .Define("costheta", "cos(atan2(pt, mcPz))")
    df_gen = df.Filter("!isSimulated && !isOverlay").Filter("abs(pdg) == 2212")
    df_reco = df_gen.Filter("hasTrack")

    # h_gen_test = df_gen.Histo2D((get_rand_string(), "Generated;pt;pz", 1000, 0, 4, 1000, -2, 2), "pt", "mcPz" )
    # h_reco_test = df_reco.Histo2D((get_rand_string(), "Has track;pt;pz", 1000, 0, 4, 1000, -2, 2), "pt", "mcPz" )
    ROOT.gStyle.SetPalette(1)
    c1 = create_canvas()
    h_gen_test = df_gen.Histo2D((get_rand_string(), "Generated;cos theta;pt", 200,-1, 1, 200, 0, 4), "costheta", "pt" )
    h_reco_test = df_reco.Histo2D((get_rand_string(), "Has track;cos theta;pt", 200,-1, 1, 200, 0, 4), "costheta", "pt" )
    h_reco_test.Divide(h_gen_test.GetPtr())
    h_reco_test.Draw("colz")
    c1.Update()
    c2 = create_canvas()
    h_gen_test2 = df_gen.Histo2D((get_rand_string(), "Generated;cos theta;p", 200,-1, 1, 200, 0, 4), "costheta", "mom" )
    h_reco_test2 = df_reco.Histo2D((get_rand_string(), "Has track;cos theta;p", 200,-1, 1, 200, 0, 4), "costheta", "mom" )
    h_reco_test2.Divide(h_gen_test2.GetPtr())
    h_reco_test2.Draw("colz")
    c2.Update()

    input("wait")

def draw_hadr_vs_secondary():
    # df analysis
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root")
    df = df.Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    df_gen = df.Filter("!isSimulated && !isOverlay").Filter("quarksToPythia == 55").Filter("abs(pdg) == 321")
    df1 = df_gen.Filter("isHadronisationDecay")
    df2 = df_gen.Filter("isBottomQuarkDecay")

    #creation and styling
    h1 = df_gen.Histo1D((get_rand_string(), "Generated;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    h2 = df1.Histo1D((get_rand_string(), "Not from b-quark decay;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    h3 = df2.Histo1D((get_rand_string(), "From b quark decay;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    n_gen_low = df_gen.Filter("mom < 3").Count()
    n_reco_low = df1.Filter("mom < 3").Count()
    n_reco_shower_low = df2.Filter("mom < 3").Count()
    n_gen_total = df_gen.Count()
 
    h2.SetLineColor(ROOT.TColor.GetColor("#20bf55"))
    h3.SetLineColor(ROOT.TColor.GetColor("#01baef"))
    h1.SetLineWidth(4)
    h2.SetLineWidth(4)
    h3.SetLineWidth(4)

    #Calculating
    print(f"Fraction of generated below 3 GeV {100*n_gen_low.GetValue()/n_gen_total.GetValue()}")
    print(f"Fraction of track below 3 GeV {100*n_reco_low.GetValue()/n_gen_total.GetValue()}")
    print(f"Fraction of track+shower below 3 GeV {100*n_reco_shower_low.GetValue()/n_gen_total.GetValue()}")


    # plotting
    c1 = create_canvas()
    h1.GetYaxis().SetTitleOffset(1.3)
    h1.GetYaxis().SetMaxDigits(3)
    h1.GetYaxis().SetNdivisions(512)
    h1.DrawClone()
    h2.DrawClone("same")
    h3.DrawClone("same")
    leg = create_legend()
    leg.AddEntry(h1.GetPtr(), h1.GetTitle(), "l")
    leg.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    leg.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    leg.Draw()
    c1.Update()

    c2 = create_canvas()
    h2.Divide(h1.GetPtr())
    h3.Divide(h1.GetPtr())
    h2.Draw("histo")
    h2.GetYaxis().SetTitleOffset(1.3)
    h2.GetYaxis().SetMaxDigits(3)
    h2.GetYaxis().SetNdivisions(512)
    h3.Draw("histo same")
    leg2 = create_legend()
    leg2.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    leg2.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    leg2.Draw()
    c2.Update()
    #splitline{Kaons from Z#rightarrowb#bar{b}}{E_{cm} = 250 GeV/c^{2}}

    input("wait")


def draw_prim_vs_secondary():
    # df analysis
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root")
    df = df.Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    df_gen = df.Filter("!isSimulated && !isOverlay").Filter("quarksToPythia == 55").Filter("abs(pdg) == 321").Filter("hasTrack")
    df1 = df_gen.Filter("isInRecoPrimaryVertex")
    df2 = df_gen.Filter("isInRecoSecondaryVertex")

    #creation and styling
    h1 = df_gen.Histo1D((get_rand_string(), "Total reconstructed (has track);Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    h2 = df1.Histo1D((get_rand_string(), "Reconstructed in primary vertex;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    h3 = df2.Histo1D((get_rand_string(), "Reconstructed in secondary vertex;Momentum (GeV/c);N entries", 200, 0, 10), "mom" )
    n_gen_low = df_gen.Filter("mom < 3").Count()
    n_reco_low = df1.Filter("mom < 3").Count()
    n_reco_shower_low = df2.Filter("mom < 3").Count()
    n_gen_total = df_gen.Count()
 
    h2.SetLineColor(ROOT.TColor.GetColor("#20bf55"))
    h3.SetLineColor(ROOT.TColor.GetColor("#01baef"))
    h1.SetLineWidth(4)
    h2.SetLineWidth(4)
    h3.SetLineWidth(4)

    #Calculating
    print(f"Fraction of generated below 3 GeV {100*n_gen_low.GetValue()/n_gen_total.GetValue()}")
    print(f"Fraction of track below 3 GeV {100*n_reco_low.GetValue()/n_gen_total.GetValue()}")
    print(f"Fraction of track+shower below 3 GeV {100*n_reco_shower_low.GetValue()/n_gen_total.GetValue()}")


    # plotting
    c1 = create_canvas()
    h1.GetYaxis().SetTitleOffset(1.3)
    h1.GetYaxis().SetMaxDigits(3)
    h1.GetYaxis().SetNdivisions(512)
    h1.DrawClone()
    h2.DrawClone("same")
    h3.DrawClone("same")
    leg = create_legend()
    leg.AddEntry(h1.GetPtr(), h1.GetTitle(), "l")
    leg.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    leg.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    leg.Draw()
    c1.Update()

    c2 = create_canvas()
    h2.Divide(h1.GetPtr())
    h3.Divide(h1.GetPtr())
    h2.Draw("histo")
    h2.GetYaxis().SetTitleOffset(1.3)
    h2.GetYaxis().SetMaxDigits(3)
    h2.GetYaxis().SetNdivisions(512)
    h3.Draw("histo same")
    leg2 = create_legend()
    leg2.AddEntry(h2.GetPtr(), h2.GetTitle(), "l")
    leg2.AddEntry(h3.GetPtr(), h3.GetTitle(), "l")
    leg2.Draw()
    c2.Update()
    #splitline{Kaons from Z#rightarrowb#bar{b}}{E_{cm} = 250 GeV/c^{2}}

    input("wait")


def get_particles_fractions():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_hadronic_dst.root").Filter("!isSimulated && !isOverlay").Filter("quarksToPythia == 55")
    leg = create_legend()
    n_total = df.Count()
    n_pi = df.Filter("abs(pdg) == 211").Count()
    n_k = df.Filter("abs(pdg) == 321").Count()
    n_p = df.Filter("abs(pdg) == 2212").Count()

    print("Pi fraction:", 100*n_pi.GetValue()/n_total.GetValue())
    print("K fraction:", 100*n_k.GetValue()/n_total.GetValue())
    print("P fraction:", 100*n_p.GetValue()/n_total.GetValue())
    input("wait")

# get_particles_fractions()
# draw_gen_vs_reco()
# draw_hadr_vs_secondary()
# test()
c1 = get_efficiency_graphs()

# c2 = get_efficiency_graphs(tof_only=True)
input("wait")