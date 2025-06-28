import ROOT
import numpy as np
from utils import *
ROOT.EnableImplicitMT()

dark24 = ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038']

colors = [ROOT.TColor.GetColor(c) for c in dark24]

def general():
    histos = {}

    n_events_z = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root").Count()
    n_events_ww = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/4f_WW_semileptonic_eLpR_dst.root").Count()
    n_events_higgs = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h.root").Count()

    df_z = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 321")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    df_ww = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/4f_WW_semileptonic_eLpR_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 321")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    df_higgs = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 321")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    histos["Zqq"] = df_z.Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["WWqq"] = df_ww.Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["H"] = df_higgs.Histo1D((get_rand_string(), "", 500, 0, 10), "mom")

    of_z = df_z.Filter("mom > 10").Count()
    of_ww = df_ww.Filter("mom > 10").Count()
    of_h = df_higgs.Filter("mom > 10").Count()

    z_overflow = 100.*of_z.GetValue()/histos["Zqq"].GetEntries()
    ww_overflow = 100.*of_ww.GetValue()/histos["WWqq"].GetEntries()
    h_overflow = 100.*of_h.GetValue()/histos["H"].GetEntries()

    histos["Zqq"].Scale(1./n_events_z.GetValue())
    histos["WWqq"].Scale(1./n_events_ww.GetValue())
    histos["H"].Scale(1./n_events_higgs.GetValue())

    histos["Zqq"].SetTitle("Z #rightarrow q#bar{q}" + f" (overflow: {z_overflow:.1f} %)" + ";Momentum (GeV/c);N kaons per event")
    histos["WWqq"].SetTitle("WW #rightarrow l#nu_{l} q#bar{q}" + f" (overflow: {ww_overflow:.1f} %)")
    histos["H"].SetTitle("ZH #rightarrow #nu_{#mu,#tau} #bar{#nu}_{#mu,#tau} H" + f" (overflow: {h_overflow:.1f} %)")

    c = create_canvas()
    for i, h in enumerate(histos.values()):
        h.SetLineColor(colors[i])
        h.Draw("histo" if i==0 else "histo same")

    histos["Zqq"].SetMaximum(1.05*max([h.GetMaximum() for h in histos.values()]))
    histos["Zqq"].GetYaxis().SetMaxDigits(3)
    histos["Zqq"].GetXaxis().SetNdivisions(512)
    leg = c.BuildLegend(0.2, 0.75, 0.76, 0.91, "", "l")
    leg.SetFillStyle(0)
    leg.SetMargin(0.15)
    c.Update()
    input("wait")


# general()

def zqq():
    histos = {}

    df_events = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 321")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    n_events_bb = df_events.Filter("quarksToPythia == 55").Count()
    n_events_cc = df_events.Filter("quarksToPythia == 44").Count()
    n_events_ss = df_events.Filter("quarksToPythia == 33").Count()
    n_events_uu = df_events.Filter("quarksToPythia == 22").Count()
    n_events_dd = df_events.Filter("quarksToPythia == 11").Count()

    histos["Zbb"] = df.Filter("quarksToPythia == 55").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Zcc"] = df.Filter("quarksToPythia == 44").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Zss"] = df.Filter("quarksToPythia == 33").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Zuu"] = df.Filter("quarksToPythia == 22").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Zdd"] = df.Filter("quarksToPythia == 11").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")

    of_z_bb = df.Filter("quarksToPythia == 55").Filter("mom > 10").Count()
    of_z_cc = df.Filter("quarksToPythia == 44").Filter("mom > 10").Count()
    of_z_ss = df.Filter("quarksToPythia == 33").Filter("mom > 10").Count()
    of_z_uu = df.Filter("quarksToPythia == 22").Filter("mom > 10").Count()
    of_z_dd = df.Filter("quarksToPythia == 11").Filter("mom > 10").Count()

    bb_overflow = 100.*of_z_bb.GetValue()/histos["Zbb"].GetEntries()
    cc_overflow = 100.*of_z_cc.GetValue()/histos["Zcc"].GetEntries()
    ss_overflow = 100.*of_z_ss.GetValue()/histos["Zss"].GetEntries()
    uu_overflow = 100.*of_z_uu.GetValue()/histos["Zuu"].GetEntries()
    dd_overflow = 100.*of_z_dd.GetValue()/histos["Zdd"].GetEntries()

    histos["Zbb"].Scale(1./n_events_bb.GetValue())
    histos["Zcc"].Scale(1./n_events_cc.GetValue())
    histos["Zss"].Scale(1./n_events_ss.GetValue())
    histos["Zuu"].Scale(1./n_events_uu.GetValue())
    histos["Zdd"].Scale(1./n_events_dd.GetValue())

    histos["Zbb"].SetTitle("Z #rightarrow b#bar{b}" + f" (overflow: {bb_overflow:.1f} %)" + ";Momentum (GeV/c);N kaons per event")
    histos["Zcc"].SetTitle("Z #rightarrow c#bar{c}" + f" (overflow: {cc_overflow:.1f} %)")
    histos["Zss"].SetTitle("Z #rightarrow s#bar{s}" + f" (overflow: {ss_overflow:.1f} %)")
    histos["Zuu"].SetTitle("Z #rightarrow u#bar{u}" + f" (overflow: {uu_overflow:.1f} %)")
    histos["Zdd"].SetTitle("Z #rightarrow d#bar{d}" + f" (overflow: {dd_overflow:.1f} %)")

    c = create_canvas()
    for i, h in enumerate(histos.values()):
        h.SetLineColor(colors[i])
        h.Draw("histo" if i==0 else "histo same")


    histos["Zbb"].SetMaximum(1.05*max([h.GetMaximum() for h in histos.values()]))
    histos["Zbb"].GetYaxis().SetMaxDigits(3)
    histos["Zbb"].GetXaxis().SetNdivisions(512)
    leg = c.BuildLegend(0.2, 0.75, 0.76, 0.91, "", "l")
    leg.SetMargin(0.15)
    leg.SetFillStyle(0)
    c.Update()
    input("wait")

# zqq()
    
def ww():
    histos = {}

    df_events = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/4f_WW_semileptonic_eLpR_dst.root")
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/4f_WW_semileptonic_eLpR_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 2212")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    n_events_ud = df_events.Filter("quarksToPythia == 12 || quarksToPythia == 21").Count()
    n_events_us = df_events.Filter("quarksToPythia == 32 || quarksToPythia == 23").Count()
    n_events_cd = df_events.Filter("quarksToPythia == 14 || quarksToPythia == 41").Count()
    n_events_cs = df_events.Filter("quarksToPythia == 34 || quarksToPythia == 43").Count()

    histos["WWud"] = df.Filter("quarksToPythia == 12 || quarksToPythia == 21").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["WWus"] = df.Filter("quarksToPythia == 32 || quarksToPythia == 23").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["WWcd"] = df.Filter("quarksToPythia == 14 || quarksToPythia == 41").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["WWcs"] = df.Filter("quarksToPythia == 34 || quarksToPythia == 43").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")

    of_ww_ud = df.Filter("quarksToPythia == 12 || quarksToPythia == 21").Filter("mom > 10").Count()
    of_ww_us = df.Filter("quarksToPythia == 32 || quarksToPythia == 23").Filter("mom > 10").Count()
    of_ww_cd = df.Filter("quarksToPythia == 14 || quarksToPythia == 41").Filter("mom > 10").Count()
    of_ww_cs = df.Filter("quarksToPythia == 34 || quarksToPythia == 43").Filter("mom > 10").Count()


    ud_overflow = 100.*of_ww_ud.GetValue()/histos["WWud"].GetEntries()
    us_overflow = 100.*of_ww_us.GetValue()/histos["WWus"].GetEntries()
    cd_overflow = 100.*of_ww_cd.GetValue()/histos["WWcd"].GetEntries()
    cs_overflow = 100.*of_ww_cs.GetValue()/histos["WWcs"].GetEntries()

    histos["WWud"].Scale(1./n_events_ud.GetValue())
    histos["WWus"].Scale(1./n_events_us.GetValue())
    histos["WWcd"].Scale(1./n_events_cd.GetValue())
    histos["WWcs"].Scale(1./n_events_cs.GetValue())

    histos["WWud"].SetTitle("WW #rightarrow l#nu_{l} ud" + f" (overflow: {ud_overflow:.1f} %)" + ";Momentum (GeV/c);N protons per event")
    histos["WWus"].SetTitle("WW #rightarrow l#nu_{l} us" + f" (overflow: {us_overflow:.1f} %)")
    histos["WWcd"].SetTitle("WW #rightarrow l#nu_{l} cd" + f" (overflow: {cd_overflow:.1f} %)")
    histos["WWcs"].SetTitle("WW #rightarrow l#nu_{l} cs" + f" (overflow: {cs_overflow:.1f} %)")

    c = create_canvas()
    for i, h in enumerate(histos.values()):
        h.SetLineColor(colors[i])
        h.Draw("histo" if i==0 else "histo same")


    histos["WWud"].SetMaximum(1.05*max([h.GetMaximum() for h in histos.values()]))
    histos["WWud"].GetYaxis().SetMaxDigits(3)
    histos["WWud"].GetXaxis().SetNdivisions(512)
    leg = c.BuildLegend(0.2, 0.75, 0.76, 0.91, "", "l")
    leg.SetMargin(0.15)
    leg.SetFillStyle(0)
    c.Update()
    input("wait")
# ww()ll 


def higgs():
    histos = {}

    df_events = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h_dst.root")
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 2212")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    df_events_ss = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h_ss_dst.root")
    df_ss = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h_ss_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 2212")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    n_events_bb = df_events.Filter("higgsDaughters == 505").Count()
    n_events_ww = df_events.Filter("higgsDaughters == 2424").Count()
    n_events_gg = df_events.Filter("higgsDaughters == 2121 || higgsDaughters == 909").Count()
    n_events_cc = df_events.Filter("higgsDaughters == 404").Count()
    n_events_zz = df_events.Filter("higgsDaughters == 2323").Count()
    n_events_ss = df_events_ss.Filter("higgsDaughters == 303").Count()

    histos["Hbb"] = df.Filter("higgsDaughters == 505").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Hww"] = df.Filter("higgsDaughters == 2424").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Hgg"] = df.Filter("higgsDaughters == 2121 || higgsDaughters == 909").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Hcc"] = df.Filter("higgsDaughters == 404").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Hzz"] = df.Filter("higgsDaughters == 2323").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")
    histos["Hss"] = df_ss.Filter("higgsDaughters == 303").Histo1D((get_rand_string(), "", 500, 0, 10), "mom")

    #overflows
    of_hbb = df.Filter("higgsDaughters == 505").Filter("mom > 10").Count()
    of_hww = df.Filter("higgsDaughters == 2424").Filter("mom > 10").Count()
    of_hgg = df.Filter("higgsDaughters == 2121 || higgsDaughters == 909").Filter("mom > 10").Count()
    of_hcc = df.Filter("higgsDaughters == 404").Filter("mom > 10").Count()
    of_hzz = df.Filter("higgsDaughters == 2323").Filter("mom > 10").Count()
    of_hss = df_ss.Filter("higgsDaughters == 303").Filter("mom > 10").Count()

    hbb_overflow = 100.*of_hbb.GetValue()/histos["Hbb"].GetEntries()
    hww_overflow = 100.*of_hww.GetValue()/histos["Hww"].GetEntries()
    hgg_overflow = 100.*of_hgg.GetValue()/histos["Hgg"].GetEntries()
    hcc_overflow = 100.*of_hcc.GetValue()/histos["Hcc"].GetEntries()
    hzz_overflow = 100.*of_hzz.GetValue()/histos["Hzz"].GetEntries()
    hss_overflow = 100.*of_hss.GetValue()/histos["Hss"].GetEntries()

    histos["Hbb"].Scale(1./n_events_bb.GetValue())
    histos["Hww"].Scale(1./n_events_ww.GetValue())
    histos["Hgg"].Scale(1./n_events_gg.GetValue())
    histos["Hcc"].Scale(1./n_events_cc.GetValue())
    histos["Hzz"].Scale(1./n_events_zz.GetValue())
    histos["Hss"].Scale(1./n_events_ss.GetValue())

    histos["Hbb"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} b#bar{b}" + f" (overflow: {hbb_overflow:.1f} %)" + ";Momentum (GeV/c);N protons per event")
    histos["Hww"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} WW" + f" (overflow: {hww_overflow:.1f} %)")
    histos["Hgg"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} gg" + f" (overflow: {hgg_overflow:.1f} %)")
    histos["Hcc"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} c#bar{c}" + f" (overflow: {hcc_overflow:.1f} %)")
    histos["Hzz"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} ZZ" + f" (overflow: {hzz_overflow:.1f} %)")
    histos["Hss"].SetTitle("ZH #rightarrow #nu_{#mu,#tau}#bar{#nu}_{#mu,#tau} s#bar{s}" + f" (overflow: {hss_overflow:.1f} %)")

    c = create_canvas()
    for i, h in enumerate(histos.values()):
        h.SetLineColor(colors[i])
        h.Draw("histo" if i==0 else "histo same")


    histos["Hbb"].SetMaximum(1.05*max([h.GetMaximum() for h in histos.values()]))
    histos["Hbb"].GetYaxis().SetMaxDigits(3)
    histos["Hbb"].GetXaxis().SetNdivisions(512)
    leg = c.BuildLegend(0.2, 0.75, 0.76, 0.91, "", "l")
    leg.SetMargin(0.15)
    leg.SetFillStyle(0)
    c.Update()
    input("wait")
# higgs()

def overlay():
    histos = {}

    df_events = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")\
                .Filter("!isSimulated && isOverlay")\
                .Filter("abs(pdg) == 2212")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    n_events = df_events.Count()

    histos["overlay"] = df.Histo1D((get_rand_string(), "", 500, 0, 10), "mom")

    #overflows
    of_overlay = df.Filter("mom > 10").Count()

    overlay_overflow = 100.*of_overlay.GetValue()/histos["overlay"].GetEntries()

    histos["overlay"].Scale(1./n_events.GetValue())

    histos["overlay"].SetTitle("#gamma#gamma #rightarrow low p_{T} hadrons" + f" (overflow: {overlay_overflow:.1f} %)" + ";Momentum (GeV/c);N protons per event")

    c = create_canvas()
    for i, h in enumerate(histos.values()):
        h.SetLineColor(colors[i])
        h.Draw("histo" if i==0 else "histo same")


    histos["overlay"].SetMaximum(1.05*max([h.GetMaximum() for h in histos.values()]))
    histos["overlay"].GetYaxis().SetMaxDigits(3)
    histos["overlay"].GetXaxis().SetNdivisions(512)
    leg = c.BuildLegend(0.2, 0.75, 0.76, 0.91, "", "l")
    leg.SetMargin(0.15)
    leg.SetFillStyle(0)
    c.Update()
    input("wait")
# overlay()
    

def reco_impact():
    #NOTE: .Filter("abs(imidiateParentPDG) != 211 && abs(imidiateParentPDG) != 321")\
    histos = {}
    df_events = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")
    n_events = df_events.Filter("quarksToPythia == 33").Count()
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")\
                .Filter("!isSimulated && !isOverlay")\
                .Filter("abs(pdg) == 2212")\
                .Filter("quarksToPythia == 33")\
                .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    histos["gen"] = df.Histo1D((get_rand_string(), "", 500, 0, 10), "mom" )
    histos["track"] = df.Filter("hasTrack").Histo1D((get_rand_string(), "", 500, 0, 10), "mom" )
    histos["nodecay"] = df.Filter("hasTrack && !isDecayedInTracker").Histo1D((get_rand_string(), "", 500, 0, 10), "mom" )
    histos["shower"] = df.Filter("hasTrack && !isDecayedInTracker && hasShower").Histo1D((get_rand_string(), "", 500, 0, 10), "mom" )
 
    histos["gen"].Scale(1./n_events.GetValue())
    histos["track"].Scale(1./n_events.GetValue())
    histos["nodecay"].Scale(1./n_events.GetValue())
    histos["shower"].Scale(1./n_events.GetValue())

    histos["gen"].SetTitle("Total generated;Momentum (GeV/c);N protons per event")
    histos["track"].SetTitle("Have a track")
    histos["nodecay"].SetTitle("Have a track + don't decay in a tracker")
    histos["shower"].SetTitle("Have a track + don't decay in a tracker + have a shower")

    c = create_canvas()
    for i, h in enumerate(histos.values()):
        h.SetLineWidth(4)
        h.SetLineColor(colors[i])
        h.Draw("histo" if i==0 else "histo same")


    histos["gen"].SetMaximum(1.05*max([h.GetMaximum() for h in histos.values()]))
    histos["gen"].GetYaxis().SetMaxDigits(3)
    histos["gen"].GetXaxis().SetNdivisions(512)
    leg = c.BuildLegend(0.28, 0.80, 0.93, 0.92, "", "l")
    leg.SetMargin(0.15)
    leg.SetFillStyle(0)

    latex.SetTextSize(0.04)
    latex.DrawLatexNDC(0.55, 0.55, "#splitline{Protons from Z#rightarrows#bar{s}}{E_{cm} = 250 GeV/c^{2}}")

    c.Update()
    input("wait")
# reco_impact()



def apply_efficiency(h_mom, gr_eff, mode="eff"):
    '''Return histogram where each bin is multiplied by the interpolated graph value (mode=eff) or 1-graph value (mode=mis-id)'''
    h = h_mom.Clone()
    for i in range( 1, h_mom.GetXaxis().GetNbins() + 1 ):
        x = h_mom.GetXaxis().GetBinCenter(i)
        eff = gr_eff.Eval(x)
        if mode == "eff":
            h.SetBinContent(i, h_mom.GetBinContent(i)*eff )
        elif mode == "misid":
            h.SetBinContent(i, h_mom.GetBinContent(i)*(1-eff) )
    return h


def tof_impact(use_dedx=True, process="zss", particle=kaon):
    df_eff = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root").Define("mom", "harmonicMomToEcal_IKF_zedLambda")
    df_dedx = df_eff.Filter("dEdx > 0. && tofClosest0 > 6.")
    df_dedx_pi = df_dedx.Filter("abs(pdg) == 211")
    df_dedx_k = df_dedx.Filter("abs(pdg) == 321")
    df_dedx_p = df_dedx.Filter("abs(pdg) == 2212")

    N_LOG_MOMENTUM_BINS, MIN_LOG_MOMENTUM, MAX_LOG_MOMENTUM = 30, -0.3, 1.3 # (0.5 to ~20) GeV/c
    LOG_MOMENTUM_BINS = np.logspace(MIN_LOG_MOMENTUM, MAX_LOG_MOMENTUM, N_LOG_MOMENTUM_BINS+1)
    N_DEDX_BINS, MIN_DEDX, MAX_DEDX = 3000, 0, 1e-6
    DEDX_BINS = np.linspace(MIN_DEDX, MAX_DEDX, N_DEDX_BINS+1)

    h2d_dedx_pi = df_dedx_pi.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_DEDX_BINS, DEDX_BINS), "mom", "dEdx")
    h2d_dedx_k = df_dedx_k.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_DEDX_BINS, DEDX_BINS), "mom", "dEdx")
    h2d_dedx_p = df_dedx_p.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_DEDX_BINS, DEDX_BINS), "mom", "dEdx")

    df_tof = df_eff.Filter("tofClosest0 > 6.")
    df_tof30 = df_tof.Define("beta", "trackLengthToEcal_IKF_zedLambda/(tofClosest30*299.792458)")\
                .Define("mass2", "mom*mom*( 1./(beta*beta) - 1.)")
    df_tof30_pi = df_tof30.Filter("abs(pdg) == 211")
    df_tof30_k = df_tof30.Filter("abs(pdg) == 321")
    df_tof30_p = df_tof30.Filter("abs(pdg) == 2212")

    N_MASS2_BINS, MIN_MASS2, MAX_MASS2 = 3000, -3, 3  # GeV^2/c^4
    h2d_tof30_pi = df_tof30_pi.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), "mom", "mass2")
    h2d_tof30_k = df_tof30_k.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), "mom", "mass2")
    h2d_tof30_p = df_tof30_p.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), "mom", "mass2")

    df_tof10 = df_tof.Define("beta", "trackLengthToEcal_IKF_zedLambda/(tofClosest10*299.792458)")\
                .Define("mass2", "mom*mom*( 1./(beta*beta) - 1.)")
    df_tof10_pi = df_tof10.Filter("abs(pdg) == 211")
    df_tof10_k = df_tof10.Filter("abs(pdg) == 321")
    df_tof10_p = df_tof10.Filter("abs(pdg) == 2212")

    h2d_tof10_pi = df_tof10_pi.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), "mom", "mass2")
    h2d_tof10_k = df_tof10_k.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), "mom", "mass2")
    h2d_tof10_p = df_tof10_p.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), "mom", "mass2")

    ################################################################ ENF OF RDATAFRAME

    _, eff_dedx_pik, _, sp_dedx_pik = calculate_pid_graphs(h2d_dedx_pi, h2d_dedx_k)
    _, eff_tof30_pik, _, sp_tof30_pik = calculate_pid_graphs(h2d_tof30_pi, h2d_tof30_k)
    sp_dedx_tof30_pik = combine_two_graphs(sp_dedx_pik, sp_tof30_pik)
    eff_dedx_tof30_pik = convert_graph_sp_to_eff( sp_dedx_tof30_pik )

    _, eff_tof10_pik, _, sp_tof10_pik = calculate_pid_graphs(h2d_tof10_pi, h2d_tof10_k)
    sp_dedx_tof10_pik = combine_two_graphs(sp_dedx_pik, sp_tof10_pik)
    eff_dedx_tof10_pik = convert_graph_sp_to_eff ( sp_dedx_tof10_pik )

    _, eff_dedx_kp, _, sp_dedx_kp = calculate_pid_graphs(h2d_dedx_k, h2d_dedx_p)
    _, eff_tof30_kp, _, sp_tof30_kp = calculate_pid_graphs(h2d_tof30_k, h2d_tof30_p)
    sp_dedx_tof30_kp = combine_two_graphs(sp_dedx_kp, sp_tof30_kp)
    eff_dedx_tof30_kp = convert_graph_sp_to_eff( sp_dedx_tof30_kp )

    _, eff_tof10_kp, _, sp_tof10_kp = calculate_pid_graphs(h2d_tof10_k, h2d_tof10_p)
    sp_dedx_tof10_kp = combine_two_graphs(sp_dedx_kp, sp_tof10_kp)
    eff_dedx_tof10_kp = convert_graph_sp_to_eff ( sp_dedx_tof10_kp )

    n_events_ss = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root").Filter("quarksToPythia == 33").Count()
    df_ss = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")\
            .Filter("quarksToPythia == 33")\
            .Filter("!isSimulated && !isOverlay")\
            .Filter("hasTrack && !isDecayedInTracker")\
            .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    n_events_gg = ROOT.RDataFrame("treeEvents", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h_dst.root").Filter("higgsDaughters == 909 || higgsDaughters == 2121").Count()
    df_gg = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h_dst.root")\
            .Filter("higgsDaughters == 909 || higgsDaughters == 2121")\
            .Filter("!isSimulated && !isOverlay")\
            .Filter("hasTrack && !isDecayedInTracker")\
            .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")

    if process == "Zss":
        n_events = n_events_ss
        df = df_ss
    elif process == "Hgg":
        n_events = n_events_gg
        df = df_gg
    else:
        raise("ERROR: Unknown process")
    df_pi = df.Filter("abs(pdg) == 211")
    df_k = df.Filter("abs(pdg) == 321")
    df_p = df.Filter("abs(pdg) == 2212")

    h_pi = df_pi.Histo1D((get_rand_string(), "Total;Momentum (GeV/c);Particles / event / 0.1 GeV/c", 200, 0, 20), "mom" )
    h_k = df_k.Histo1D((get_rand_string(), "Total;Momentum (GeV/c);Particles / event / 0.1 GeV/c", 200, 0, 20), "mom" )
    h_p = df_p.Histo1D((get_rand_string(), "Total;Momentum (GeV/c);Particles / event / 0.1 GeV/c", 200, 0, 20), "mom" )

    h_pi.Scale(1./n_events.GetValue())
    h_k.Scale(1./n_events.GetValue())
    h_p.Scale(1./n_events.GetValue())

    if particle == kaon:
        h_total_signal = h_k
    elif particle == proton:
        h_total_signal = h_p
    else:
        raise("Wrong particle")

    canvas = create_canvas()
    if use_dedx:
        # dE/dx kaon efficiency
        if particle == kaon:
            h_dedx_signal = apply_efficiency(h_k, eff_dedx_pik)
        elif particle == proton:
            h_dedx_signal = apply_efficiency(h_p, eff_dedx_kp)
        h_dedx_signal.SetLineColor( ROOT.TColor.GetColor("#ff0000") )
        h_dedx_signal.SetLineWidth(2)
        h_dedx_signal.SetLineStyle(2)

        # dE/dx + TOF 30 kaon efficiency
        if particle == kaon:
            h_dedx_tof30_signal = apply_efficiency(h_k, eff_dedx_tof30_pik)
        elif particle == proton:
            h_dedx_tof30_signal = apply_efficiency(h_p, eff_dedx_tof30_kp)
        h_dedx_tof30_signal.SetLineColor( ROOT.TColor.GetColor("#ff0000") )
        h_dedx_tof30_signal.SetLineWidth(2)
        h_dedx_tof30_signal.SetLineStyle(7)

        # dE/dx + TOF 10 kaon efficiency
        if particle == kaon:
            h_dedx_tof10_signal = apply_efficiency(h_k, eff_dedx_tof10_pik)
        elif particle == proton:
            h_dedx_tof10_signal = apply_efficiency(h_p, eff_dedx_tof10_kp)        
        h_dedx_tof10_signal.SetLineColor( ROOT.TColor.GetColor("#ff0000") )
        h_dedx_tof10_signal.SetLineWidth(2)

        # dE/dx pion mis-id
        if particle == kaon:
            h_dedx_bkg = apply_efficiency(h_pi, eff_dedx_pik, mode="misid")
        elif particle == proton:
            h_dedx_bkg = apply_efficiency(h_k, eff_dedx_kp, mode="misid")        
        h_dedx_bkg.SetLineColor( ROOT.TColor.GetColor("#0066ff") )
        h_dedx_bkg.SetLineWidth(2)
        h_dedx_bkg.SetLineStyle(2)


        # dE/dx + TOF 30 pion mis-id
        if particle == kaon:
            h_dedx_tof30_bkg = apply_efficiency(h_pi, eff_dedx_tof30_pik, mode="misid")
        elif particle == proton:
            h_dedx_tof30_bkg = apply_efficiency(h_k, eff_dedx_tof30_kp, mode="misid")
        h_dedx_tof30_bkg.SetLineColor( ROOT.TColor.GetColor("#0066ff") )
        h_dedx_tof30_bkg.SetLineWidth(2)
        h_dedx_tof30_bkg.SetLineStyle(7)


        # dE/dx + TOF 10 pion mis-id
        if particle == kaon:
            h_dedx_tof10_bkg = apply_efficiency(h_pi, eff_dedx_tof10_pik, mode="misid")
        elif particle == proton:
            h_dedx_tof10_bkg = apply_efficiency(h_k, eff_dedx_tof10_kp, mode="misid")
        h_dedx_tof10_bkg.SetLineColor( ROOT.TColor.GetColor("#0066ff") )
        h_dedx_tof10_bkg.SetLineWidth(2)

        h_total_signal.GetYaxis().SetMaxDigits(3)
        if particle == kaon:
            h_total_signal.GetXaxis().SetRangeUser(0, 8)
        elif particle == proton:
            h_total_signal.GetXaxis().SetRangeUser(0, 15)
        h_total_signal.DrawCopy("histo")
        h_dedx_signal.DrawCopy("histo same")
        h_dedx_tof30_signal.DrawCopy("histo same")
        h_dedx_tof10_signal.DrawCopy("histo same")
        h_dedx_bkg.DrawCopy("histo same")
        h_dedx_tof30_bkg.DrawCopy("histo same")
        h_dedx_tof10_bkg.DrawCopy("histo same")
        
        h_total_signal.SetMaximum(1.05*max([h.GetMaximum() for h in [h_total_signal, h_dedx_signal, h_dedx_tof30_signal, h_dedx_tof10_signal, h_dedx_bkg, h_dedx_tof30_bkg, h_dedx_tof10_bkg]]))

        leg = create_legend(0.2, 0.75, 0.76, 0.91)
        leg.SetNColumns(2)
        leg.SetMargin(0.15)

        if particle == kaon:
            signal_total_title = "K total (#varepsilon=100%)"
            signal_eff_title = "K identified"
            bkg_misid_title = "#pi misidentified"
        elif particle == proton:
            signal_total_title = "p total (#varepsilon=100%)"
            signal_eff_title = "p identified"
            bkg_misid_title = "K misidentified"
        else:
            raise("Wrong particle")

        h_leg1 = ROOT.TH1F(get_rand_string(), signal_total_title, 1, 0, 1)
        leg.AddEntry(h_leg1, h_leg1.GetTitle(), "l")

        h_leg2 = ROOT.TH1F(get_rand_string(), "dE/dx only", 1, 0, 1)
        h_leg2.SetLineStyle(2)
        leg.AddEntry(h_leg2, h_leg2.GetTitle(), "l")

        h_leg3 = ROOT.TH1F(get_rand_string(), signal_eff_title, 1, 0, 1)
        h_leg3.SetLineColor(ROOT.TColor.GetColor("#ff0000"))
        leg.AddEntry(h_leg3, h_leg3.GetTitle(), "l")

        h_leg4 = ROOT.TH1F(get_rand_string(), "dE/dx + TOF 30 ps", 1, 0, 1)
        h_leg4.SetLineStyle(7)
        leg.AddEntry(h_leg4, h_leg4.GetTitle(), "l")

        h_leg5 = ROOT.TH1F(get_rand_string(), bkg_misid_title, 1, 0, 1)
        h_leg5.SetLineColor(ROOT.TColor.GetColor("#0066ff"))
        leg.AddEntry(h_leg5, h_leg5.GetTitle(), "l")

        h_leg6 = ROOT.TH1F(get_rand_string(), "dE/dx + TOF 10 ps", 1, 0, 1)
        h_leg6.SetLineStyle(1)
        leg.AddEntry(h_leg6, h_leg6.GetTitle(), "l")
        leg.DrawCopy()

        latex.SetTextSize(0.04)
        if process == "Zss":
            latex.DrawLatexNDC(0.55, 0.55, "Z#rightarrows#bar{s} at E_{cm} = 250 GeV/c^{2}")
        elif process == "Hgg":
            latex.DrawLatexNDC(0.55, 0.55, "H#rightarrowgg at E_{cm} = 250 GeV/c^{2}")
    elif not use_dedx:
        # TOF 30 kaon efficiency
        if particle == kaon:
            h_tof30_signal = apply_efficiency(h_k, eff_tof30_pik)
        elif particle == proton:
            h_tof30_signal = apply_efficiency(h_p, eff_tof30_kp)
        h_tof30_signal.SetLineColor( ROOT.TColor.GetColor("#ff0000") )
        h_tof30_signal.SetLineWidth(2)
        h_tof30_signal.SetLineStyle(7)

        # TOF 10 kaon efficiency
        if particle == kaon:
            h_tof10_signal = apply_efficiency(h_k, eff_tof10_pik)
        elif particle == proton:
            h_tof10_signal = apply_efficiency(h_p, eff_tof10_kp)
        h_tof10_signal.SetLineColor( ROOT.TColor.GetColor("#ff0000") )
        h_tof10_signal.SetLineWidth(2)

        # TOF 30 pion mis-id
        if particle == kaon:
            h_tof30_bkg = apply_efficiency(h_pi, eff_tof30_pik, mode="misid")
        elif particle == proton:
            h_tof30_bkg = apply_efficiency(h_k, eff_tof30_kp, mode="misid")        
        h_tof30_bkg.SetLineColor( ROOT.TColor.GetColor("#0066ff") )
        h_tof30_bkg.SetLineWidth(2)
        h_tof30_bkg.SetLineStyle(7)

        # TOF 10 pion mis-id
        if particle == kaon:
            h_tof10_bkg = apply_efficiency(h_pi, eff_tof10_pik, mode="misid")
        elif particle == proton:
            h_tof10_bkg = apply_efficiency(h_k, eff_tof10_kp, mode="misid")        
        h_tof10_bkg.SetLineColor( ROOT.TColor.GetColor("#0066ff") )
        h_tof10_bkg.SetLineWidth(2)

        h_total_signal.GetYaxis().SetMaxDigits(3)
        if particle == kaon:
            h_total_signal.GetXaxis().SetRangeUser(0, 8)
        elif particle == proton:
            h_total_signal.GetXaxis().SetRangeUser(0, 15)
        h_total_signal.DrawCopy("histo")
        h_tof30_signal.DrawCopy("histo same")
        h_tof10_signal.DrawCopy("histo same")
        h_tof30_bkg.DrawCopy("histo same")
        h_tof10_bkg.DrawCopy("histo same")

        leg = create_legend(0.2, 0.75, 0.76, 0.91)
        leg.SetNColumns(2)
        leg.SetMargin(0.15)

        if particle == kaon:
            signal_total_title = "K total (#varepsilon=100%)"
            signal_eff_title = "K identified"
            bkg_misid_title = "#pi misidentified"
        elif particle == proton:
            signal_total_title = "p total (#varepsilon=100%)"
            signal_eff_title = "p identified"
            bkg_misid_title = "K misidentified"

        h_leg1 = ROOT.TH1F(get_rand_string(), signal_total_title, 1, 0, 1)
        leg.AddEntry(h_leg1, h_leg1.GetTitle(), "l")

        h_leg2 = ROOT.TH1F(get_rand_string(), "", 1, 0, 1)
        h_leg2.SetLineStyle(2)
        leg.AddEntry(0, "", "")

        h_leg3 = ROOT.TH1F(get_rand_string(), signal_eff_title, 1, 0, 1)
        h_leg3.SetLineColor(ROOT.TColor.GetColor("#ff0000"))
        leg.AddEntry(h_leg3, h_leg3.GetTitle(), "l")

        h_leg4 = ROOT.TH1F(get_rand_string(), "TOF 30 ps", 1, 0, 1)
        h_leg4.SetLineStyle(7)
        leg.AddEntry(h_leg4, h_leg4.GetTitle(), "l")

        h_leg5 = ROOT.TH1F(get_rand_string(), bkg_misid_title, 1, 0, 1)
        h_leg5.SetLineColor(ROOT.TColor.GetColor("#0066ff"))
        leg.AddEntry(h_leg5, h_leg5.GetTitle(), "l")

        h_leg6 = ROOT.TH1F(get_rand_string(), "TOF 10 ps", 1, 0, 1)
        h_leg6.SetLineStyle(1)
        leg.AddEntry(h_leg6, h_leg6.GetTitle(), "l")
        leg.DrawCopy()
        latex.SetTextSize(0.04)
        if process == "Zss":
            latex.DrawLatexNDC(0.55, 0.55, "Z#rightarrows#bar{s} at E_{cm} = 250 GeV/c^{2}")
        elif process == "Hgg":
            latex.DrawLatexNDC(0.55, 0.55, "H#rightarrowgg at E_{cm} = 250 GeV/c^{2}")
    canvas.Update()
    # input("wait")
    return canvas





c1 = tof_impact(use_dedx=True, process="Zss", particle=kaon)
c2 = tof_impact(use_dedx=True, process="Hgg", particle=kaon)
c3 = tof_impact(use_dedx=True, process="Zss", particle=proton)
c4 = tof_impact(use_dedx=True, process="Hgg", particle=proton)
c5 = tof_impact(use_dedx=False, process="Zss", particle=kaon)
c6 = tof_impact(use_dedx=False, process="Hgg", particle=kaon)
c7 = tof_impact(use_dedx=False, process="Zss", particle=proton)
c8 = tof_impact(use_dedx=False, process="Hgg", particle=proton)
input("wait")

def count_particle_species():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/2f_Z_hadronic_dst.root")\
            .Filter("quarksToPythia == 33")\
            .Filter("!isSimulated && !isOverlay")\
            .Filter("hasTrack && !isDecayedInTracker")\
            .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    n_total = df.Count()
    n_pi = df.Filter("abs(pdg) == 211").Count()
    n_k = df.Filter("abs(pdg) == 321").Count()
    n_p = df.Filter("abs(pdg) == 2212").Count()

    print("Pi fraction:", 100*n_pi.GetValue()/n_total.GetValue())
    print("K fraction:", 100*n_k.GetValue()/n_total.GetValue())
    print("P fraction:", 100*n_p.GetValue()/n_total.GetValue())

    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/higgs_Pn23n23h_dst.root")\
            .Filter("higgsDaughters == 909 || higgsDaughters == 2121")\
            .Filter("!isSimulated && !isOverlay")\
            .Filter("hasTrack && !isDecayedInTracker")\
            .Define("mom", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")
    n_total = df.Count()
    n_pi = df.Filter("abs(pdg) == 211").Count()
    n_k = df.Filter("abs(pdg) == 321").Count()
    n_p = df.Filter("abs(pdg) == 2212").Count()

    print("Pi fraction:", 100*n_pi.GetValue()/n_total.GetValue())
    print("K fraction:", 100*n_k.GetValue()/n_total.GetValue())
    print("P fraction:", 100*n_p.GetValue()/n_total.GetValue())

    input("wait")
# count_particle_species()