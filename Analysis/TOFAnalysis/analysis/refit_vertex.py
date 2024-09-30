from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.EnableImplicitMT()


def lcfiplus_good_track_selection():
    #TODO: THIS IS NOT A GOOD CODE, but it does what it should

    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
            .Filter("tofClosest0 > 6.")\
            .Define("pTrue", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")\
            .Filter("pTrue < 1.").Filter("abs(pdg) == 321")\
            .Define("recoPt", "sqrt(recoIpPx*recoIpPx + recoIpPy*recoIpPy)")\
            .Define("refittedRecoPt", "sqrt(refittedRecoIpPx*refittedRecoIpPx + refittedRecoIpPy*refittedRecoIpPy)")

    n = {}
    n["total"] = df.Count()

    # default pi mass fit
    n["d0_before"] = df.Filter("abs(d0IP) < 10").Count()
    n["z0_before"] = df.Filter("abs(d0IP) < 10 && abs(z0IP) < 20").Count()
    n["d0_err_before"] = df.Filter("abs(d0IP) < 10 && abs(z0IP) < 20 && d0ErrIP < 0.1").Count()
    n["z0_err_before"] = df.Filter("abs(d0IP) < 10 && abs(z0IP) < 20 && d0ErrIP < 0.1 && z0ErrIP < 0.1").Count()
    n["pt_before"] = df.Filter("abs(d0IP) < 10 && abs(z0IP) < 20 && d0ErrIP < 0.1 && z0ErrIP < 0.1 && recoPt > 0.1").Count()
    # n["from_ip_before"] = df.Filter("abs(d0IP) < 10 && abs(z0IP) < 20 && d0ErrIP < 0.1 && z0ErrIP < 0.1 && recoPt > 0.1 && !isInRecoPrimaryVertex").Count()

    n["d0_after"] = df.Filter("abs(refittedD0IP) < 10").Count()
    n["z0_after"] = df.Filter("abs(refittedD0IP) < 10 && abs(refittedZ0IP) < 20").Count()
    n["d0_err_after"] = df.Filter("abs(refittedD0IP) < 10 && abs(refittedZ0IP) < 20 && refittedD0ErrIP < 0.1").Count()
    n["z0_err_after"] = df.Filter("abs(refittedD0IP) < 10 && abs(refittedZ0IP) < 20 && refittedD0ErrIP < 0.1 && refittedZ0ErrIP < 0.1").Count()
    n["pt_after"] = df.Filter("abs(refittedD0IP) < 10 && abs(refittedZ0IP) < 20 && refittedD0ErrIP < 0.1 && refittedZ0ErrIP < 0.1 && refittedRecoPt > 0.1").Count()
    # n["from_ip_after"] = df.Filter("abs(refittedD0IP) < 10 && abs(refittedZ0IP) < 20 && abs(refittedD0ErrIP) < 0.1 && abs(refittedZ0ErrIP) < 0.1 && refittedRecoPt > 0.1 && !isInRecoPrimaryRefitVertex").Count()

    h_before = ROOT.TH1F("h_before", "default", 6, 0, 6)
    h_before.GetXaxis().SetBinLabel(1, "no cut")
    h_before.GetXaxis().SetBinLabel(2, "d_{0} cut")
    h_before.GetXaxis().SetBinLabel(3, "z_{0} cut")
    h_before.GetXaxis().SetBinLabel(4, "#sigma_{d0} cut")
    h_before.GetXaxis().SetBinLabel(5, "#sigma_{z0} cut")
    h_before.GetXaxis().SetBinLabel(6, "p_{T} cut")
    # h_before.GetXaxis().SetBinLabel(7, "prim.vtx")
    h_before.LabelsDeflate("X")
    h_before.LabelsOption("v")
    h_before.SetBinContent(1, n["total"].GetValue() )
    h_before.SetBinContent(2, n["d0_before"].GetValue() )
    h_before.SetBinContent(3, n["z0_before"].GetValue() )
    h_before.SetBinContent(4, n["d0_err_before"].GetValue() )
    h_before.SetBinContent(5, n["z0_err_before"].GetValue() )
    h_before.SetBinContent(6, n["pt_before"].GetValue() )
    # h_before.SetBinContent(7, n["from_ip_before"].GetValue() )
    h_after = ROOT.TH1F("h_after", "refitted", 6, 0, 6)
    h_after.SetBinContent(1, n["total"].GetValue() )
    h_after.SetBinContent(2, n["d0_after"].GetValue() )
    h_after.SetBinContent(3, n["z0_after"].GetValue() )
    h_after.SetBinContent(4, n["d0_err_after"].GetValue() )
    h_after.SetBinContent(5, n["z0_err_after"].GetValue() )
    h_after.SetBinContent(6, n["pt_after"].GetValue() )
    # h_after.SetBinContent(7, n["from_ip_after"].GetValue() )
    h_before.Scale(100./n["total"].GetValue())
    h_after.Scale(100./n["total"].GetValue())

    c = create_canvas()
    h_before.Draw("PLhist")
    h_before.GetYaxis().SetTitle("Fraction of kaons passed (%)")
    h_before.GetYaxis().SetTitleOffset(1.2)
    h_before.GetXaxis().SetLabelOffset(0.005)
    h_before.GetYaxis().SetRangeUser(61, 104)
    h_before.GetYaxis().SetNdivisions(512)
    h_after.Draw("PLhist same")

    h_after.SetLineColor(ROOT.kRed)
    h_after.SetMarkerColor(ROOT.kRed)
    h_after.SetMarkerStyle(29)
    h_before.SetMarkerSize(1.7)
    h_after.SetMarkerSize(1.5)
    leg = c.BuildLegend()
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)

    latex.DrawLatexNDC(0.23, 0.67, "p < 1 GeV/c")
    c.Update()
    input("wait")


# lcfiplus_good_track_selection()


def vertex_position_resolution():
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
        .Define("pTrue", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")\
        .Filter("pTrue < 2.")\
        .Filter("abs(pdg) == 321")\
        .Filter("isInTrueSecondaryVertex && isInRecoSecondaryVertex && isInRecoSecondaryRefitVertex")\
        .Filter("tofClosest0 > 6.")\
        .Filter("nTracksAtRecoVertex == nTracksAtRecoRefitVertex && nTracksAtRecoVertex == nTracksAtTrueVertex")

    h_default = df.Define("diff", "(recoVertexPosZ - trueVertexPosZ)/recoVertexPosErrZ")\
                        .Histo1D((get_rand_string(), "Default;(_{}z_{reco} - z_{true}) / #sigma_{z};N entries", 200, -5, 5), "diff")
    h_refitted = df.Define("diff", "(recoRefitVertexPosZ - trueVertexPosZ)/recoRefitVertexPosErrZ")\
                        .Histo1D((get_rand_string(), "Refitted;(_{}z_{reco} - z_{true}) / #sigma_{z};N entries", 200, -5, 5), "diff")

    c = create_canvas()
    h_default.SetTitle(f"Default  Std Dev: {h_default.GetStdDev():.3f}")
    h_refitted.SetTitle(f"Refitted  Std Dev: {h_refitted.GetStdDev():.3f}")
    h_refitted.SetLineColor(ROOT.kRed)
    h_default.GetXaxis().SetTitleOffset(1.)
    h_default.GetYaxis().SetMaxDigits(3)
    h_default.SetMinimum( 0. )
    h_default.SetMaximum( 1.3*max(h_default.GetMaximum(), h_refitted.GetMaximum()) )
    h_default.Draw()
    h_refitted.Draw("same")
    leg = create_legend(0.19, 0.77, 0.999, 0.91)
    leg.AddEntry(h_default.GetPtr(), h_default.GetTitle(), "l")
    leg.AddEntry(h_refitted.GetPtr(), h_refitted.GetTitle(), "l")
    leg.Draw()
    latex.DrawLatexNDC(0.23, 0.67, "p < 2 GeV/c")
    c.Update()
    input("wait")
vertex_position_resolution()