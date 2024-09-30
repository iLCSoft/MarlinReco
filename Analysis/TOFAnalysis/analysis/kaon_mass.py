from utils import *
ROOT.EnableImplicitMT()

colors = ["#0099ff", "#128ce9", "#1d77c5", "#225d98", "#214772", "#1f3c5f", "#1c314d"]
colors = [ ROOT.TColor.GetColor(c) for c in colors]


def main():
    df_init = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
                  .Filter("abs(pdg) == 321")\
                  .Filter("tofClosest0 > 6.")

    resolutions = [1, 5, 10, 20, 30]
    histos = {}
    for res in resolutions:
        df_res = df_init.Define("beta", f"trackLengthToEcal_IKF_zedLambda/(tofClosest{res}*299.792458)")\
                    .Define("mass2", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda*( 1./(beta*beta) - 1.)")\
                    .Define("mass", "sqrt(mass2)*1000")

        h = df_res.Histo1D((get_rand_string(), "", 200, 452, 538),"mass")
        histos[res] = h


    canvas = create_canvas()
    legend = create_legend(0.4, 0.63, 0.98, 0.94)
    for i, (res, h) in enumerate(histos.items()):
        if i == 0:
            h.Draw("L")
            h.SetTitle(";Mass (MeV/c^{2}); N entries")
            h.GetYaxis().SetMaxDigits(3)
            h.GetXaxis().SetNdivisions(508)
        else:
            h.Draw("Lsame")
        h.SetLineColor(colors[i])
        h.SetMarkerColor(colors[i])
        h.SetLineWidth(4)
        legend.AddEntry(h.GetPtr(), f"{res} ps","l")

    legend.Draw()

    latex = ROOT.TLatex()
    latex.SetTextFont(52)
    # latex.DrawLatex(12, 6.06, "ILD preliminary")

    line = ROOT.TLine(kaon.mass*1e3, 0., kaon.mass*1e3, histos[1].GetMaximum())
    line.SetLineColor(15)
    line.SetLineWidth(2)
    line.SetLineStyle(9)
    line.Draw()

    canvas.Update()
    input("wait")


def simulate_bias():
    particle = pion
    momentum = 1 # GeV/c
    track_length = 2000. # mm
    tof_true = track_length/SPEED_OF_LIGHT * np.sqrt( 1 + particle.mass2/(momentum*momentum) ) * 1000 # ps

    canvas = create_canvas()
    for i, tof_res in enumerate([0, 1, 5, 100, 300]):
        tofs = np.random.normal(tof_true, tof_res, 10000000)/1000. # in ns
        masses2 = momentum * momentum * ( (SPEED_OF_LIGHT*tofs/track_length)**2 -1 )
        masses = np.sqrt(masses2)*1000 # in MeV/c^{2}
        h = ROOT.TH1F(f"h{tof_res}", f"Resolution {tof_res} ps", 1000, particle.mass*1000-5, particle.mass*1000+5)
        h.FillN(len(masses), masses, np.ones_like(masses))
        h.SetLineColor(colors[i])
        h.Scale(1./h.GetEntries())
        h.DrawCopy("histo" if i == 0 else "histo same")

    canvas.BuildLegend(0.5, 0.5, 0.8, 0.8, "", "l")
    canvas.Update()
    input("wait")


simulate_bias()
# main()
