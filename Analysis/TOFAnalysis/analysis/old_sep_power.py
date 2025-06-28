import ROOT
import numpy as np
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(1)
ROOT.EnableImplicitMT()

colors = ['#1b9e77', '#d95f02', '#7570b3']
colors = [ ROOT.TColor.GetColor(c) for c in colors]

def get_separation_power(mu1, mu2, sigma1, sigma2, mu1_error=0, mu2_error=0):
    dm = abs( mu1 - mu2 )
    ds = sigma1*sigma1 + sigma2*sigma2
    if ds == 0:
        return -1
    sds = np.sqrt(0.5*ds)
    spe_1 = 2*mu1_error*mu1_error/ds
    spe_2 = 2*mu2_error*mu2_error/ds
    spe_3 = 2*( (dm*sigma1*mu1_error)**2 )/(ds*ds*ds)
    spe_4 = 2*( (dm*sigma2*mu2_error)**2 )/(ds*ds*ds)

    sep_power = 0
    sep_power_error = 0
    sep_power = dm/sds
    sep_power_error = np.sqrt( spe_1 + spe_2 + spe_3 + spe_4 )
    return sep_power, sep_power_error


df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
df = df.Filter("tofClosest0 > 6.")

#linear binning
n_mom_bins = 50
min_x = 0
max_x = 20
# n_mom_bins = 30
# min_x = 0
# max_x = 15

#log binning
# n_mom_bins = 50
# min_x = -1 # 0.1 GeV
# max_x = 1.3 # ~ 20 GeV
# x_bins = np.array( [ 10**(min_x + (max_x-min_x)*i/n_mom_bins ) for i in range(n_mom_bins+1) ] )


def get_sep_power_graph(df, tof_column="tofClosest0"):
    df_total = df.Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")\
                .Filter("beta >= 0 && beta <= 1")\
                .Define("mass", "harmonicMomToEcal_IKF_zedLambda*sqrt( 1./(beta*beta) - 1.)*1000")
    h_2d_total = df_total.Histo2D(("h_total", "All; momentum [GeV]; Mass [MeV]", n_mom_bins, min_x, max_x, 1000, 0, 6000), "harmonicMomToEcal_IKF_zedLambda","mass")

    df_pion = df.Filter("abs(pdg) == 211")\
                .Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")\
                .Filter("beta >= 0 && beta <= 1")\
                .Define("mass", "harmonicMomToEcal_IKF_zedLambda*sqrt( 1./(beta*beta) - 1.)*1000")
    h_2d_pion = df_pion.Histo2D(("h_pion", "#pi; momentum [GeV]; Mass [MeV]", n_mom_bins, min_x, max_x, 1000, 0, 6000), "harmonicMomToEcal_IKF_zedLambda","mass")

    df_kaon = df.Filter("abs(pdg) == 321")\
                .Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")\
                .Filter("beta >= 0 && beta <= 1")\
                .Define("mass", "harmonicMomToEcal_IKF_zedLambda*sqrt( 1./(beta*beta) - 1.)*1000")
    h_2d_kaon = df_kaon.Histo2D(("h_kaon", "K; momentum [GeV]; Mass [MeV]", n_mom_bins, min_x, max_x, 1000, 0, 6000), "harmonicMomToEcal_IKF_zedLambda","mass")

    df_proton = df.Filter("abs(pdg) == 2212")\
                .Define("beta", f"trackLengthToEcal_IKF_zedLambda/({tof_column}*299.792458)")\
                .Filter("beta >= 0 && beta <= 1")\
                .Define("mass", "harmonicMomToEcal_IKF_zedLambda*sqrt( 1./(beta*beta) - 1.)*1000")
    h_2d_proton = df_proton.Histo2D(("h_proton", "proton; momentum [GeV]; Mass [MeV]", n_mom_bins, min_x, max_x, 1000, 0, 6000), "harmonicMomToEcal_IKF_zedLambda","mass")

    # DRAW FANCY 2D plot
    # ROOT.gStyle.SetPadRightMargin(0.12)
    # canvas = ROOT.TCanvas("c_2d_m_vs_p_total",
    #                         "",
    #                         int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
    #                         int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
    # h_2d_total.Draw("colz")
    # h_2d_total.SetMinimum(1)
    # h_2d_total.SetMaximum(100000)
    # canvas.SetLogz()
    # canvas.SetGridx(0)
    # canvas.SetGridy(0)
    # canvas.Update()
    # palette = h_2d_total.GetListOfFunctions().FindObject("palette")
    # palette.SetX1NDC(0.89)
    # palette.SetX2NDC(0.91)
    # canvas.Modified()
    # canvas.Update()

    ### FIT HISTOGRAMS IN SLICES ###
    gaus_tof = ROOT.TF1("gaus_tof", "gaus")
    h_slices_pions = ROOT.TObjArray()
    h_2d_pion.FitSlicesY(gaus_tof, 0, -1, 0, "QN", h_slices_pions)
    h_slices_kaons = ROOT.TObjArray()
    h_2d_kaon.FitSlicesY(gaus_tof, 0, -1, 0, "QN", h_slices_kaons)
    h_slices_protons = ROOT.TObjArray()
    h_2d_proton.FitSlicesY(gaus_tof, 0, -1, 0, "QN", h_slices_protons)
    ### DONE ###

    ### Get mean and sigmas ###
    h_pion_mean = h_slices_pions.At(1)
    h_pion_sigma = h_slices_pions.At(2)
    h_kaon_mean = h_slices_kaons.At(1)
    h_kaon_sigma = h_slices_kaons.At(2)
    h_proton_mean = h_slices_protons.At(1)
    h_proton_sigma = h_slices_protons.At(2)
    ### DONE ###

    ### Calculate Sep. Power and plot Graphs ###
    gr_sp_pik = ROOT.TGraph()
    gr_sp_kp = ROOT.TGraph()
    gr_sp_pik.SetTitle("; momentum [GeV]; #pi/K separation power")
    gr_sp_kp.SetTitle("; momentum [GeV]; K/p separation power")

    for i in range(1, n_mom_bins+1):
        x = h_pion_mean.GetBinCenter(i)
        mu_pion = h_pion_mean.GetBinContent(i)
        mu_kaon = h_kaon_mean.GetBinContent(i)
        mu_proton = h_proton_mean.GetBinContent(i)

        sigma_pion = h_pion_sigma.GetBinContent(i)
        sigma_kaon = h_kaon_sigma.GetBinContent(i)
        sigma_proton = h_proton_sigma.GetBinContent(i)

        sep_power_pik, _ = get_separation_power(mu_pion, mu_kaon, sigma_pion, sigma_kaon)
        sep_power_kp, _ = get_separation_power(mu_kaon, mu_proton, sigma_kaon, sigma_proton)

        gr_sp_pik.SetPoint(i, x, sep_power_pik)
        gr_sp_kp.SetPoint(i, x, sep_power_kp)

    # gr_sp_pik.Draw("APL")
    # gr_sp_kp.Draw("PLsame")
    # gr_sp_kp.SetLineColor(4)
    # canvas2.Update()
    # input("wait")
    return gr_sp_kp

canvas = ROOT.TCanvas("c",
                        "",
                        int(600/(1. - ROOT.gStyle.GetPadLeftMargin() - ROOT.gStyle.GetPadRightMargin())),
                        int(600/(1. - ROOT.gStyle.GetPadTopMargin() - ROOT.gStyle.GetPadBottomMargin())) )
canvas.SetGridx(True)
canvas.SetGridy(True)

gr_sp={}
colors = ["#03045e","#023e8a","#0077b6","#0096c7","#00b4d8","#48cae4"]
colors = [ ROOT.TColor.GetColor(c) for c in colors[::-1]]

legend = ROOT.TLegend()
for i, res in enumerate( [0, 10, 30, 50, 70, 90] ):
    gr_sp[res] = get_sep_power_graph(df, tof_column=f"tofClosest{int(res)}")
    gr_sp[res].Draw("AL" if i == 0 else "Lsame")
    gr_sp[res].SetLineColor(colors[i])
    gr_sp[res].SetLineWidth(4)
    legend.AddEntry(gr_sp[res], f"{res} ps","l")

gr_sp[0].GetXaxis().SetRangeUser(0, 19)
gr_sp[0].GetYaxis().SetRangeUser(0, 12)
legend.Draw()
canvas.Update()
input("wait")













# get histograms in momentum slices


# n_y_bins = 500
# min_y = 0
# max_y = 1e-6
# y_bins = [ min_y + (max_y - min_y)*i/n_y_bins for i in range(n_y_bins+1) ]

# h_bb_pions = df_pions.Filter("dEdx > 0").Histo2D(("h_bb_pions", "title; momentum (GeV); dEdx (GeV/mm)", n_x_bins, np.array(x_bins), n_y_bins, np.array(y_bins)), "momentumIP", "dEdx")
# h_bb_kaons = df_kaons.Filter("dEdx > 0").Histo2D(("h_bb_kaons", "title; momentum (GeV); dEdx (GeV/mm)", n_x_bins, np.array(x_bins), n_y_bins, np.array(y_bins)), "momentumIP", "dEdx")
# h_bb_protons = df_protons.Filter("dEdx > 0").Histo2D(("h_bb_protons", "title; momentum (GeV); dEdx (GeV/mm)", n_x_bins, np.array(x_bins), n_y_bins, np.array(y_bins)), "momentumIP", "dEdx")

# h_tof_pions = df_pions.Filter("tof0 > 6. && nClusters == 1").Define("beta", "trackLength/(tof30*299.792458)").Filter("beta >= 0 && beta <= 1").Define("mass", "momentumHM*sqrt( 1./(beta*beta) - 1.)*1000").Histo2D(("h_tof_pions", "title; momentum (GeV); mass (MeV)", n_x_bins, np.array(x_bins), 500, -100, 1300.), "momentumIP", "mass")
# h_tof_kaons = df_kaons.Filter("tof0 > 6. && nClusters == 1").Define("beta", "trackLength/(tof30*299.792458)").Filter("beta >= 0 && beta <= 1").Define("mass", "momentumHM*sqrt( 1./(beta*beta) - 1.)*1000").Histo2D(("h_tof_kaons", "title; momentum (GeV); mass (MeV)", n_x_bins, np.array(x_bins), 500, -100, 1300.), "momentumIP", "mass")
# h_tof_protons = df_protons.Filter("tof0 > 6. && nClusters == 1").Define("beta", "trackLength/(tof30*299.792458)").Filter("beta >= 0 && beta <= 1").Define("mass", "momentumHM*sqrt( 1./(beta*beta) - 1.)*1000").Histo2D(("h_tof_protons", "title; momentum (GeV); mass (MeV)", n_x_bins, np.array(x_bins), 500, -100, 1300.), "momentumIP", "mass")

# h_tof_2d = df.Filter("tof30 > 6. && nClusters == 1 && (abs(pdg) == 211 ||abs(pdg) == 321 || abs(pdg) == 2212)")\
#              .Define("beta", "trackLength/(tof30*299.792458)").Filter("beta >= 0 && beta <= 1")\
#              .Define("mass", "momentumHM*sqrt( 1./(beta*beta) - 1.)*1000")\
#              .Histo2D(("h_tof", "Mass reconstructed with time-of-flight assuming 30 ps per particle; momentum (GeV); mass (MeV)", n_x_bins, np.array(x_bins), 500, -100, 1300.), "momentumIP", "mass")


# get bin centers
# h_tmp = ROOT.TH1F("h_tmp", "title", n_x_bins, np.array(x_bins))
# x_bin_centers = [ h_tmp.GetBinCenter(i) for i in range(1, n_x_bins+1) ]
# h_tof_2d.Draw("colz")
# style_2dhistogram(h_tof_2d)
# input("wait")

#experimenting with linear scale
# h_bb_pions = df.Histo2D(("h_bb_pions", "title; momentum (GeV); dEdx (GeV/mm)", 1000, 0, 20, 300, 0, 1e-6), "momentumIP", "dEdx")

# h_mom = df.Histo1D(("h_mom", "title; momentum (GeV); N particles", 1000, 0, 20), "momentumIP")
# # MOMENTUM PLOT
# h_mom_pions = df_pions.Histo1D(("h_mom_pions", "Pions; momentum (GeV); N particles",300, 0, 10), "momentumIP")
# h_mom_kaons = df_kaons.Histo1D(("h_mom_kaons", "Kaons; momentum (GeV); N particles",300, 0, 10), "momentumIP")
# h_mom_protons = df_protons.Histo1D(("h_mom_protons", "Protons; momentum (GeV); N particles",300, 0, 10), "momentumIP")
# h_mom_pions.Scale( 1. / h_mom_pions.GetEntries() )
# h_mom_kaons.Scale( 1. / h_mom_kaons.GetEntries() )
# h_mom_protons.Scale( 1. / h_mom_protons.GetEntries() )
# h_mom_pions.Draw("histo ")
# h_mom_pions.SetLineColor(colors[0])
# h_mom_pions.SetLineWidth(5)
# h_mom_protons.Draw("histo same")
# h_mom_protons.SetLineColor(colors[2])
# h_mom_protons.SetLineWidth(5)
# h_mom_kaons.Draw("histo same")
# h_mom_kaons.SetLineColor(colors[1])
# h_mom_kaons.SetLineWidth(5)

# print("Bin number for 3 GeV: ", h_mom_pions.FindBin(3.))
# print("Pion integral below 3 GeV: ", h_mom_pions.Integral(0, h_mom_pions.FindBin(3.)))
# print("Kaon integral below 3 GeV: ", h_mom_kaons.Integral(0, h_mom_kaons.FindBin(3.)))
# print("Proton integral below 3 GeV: ", h_mom_protons.Integral(0, h_mom_protons.FindBin(3.)))
# input("wait")
# # END OF MOMENTUM PLOT

# gaus = ROOT.TF1("gaus", "gaus")
# gaus.SetParameters(100, max_y/3.,max_y/10.)
# gaus.SetParLimits(0, 0, 1e5)
# gaus.SetParLimits(1, 0, max_y)
# gaus.SetParLimits(2, 0, max_y/2.)

# fits_pions = ROOT.TObjArray()
# h_bb_pions.FitSlicesY(gaus, 0, -1, 0, "QN", fits_pions)
# fits_kaons = ROOT.TObjArray()
# h_bb_kaons.FitSlicesY(gaus, 0, -1, 0, "QN", fits_kaons)
# fits_protons = ROOT.TObjArray()
# h_bb_protons.FitSlicesY(gaus, 0, -1, 0, "QN", fits_protons)

# gaus_tof = ROOT.TF1("gaus_tof", "gaus")
# # gaus_tof.SetParameters(100, max_y/3.,max_y/10.)
# # gaus_tof.SetParLimits(0, 0, 1e5)
# # gaus_tof.SetParLimits(1, 0, max_y)
# # gaus_tof.SetParLimits(2, 0, max_y/2.)
# fits_tof_pions = ROOT.TObjArray()
# h_tof_pions.FitSlicesY(gaus_tof, 0, -1, 0, "QN", fits_tof_pions)
# fits_tof_kaons = ROOT.TObjArray()
# h_tof_kaons.FitSlicesY(gaus_tof, 0, -1, 0, "QN", fits_tof_kaons)
# fits_tof_protons = ROOT.TObjArray()
# h_tof_protons.FitSlicesY(gaus_tof, 0, -1, 0, "QN", fits_tof_protons)

 
# #pions vs kaons dE/dx
# # h_sp_pik = ROOT.TH1F("h_sp_pik", "title; momentum (GeV); Separation power", n_x_bins, np.array(x_bins))
# # h_sp_kp = ROOT.TH1F("h_sp_kp", "title; momentum (GeV); Separation power", n_x_bins, np.array(x_bins))
# # Trying graphs
# gr_sp_pik = ROOT.TGraph()
# gr_sp_kp = ROOT.TGraph()
# gr_sp_pik.SetTitle("title; momentum (GeV); Separation power")
# gr_sp_kp.SetTitle("title; momentum (GeV); Separation power")
# gr_sp_pik.SetPoint(0, x_bin_centers[0], 0)
# gr_sp_kp.SetPoint(0, x_bin_centers[0], 0)

# fit_pion_mean = fits_pions.At(1)
# fit_pion_sigma = fits_pions.At(2)
# fit_kaon_mean = fits_kaons.At(1)
# fit_kaon_sigma = fits_kaons.At(2)
# fit_proton_mean = fits_protons.At(1)
# fit_proton_sigma = fits_protons.At(2)

# for i in range(1, n_x_bins+1):
#     mu_pion = fit_pion_mean.GetBinContent(i)
#     mu_kaon = fit_kaon_mean.GetBinContent(i)
#     mu_proton = fit_proton_mean.GetBinContent(i)

#     sigma_pion = fit_pion_sigma.GetBinContent(i)
#     sigma_kaon = fit_kaon_sigma.GetBinContent(i)
#     sigma_proton = fit_proton_sigma.GetBinContent(i)

#     sep_power, sep_power_error = get_separation_power(mu_pion, mu_kaon, sigma_pion, sigma_kaon)
#     # h_sp_pik.SetBinContent(i, sep_power)
#     gr_sp_pik.SetPoint(i, x_bin_centers[i-1], sep_power)

#     sep_power, sep_power_error = get_separation_power(mu_kaon, mu_proton, sigma_kaon, sigma_proton)
#     # h_sp_kp.SetBinContent(i, sep_power)
#     gr_sp_kp.SetPoint(i, x_bin_centers[i-1], sep_power)

# gr_sp_pik.SetPoint(n_x_bins+1, x_bin_centers[n_x_bins-1], 0)
# gr_sp_kp.SetPoint(n_x_bins+1, x_bin_centers[n_x_bins-1], 0)



# ################ DO TOF separation power in the same bins
# # h_sp_pik_tof = ROOT.TH1F("h_sp_pik_tof", "title; momentum (GeV); Separation power", n_x_bins, np.array(x_bins))
# # h_sp_kp_tof = ROOT.TH1F("h_sp_kp_tof", "title; momentum (GeV); Separation power", n_x_bins, np.array(x_bins))
# gr_sp_pik_tof = ROOT.TGraph()
# gr_sp_kp_tof = ROOT.TGraph()
# gr_sp_pik_tof.SetTitle("#pi / K; momentum (GeV); Separation power")
# gr_sp_kp_tof.SetTitle("K / p; momentum (GeV); Separation power")
# gr_sp_pik_tof.SetPoint(0, x_bin_centers[0], 0)
# gr_sp_kp_tof.SetPoint(0, x_bin_centers[0], 0)

# fit_pion_mean = fits_tof_pions.At(1)
# fit_pion_sigma = fits_tof_pions.At(2)
# fit_kaon_mean = fits_tof_kaons.At(1)
# fit_kaon_sigma = fits_tof_kaons.At(2)
# fit_proton_mean = fits_tof_protons.At(1)
# fit_proton_sigma = fits_tof_protons.At(2)

# for i in range(1, n_x_bins+1):
#     mu_pion = fit_pion_mean.GetBinContent(i)
#     mu_kaon = fit_kaon_mean.GetBinContent(i)
#     mu_proton = fit_proton_mean.GetBinContent(i)

#     sigma_pion = fit_pion_sigma.GetBinContent(i)
#     sigma_kaon = fit_kaon_sigma.GetBinContent(i)
#     sigma_proton = fit_proton_sigma.GetBinContent(i)

#     sep_power, sep_power_error = get_separation_power(mu_pion, mu_kaon, sigma_pion, sigma_kaon)
#     # h_sp_pik_tof.SetBinContent(i, sep_power)
#     gr_sp_pik_tof.SetPoint(i, x_bin_centers[i-1], sep_power)

#     sep_power, sep_power_error = get_separation_power(mu_kaon, mu_proton, sigma_kaon, sigma_proton)
#     # h_sp_kp_tof.SetBinContent(i, sep_power)
#     gr_sp_kp_tof.SetPoint(i, x_bin_centers[i-1], sep_power)

# gr_sp_pik_tof.SetPoint(n_x_bins+1, x_bin_centers[n_x_bins-1], 0)
# gr_sp_kp_tof.SetPoint(n_x_bins+1, x_bin_centers[n_x_bins-1], 0)


# # For combined sep power
# # h_sp_pik_combined = ROOT.TH1F("h_sp_pik_comb", "#pi / K combined; momentum (GeV); Separation power", n_x_bins, np.array(x_bins))
# # h_sp_kp_combined = ROOT.TH1F("h_sp_kp_comb", "K /p combined; momentum (GeV); Separation power", n_x_bins, np.array(x_bins))
# gr_sp_pik_combined = ROOT.TGraph()
# gr_sp_kp_combined = ROOT.TGraph()
# gr_sp_pik_combined.SetTitle("#pi / K combined; momentum (GeV); Separation power")
# gr_sp_kp_combined.SetTitle("K / p combined; momentum (GeV); Separation power")

# for i in range(0, n_x_bins+2):
#     # sp_dedx = h_sp_pik.GetBinContent(i)
#     # sp_tof = h_sp_pik_tof.GetBinContent(i)
#     # h_sp_pik_combined.SetBinContent(i, np.sqrt(sp_dedx*sp_dedx + sp_tof*sp_tof) )
#     x = gr_sp_pik.GetPointX(i)

#     sp_dedx = gr_sp_pik.GetPointY(i)
#     sp_tof = gr_sp_pik_tof.GetPointY(i)
#     gr_sp_pik_combined.SetPoint(i, x, np.sqrt(sp_dedx*sp_dedx + sp_tof*sp_tof))

#     # sp_dedx = h_sp_kp.GetBinContent(i)
#     # sp_tof = h_sp_kp_tof.GetBinContent(i)
#     # h_sp_kp_combined.SetBinContent(i, np.sqrt(sp_dedx*sp_dedx + sp_tof*sp_tof) )

#     sp_dedx = gr_sp_kp.GetPointY(i)
#     sp_tof = gr_sp_kp_tof.GetPointY(i)
#     gr_sp_kp_combined.SetPoint(i, x, np.sqrt(sp_dedx*sp_dedx + sp_tof*sp_tof))

# canvas = ROOT.TCanvas()
# # TOF ONLY SEPARATION POWERS
# gr_sp_pik_tof.Draw("AL")
# gr_sp_pik_tof.SetLineColor(colors[1])
# gr_sp_pik_tof.SetLineStyle(2)
# gr_sp_pik_tof.SetLineWidth(7)

# gr_sp_kp_tof.Draw("LP same")
# gr_sp_kp_tof.SetLineColor(colors[1])
# gr_sp_kp_tof.SetLineStyle(10)
# gr_sp_kp_tof.SetLineWidth(7)

# # TOF AND dEdX combined

# # gr_sp_kp.Draw("ALF")
# # gr_sp_kp.SetLineColor(colors[0])
# # gr_sp_kp.SetFillColor(colors[0])
# # gr_sp_kp.SetLineWidth(3)
# # gr_sp_kp.SetFillStyle(3002)

# # gr_sp_kp_tof.Draw("LF same")
# # gr_sp_kp_tof.SetLineColor(colors[1])
# # gr_sp_kp_tof.SetFillColor(colors[1])
# # gr_sp_kp_tof.SetLineWidth(3)
# # gr_sp_kp_tof.SetFillStyle(3011)

# # gr_sp_kp_combined.Draw("L")
# # gr_sp_kp_combined.SetLineColor(colors[2])
# # gr_sp_kp_combined.SetLineWidth(5)


# # h_sp_kp.Draw("histo")
# # h_sp_kp.SetLineWidth(3)
# # h_sp_kp.SetLineColor(2)
# # h_sp_kp_tof.SetLineWidth(5)
# # h_sp_kp_tof.SetLineColor(4)
# # h_sp_kp_combined.Draw("histo same")
# # h_sp_kp_combined.SetLineWidth(3)
# # h_sp_kp_combined.SetLineColor(6)

# # h_sp_pik_tof.SetMarkerStyle(20)
# # h_sp_pik_tof.SetMarkerColor(4)

# # h_sp_kp_tof.Draw("histo same")
# # h_sp_kp_tof.SetLineWidth(3)
# # h_sp_kp_tof.SetLineColor(8)
# # h_sp_kp_tof.SetMarkerStyle(20)
# # h_sp_kp_tof.SetMarkerColor(8)

# canvas.BuildLegend()
# canvas.SetLogx()
# input("wait")













# canvas_bb = ROOT.TCanvas()
# canvas_bb.SetLogx()
# h_tof_pions.Draw("colz")
# canvas_bb.Update()
# # style_2dhistogram(h_bb_pions)
# style_2dhistogram(h_tof_pions)
# canvas_bb.Modified()
# canvas_bb.Update()

# canvas_fit = ROOT.TCanvas()
# canvas_fit.Divide(2, 2)
# for i in range(4):
#     canvas_fit.cd(i+1)
#     ROOT.gPad.SetLogx()
#     fits_tof_pions.At(i).SetStats(0)
#     fits_tof_pions.At(i).SetMinimum(0)
#     # if i==0 :
#     #     fits_tof_pions.At(i).SetMaximum(10000)
#     # elif i==1:
#     #     fits_tof_pions.At(i).SetMaximum(1e-6)
#     # elif i==2:
#     #     fits_tof_pions.At(i).SetMaximum(1e-7)
#     # elif i==3:
#     #     fits_tof_pions.At(i).SetMaximum(10)


#     fits_tof_pions.At(i).Draw()

# canvas_fit.Update()







# input("wait")




