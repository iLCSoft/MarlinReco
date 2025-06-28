from utils import *
ROOT.gStyle.SetPalette(ROOT.kBird)
ROOT.gStyle.SetNumberContours(256)
ROOT.EnableImplicitMT()

# NOTE: TIMING CUT IS PRESENT !!!!!!!
df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
        .Filter("tofClosest0 > 6.")\
        .Define("pTrue", "sqrt(mcPx*mcPx + mcPy*mcPy + mcPz*mcPz)")\
        .Filter("pTrue < 2.")


def refit_impact(df, variable="omega", particle=kaon):
    if particle == kaon:
        df = df.Filter("abs(pdg) == 321")
    elif particle == proton:
        df = df.Filter("abs(pdg) == 2212")
    else:
        return 0

    if variable == "omega":
        title_default = ";Momentum (GeV/c);#Omega_{reco} - #Omega_{true} (1/mm)"
        title_refitted = ";Momentum (GeV/c);#Omega_{reco} - #Omega_{true} (1/mm)"
        title_pull = ";(_{}#Omega_{reco} - #Omega_{true}) / #sigma_{#Omega};N entries"
        n_x_bins, x_min, x_max, n_y_bins, y_min, y_max =  200, 0, 2, 200, -0.0001, 0.0001

        y_column_default = "omegaIP - omegaTrue"
        y_column_refitted = "refittedOmegaIP - omegaTrue"
        y_column_default_pull = "(omegaIP - omegaTrue)/omegaErrIP"
        y_column_refitted_pull = "(refittedOmegaIP - omegaTrue)/refittedOmegaErrIP"
    elif variable == "tanLambda":
        title_default = ";Momentum (GeV/c);tan#lambda_{reco} - tan#lambda_{true}"
        title_refitted = ";Momentum (GeV/c);tan#lambda_{reco} - tan#lambda_{true}"
        title_pull = ";(_{}tan#lambda_{reco} - tan#lambda_{true}) / #sigma_{tan#lambda};N entries"
        n_x_bins, x_min, x_max, n_y_bins, y_min, y_max =  200, 0, 2, 200, -0.01, 0.01

        y_column_default = "tanLambdaIP - tanLambdaTrue"
        y_column_refitted = "refittedTanLambdaIP - tanLambdaTrue"
        y_column_default_pull = "(tanLambdaIP - tanLambdaTrue)/tanLambdaErrIP"
        y_column_refitted_pull = "(refittedTanLambdaIP - tanLambdaTrue)/refittedTanLambdaErrIP"
    elif variable == "d0":
        title_default = ";Momentum (GeV/c);d0_{reco} - d0_{true} (mm)"
        title_refitted = ";Momentum (GeV/c);d0_{reco} - d0_{true} (mm)"
        title_pull = ";(_{}d0_{reco} - d0_{true}) / #sigma_{d0};N entries"
        n_x_bins, x_min, x_max, n_y_bins, y_min, y_max =  200, 0, 2, 200, -0.3, 0.3

        y_column_default = "d0IP - d0True"
        y_column_refitted = "refittedD0IP - d0True"
        y_column_default_pull = "(d0IP - d0True)/d0ErrIP"
        y_column_refitted_pull = "(refittedD0IP - d0True)/refittedD0ErrIP"
    elif variable == "z0":
        title_default = ";Momentum (GeV/c);z0_{reco} - z0_{true} (mm)"
        title_refitted = ";Momentum (GeV/c);z0_{reco} - z0_{true} (mm)"
        title_pull = ";(_{}z0_{reco} - z0_{true}) / #sigma_{z0};N entries"
        n_x_bins, x_min, x_max, n_y_bins, y_min, y_max =  200, 0, 2, 200, -0.3, 0.3

        y_column_default = "z0IP - z0True"
        y_column_refitted = "refittedZ0IP - z0True"
        y_column_default_pull = "(z0IP - z0True)/z0ErrIP"
        y_column_refitted_pull = "(refittedZ0IP - z0True)/refittedZ0ErrIP"
    elif variable == "phi":
        title_default = ";Momentum (GeV/c);#varphi_{reco} - #varphi_{true}"
        title_refitted = ";Momentum (GeV/c);#varphi_{reco} - #varphi_{true}"
        title_pull = ";(_{}#varphi_{reco} - #varphi_{true}) / #sigma_{#varphi};N entries"
        n_x_bins, x_min, x_max, n_y_bins, y_min, y_max =  200, 0, 2, 200, -0.01, 0.01

        y_column_default = "phiIP - phiTrue"
        y_column_refitted = "refittedPhiIP - phiTrue"
        y_column_default_pull = "(phiIP - phiTrue)/phiErrIP"
        y_column_refitted_pull = "(refittedPhiIP - phiTrue)/refittedPhiErrIP"



    h_default = df.Define("diff", y_column_default).Histo2D((get_rand_string(), title_default, n_x_bins, x_min, x_max, n_y_bins, y_min, y_max), "pTrue", "diff")
    h_refitted = df.Define("diff", y_column_refitted).Histo2D((get_rand_string(), title_refitted, n_x_bins, x_min, x_max, n_y_bins, y_min, y_max), "pTrue", "diff")
    h_default_pull = df.Define("diff", y_column_default_pull).Histo1D((get_rand_string(), title_pull, 200, -5, 5), "diff")
    h_refitted_pull = df.Define("diff", y_column_refitted_pull).Histo1D((get_rand_string(), title_pull, 200, -5, 5), "diff")
    print( h_default.ProjectionY().GetStdDev() )
    print( h_refitted.ProjectionY().GetStdDev() )

    c1 = draw_2d_plot(h_default, 0.4, 0.65, 0.65)
    c2 = draw_2d_plot(h_refitted, 0.4, 0.65, 0.65)

    h_default.GetYaxis().SetTitleOffset(1.9)
    h_refitted.GetYaxis().SetTitleOffset(1.9)
    h_default.GetYaxis().SetMaxDigits(3)
    h_refitted.GetYaxis().SetMaxDigits(3)

    h_default.SetMinimum(1)
    h_refitted.SetMinimum(1)
    common_max = 1.05*max(h_default.GetMaximum(), h_refitted.GetMaximum())
    h_default.SetMaximum(common_max)
    h_refitted.SetMaximum(common_max)
    c1.Update()
    c2.Update()

    c3 = create_canvas()
    h_default_pull.SetTitle(f"Default  Std Dev: {h_default_pull.GetStdDev():.3f}")
    h_refitted_pull.SetTitle(f"Refitted  Std Dev: {h_refitted_pull.GetStdDev():.3f}")
    h_refitted_pull.SetLineColor(ROOT.kRed)
    h_default_pull.GetXaxis().SetTitleOffset(1.)
    h_default_pull.GetYaxis().SetMaxDigits(3)
    h_default_pull.SetMinimum( 0. )
    h_default_pull.SetMaximum( 1.3*max(h_default_pull.GetMaximum(), h_refitted_pull.GetMaximum()) )
    h_default_pull.Draw()
    h_refitted_pull.Draw("same")
    leg = create_legend(0.19, 0.77, 0.999, 0.91)
    leg.AddEntry(h_default_pull.GetPtr(), h_default_pull.GetTitle(), "l")
    leg.AddEntry(h_refitted_pull.GetPtr(), h_refitted_pull.GetTitle(), "l")
    leg.Draw()
    latex.DrawLatexNDC(0.23, 0.67, "p < 2 GeV/c")
    c3.Update()

    c1.SetName(f"refit_{particle.short_name}_{variable}_default")
    c2.SetName(f"refit_{particle.short_name}_{variable}_refitted")
    c3.SetName(f"refit_{particle.short_name}_{variable}_pulls")

    c1.Print("", ".png")
    c2.Print("", ".png")
    c3.Print("", ".png")

    input("wait")



refit_impact(df, variable="phi", particle=proton)