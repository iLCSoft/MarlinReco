from utils import *
ROOT.gInterpreter.Declare('#include "tof.hpp"')

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root").Range(100)
df = df.Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")\
        .Filter("caloIDClosest == 1")\
        .Filter("tofClosest0 > 6.")\
        .Define("nHitsIn10Layers", "getNHitsInLayers(layerHit, 10)")\
        .Filter("nHitsIn10Layers > 0")\
        .Define("rImpact", "ROOT::Math::XYZVector(recoCaloX, recoCaloY, recoCaloZ)")\
        .Define("momImpact", "ROOT::Math::XYZVector(recoCaloPx, recoCaloPy, recoCaloPz)")\
        .Define("mom", "momImpact.r()")\
        .Define("hitPos", "hitPos(xHit, yHit, zHit)")\
        .Define("dToImpact", "dToImpact(hitPos, rImpact)")\
        .Define("dToLine", "dToLine(hitPos, rImpact, momImpact)")\
        .Define("dl", "sqrt(dToImpact*dToImpact - dToLine*dToLine)")\
        .Define("tSurface", "getTimeAtSurface(tHit, dToImpact)")

# LIMIT 10 LAYERS
df = df.Redefine("dl", f"dl[layerHit<{10}]")\
        .Redefine("dToImpact", f"dToImpact[layerHit<{10}]")\
        .Redefine("dToLine", f"dToLine[layerHit<{10}]")\
        .Redefine("tHit", f"tHit[layerHit<{10}]")\
        .Redefine("tSurface", f"tSurface[layerHit<{10}]")\
        .Redefine("layerHit", f"layerHit[layerHit<{10}]")

        
# R 10 mm CUT
df = df.Define("r_cut", "selectCylinderHits(dToLine, 10.)")\
        .Redefine("dl", f"dl[r_cut]")\
        .Redefine("dToImpact", f"dToImpact[r_cut]")\
        .Redefine("dToLine", f"dToLine[r_cut]")\
        .Redefine("layerHit", f"layerHit[r_cut]")\
        .Redefine("tHit", f"tHit[r_cut]")\
        .Redefine("tSurface", f"tSurface[r_cut]")

# # T_median 170 ps CUT
df = df.Define("t_cut", "selectMedianHits(tSurface, 170.)")\
        .Redefine("dl", f"dl[t_cut]")\
        .Redefine("dToImpact", f"dToImpact[t_cut]")\
        .Redefine("dToLine", f"dToLine[t_cut]")\
        .Redefine("tHit", f"tHit[t_cut]")\
        .Redefine("layerHit", f"layerHit[t_cut]")\
        .Redefine("tSurface", f"tSurface[t_cut]")

df = df.Define("tof_reco", "Mean(tSurface)")

data = df.AsNumpy(["dToImpact", "tHit", "mom", "tofClosest0", "pdg", "tof_reco"])
print(len(data["dToImpact"]))
print(len(data["tofClosest0"]))
print(len(data["tof_reco"]))

c = create_canvas(0.32, 0.8, 0.7)
frame = c.DrawFrame(-10,0,100.,20.)
frame.GetXaxis().SetTitle("d (mm)")
frame.GetYaxis().SetTitle("True T_{hit} (ns)")
frame.GetYaxis().SetTitleOffset(1.7)
for i, (d_arr,t_arr) in enumerate(zip(data["dToImpact"], data["tHit"])):
    tof_reco = data["tof_reco"][i]
    pdg = data["pdg"][i]
    mom = data["mom"][i]
    if mom > 1.:
        continue
    tof_true = data["tofClosest0"][i]
    if abs(pdg) == 211:
        mass = pion.mass
    elif abs(pdg) == 321:
        mass = kaon.mass
    elif abs(pdg) == 2212:
        mass = proton.mass
    beta = mom*np.sqrt(1/(mass*mass+mom*mom))

    frame.GetYaxis().SetRangeUser(tof_true-0.1, tof_true+0.3)

    d_arr = np.array([x for x in d_arr])    
    t_arr = np.array([x for x in t_arr])    
    gr = ROOT.TGraph(len(d_arr), d_arr, t_arr)
    gr.SetMarkerSize(2)
    gr.Draw("P")
    # frame.GetXaxis().SetRangeUser(-10, 100)
    # gr.GetYaxis().SetRangeUser(1.2)

    f_fit = ROOT.TF1("f_fit", "pol1", 0, 100)
    f_fit.SetLineColor(2)

    gr.Fit("f_fit", "QR")
    fit = gr.GetFunction("f_fit")
    tof_fit = fit.GetParameter("p0")
    f1 = ROOT.TF1("f1", "pol1", 0, 100)
    f1.SetParameters(tof_reco, 1./SPEED_OF_LIGHT)
    f1.SetLineColor(8)
    f1.Draw("same")

    gr_reco = ROOT.TGraph()
    gr_reco.SetPoint(0, 0., tof_reco)
    gr_reco.SetMarkerStyle(29)
    gr_reco.SetMarkerColor(8)
    gr_reco.SetMarkerSize(2)
    gr_reco.SetTitle("TOF reco;d (mm); T_{hit} (ns)")
    gr_reco.Draw("Psame")

    gr_fit = ROOT.TGraph()
    gr_fit.SetPoint(0, 0., tof_fit)
    gr_fit.SetMarkerStyle(29)
    gr_fit.SetMarkerColor(2)
    gr_fit.SetMarkerSize(2)
    gr_fit.SetTitle("TOF fit;d (mm); T_{hit} (ns)")
    gr_fit.Draw("Psame")

    gr_true = ROOT.TGraph()
    gr_true.SetPoint(0, 0., tof_true)
    gr_true.SetMarkerStyle(29)
    gr_true.SetMarkerColor(ROOT.kYellow+1)
    gr_true.SetMarkerSize(2)
    gr_true.SetTitle("TOF true;d (mm); T_{hit} (ns)")
    gr_true.Draw("Psame")

    print("pdg: ", pdg)
    print("mom: ", mom)
    print("tof true: ", tof_true)
    print("tof reco: ", tof_reco)
    print("tof_fit:", tof_fit)
    print("beta_true: ", beta)
    print("beta_reco: 1")
    print("beta_fit: ", 1./fit.GetParameter("p1")/SPEED_OF_LIGHT)
    print("1/c", 1./SPEED_OF_LIGHT)

    latex.DrawLatexNDC(0.6, 0.5, f"PDG: {pdg}")
    latex.DrawLatexNDC(0.6, 0.44, "p: "+f"{mom:.1f}"+" GeV/c")
    latex.DrawLatexNDC(0.6, 0.38, "#beta_{true}: "+f"{beta:.3f}")
    latex.DrawLatexNDC(0.6, 0.32, "#beta_{fit}: "+f"{1./fit.GetParameter(1)/SPEED_OF_LIGHT:.3f}")

    c.Update()
    input("wait")

# h_test = df_test.Histo1D((get_rand_string(), "; T_{reco} - T_{true} (ps); N entries", 1500, -300, 300), "dt")
# arr_test = df_test.AsNumpy(["momImpact" ,"tHit" ,"dToImpact" ,"pdg"])
