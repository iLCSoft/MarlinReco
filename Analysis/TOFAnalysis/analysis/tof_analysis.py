import sys
import ROOT
import numpy as np
from sklearn import linear_model
from utils import *
import time
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetNumberContours(256)
ROOT.gStyle.SetOptTitle(0)
ROOT.EnableImplicitMT()
ROOT.gInterpreter.Declare('#include "tof.hpp"')

def create_tof_studies_file():
    '''Create tof_studies.root that containes stripped information from BohdanAna specifically for TOF reconstruction studies'''
    df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
    # Filtering
    df = df.Filter("if (rdfentry_ % 1000000 == 0){ std::cout << rdfentry_ << std::endl; } return true;")\
            .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")\
            .Filter("caloIDClosest == 1")\
            .Filter("tofClosest0 > 6.")\
            .Define("nHitsIn10Layers", "getNHitsInLayers(layerHit, 10)")\
            .Filter("nHitsIn10Layers > 0")

    # Defining 3D vectors
    df = df.Define("rImpact", "ROOT::Math::XYZVector(recoCaloX, recoCaloY, recoCaloZ)")\
            .Define("momImpact", "ROOT::Math::XYZVector(recoCaloPx, recoCaloPy, recoCaloPz)")\
            .Define("hitPos", "hitPos(xHit, yHit, zHit)")

    df = df.Define("dToImpact", "dToImpact(hitPos, rImpact)")\
            .Define("dToLine", "dToLine(hitPos, rImpact, momImpact)")\
            .Define("dl", "sqrt(dToImpact*dToImpact - dToLine*dToLine)")

    write_columns = ["tHit", "layerHit", "dToImpact", "dToLine", "dl", "pdg", "tofClosest0"]
    df.Snapshot("treename", "/nfs/dust/ilc/user/dudarboh/tof/tof_studies_UPDATE_ME.root", write_columns)
    sys.exit("Done!")

def filter_n_layers(df, n_layers):
    '''Return df with hits only within first n_layers'''
    df = df.Redefine("dl", f"dl[layerHit<{n_layers}]")\
            .Redefine("dToImpact", f"dToImpact[layerHit<{n_layers}]")\
            .Redefine("dToLine", f"dToLine[layerHit<{n_layers}]")\
            .Redefine("tHit", f"tHit[layerHit<{n_layers}]")\
            .Redefine("layerHit", f"layerHit[layerHit<{n_layers}]")
    return df

def smear_time(df):
    '''Return df with time smeared with time_resolution (in ps)'''
    # Only with 50 ps for now
    df = df.Redefine("tHit", "smear(rdfslot_, tHit, gaus50)")\
            .Define("tSurface", "getTimeAtSurface(tHit, dToImpact)")
    return df

def scan_cuts(df):
    '''Return numpy array of the dt results using scan cuts'''
    d_perp_cuts = np.linspace(4, 14, 20)
    dt_cuts = np.linspace(100, 250, 20)
    results = []
    for i, d_perp_cut in enumerate(d_perp_cuts):
        for j, dt_cut in enumerate(dt_cuts):

            results.append(f"dt_{i}_{j}")
            df = df.Define(f"tSurface_cyl_mask_{i}_{j}", f"tSurface[selectCylinderHits(dToLine, {d_perp_cut})]")\
            .Define(f"tSurface_both_masks_{i}_{j}", f"tSurface_cyl_mask_{i}_{j}[selectMedianHits(tSurface_cyl_mask_{i}_{j}, {dt_cut})]")\
            .Define(f"dt_{i}_{j}", f"1000*(Mean(tSurface_both_masks_{i}_{j}) - tofClosest0)")
    return d_perp_cuts, dt_cuts, df.AsNumpy(results)



t = time.time()
df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/tof_studies.root")\
    .Filter("if (rdfentry_ % 1000000 == 0){ std::cout << rdfentry_ << std::endl; } return true;")\
    .Filter("rdfentry_ < 2000000")
df = filter_n_layers(df, 10)
df = smear_time(df)


canvas = create_canvas()
### ALL HITS
# df_test = df.Define("dt", "1000*(Mean( tSurface ) - tofClosest0)")
### FRANK SELECTION
# df_test = df.Define("dt", "1000*(Mean( tSurface[selectFrankHits(dToLine, layerHit)]) - tofClosest0)")
### NEW SELECTION
# df_test = df.Define("tSurface_cyl_mask", "tSurface[selectCylinderHits(dToLine, 10.)]")\
#             .Define("tSurface_both_masks", "tSurface_cyl_mask[selectMedianHits(tSurface_cyl_mask, 170.)]")\
#             .Define("dt", "1000*(Mean(tSurface_both_masks) - tofClosest0)")

# h_test = df_test.Histo1D((get_rand_string(), "; T_{reco} - T_{true} (ps); N entries", 1500, -300, 300), "dt")
# arr_test = df_test.AsNumpy(["dt"])

### CONRAD NN RESULTS
# arr_test = {"dt": np.load("./konrad_ml/konrad_nn_tof_reco.npy") - np.load("./konrad_ml/konrad_nn_tof_true.npy")}
# h_test = ROOT.TH1D(get_rand_string(), "; T_{reco} - T_{true} (ps); N entries", 1500, -300, 300)
# h_test.FillN(len(arr_test["dt"]), arr_test["dt"], np.ones_like(arr_test["dt"]) )

### FIT ATTEMPT ON NEW SELECTION
df_test = df.Define("tSurface_cyl_mask", "tSurface[selectCylinderHits(dToLine, 10.)]")\
            .Define("tSurface_both_masks", "tSurface_cyl_mask[selectMedianHits(tSurface_cyl_mask, 170.)]")\
            .Define("dt", "1000*(Mean(tSurface_both_masks) - tofClosest0)")

h_test = df_test.Histo1D((get_rand_string(), "; T_{reco} - T_{true} (ps); N entries", 1500, -300, 300), "dt")
arr_test = df_test.AsNumpy(["dt"])

_, rms90, _, _ = fit90( arr_test["dt"] )
std = np.std(arr_test["dt"])
print(" above 50 ps:", np.sum(arr_test["dt"] > 50.)/len(arr_test["dt"]) )
print("RMS:", std)
print("RMS90:", rms90)
h_test.Draw()
h_test.GetYaxis().SetMaxDigits(3)
h_test.GetXaxis().SetTitleOffset(1)
fit = ROOT.TF1("fit", "gaus", -20., 20.)
fit.SetNpx(500)
h_test.Fit("fit", "QR0")
fit.SetRange(-300, 300)
fit.Draw("same")
std_fit = h_test.GetFunction("fit").GetParameter(2)
print("RMS fit:", std_fit)
latex.DrawLatexNDC(0.23, 0.8, "RMS_{total}: " + f"{std:.1f} ps")
latex.DrawLatexNDC(0.23, 0.73, "RMS_{90}: " + f"{rms90:.1f} ps")
latex.DrawLatexNDC(0.23, 0.66, "RMS_{fit}: " + f"{std_fit:.1f} ps")
input("waut")
#################################################### COMPARE FRANK VS NEW

# canvas = create_canvas()
# df_frank = df.Define("dt", "1000*(Mean( tSurface[selectFrankHits(dToLine, layerHit)]) - tofClosest0)")
# h_frank = df_frank.Histo1D((get_rand_string(), "; T_{reco} - T_{true} (ps); N entries", 1500, -300, 300), "dt")
# df_new = df.Define("tSurface_cyl_mask", "tSurface[selectCylinderHits(dToLine, 10.)]")\
#             .Define("tSurface_both_masks", "tSurface_cyl_mask[selectMedianHits(tSurface_cyl_mask, 170.)]")\
#             .Define("dt", "1000*(Mean(tSurface_both_masks) - tofClosest0)")
# h_new = df_new.Histo1D((get_rand_string(), "; T_{reco} - T_{true} (ps); N entries", 1500, -300, 300), "dt")
# arr_frank = df_frank.AsNumpy(["dt"])
# arr_new = df_new.AsNumpy(["dt"])


# _, rms90, _, _ = fit90( arr_frank["dt"] )
# print("RMS:", np.std(arr_frank["dt"]))
# print("RMS90:", rms90)
# h_frank.Draw()
# h_frank.GetYaxis().SetMaxDigits(3)
# h_frank.GetXaxis().SetTitleOffset(1)
# fit = ROOT.TF1("fit", "gaus", -20., 20.)
# fit.SetNpx(500)
# h_frank.Fit("fit", "QR0")
# fit.SetRange(-300, 300)
# fit.Draw("same")
# print("RMS fit:", h_frank.GetFunction("fit").GetParameter(2))

# _, rms90, _, _ = fit90( arr_new["dt"] )
# print("RMS new:", np.std(arr_new["dt"]))
# print("RMS90 new:", rms90)
# h_new.Draw("same")
# h_new.SetLineColor(4)
# h_new.GetYaxis().SetMaxDigits(3)
# h_new.GetXaxis().SetTitleOffset(1)
# fit_new = ROOT.TF1("fit_new", "gaus", -20., 20.)
# fit_new.SetNpx(500)
# h_new.Fit("fit_new", "QR0")
# fit_new.SetRange(-300, 300)
# fit_new.SetLineColor(6)
# fit_new.Draw("same")
# print("RMS fit_new new:", h_new.GetFunction("fit_new").GetParameter(2))

# input("waut")



############################################ GET THE SCAN PLOT #################################

# print("Scanning...")
# d_perp_cuts, dt_cuts, results = scan_cuts(df)

# print("Plotting...")
# gr = ROOT.TGraph2D()
# scan_results_for_konrad = []
# for i, d_perp_cut in enumerate(d_perp_cuts):
#     for j, dt_cut in enumerate(dt_cuts):
#         _, rms90, _, _ = fit90( results[f"dt_{i}_{j}"] )
#         scan_results_for_konrad.append([d_perp_cut, dt_cut, rms90])
#         gr.SetPoint( i*len(d_perp_cuts) + j, d_perp_cut, dt_cut, rms90 )

# np.save("scan.npy", np.array(scan_results_for_konrad) )

# canvas = create_canvas(0.33, 0.58, 0.65)
# gr.Draw("colz")
# canvas.Update()
# h = gr.GetHistogram()
# h.GetXaxis().SetTitle("R (mm)")
# h.GetYaxis().SetTitle("T_{cut} (ps)")
# h.GetXaxis().SetTitleOffset(1.1)
# h.GetYaxis().SetTitleOffset(1.4)
# palette = h.GetListOfFunctions().FindObject("palette")
# x_min = h.GetXaxis().GetXmin()
# x_max = h.GetXaxis().GetXmax()
# palette.SetX1(x_min + 1.01*(x_max-x_min))
# palette.SetX2(x_min + 1.05*(x_max-x_min))
# palette.SetMaxDigits(3)
# palette.SetLabelOffset(0.006)

# canvas.Modified()
# canvas.Update()

# minbin = h.GetMinimumBin()
# for i, d_perp_cut in enumerate(d_perp_cuts):
#     for j, dt_cut in enumerate(dt_cuts):
#         if h.FindBin(d_perp_cut, dt_cut) == minbin:
#             print("R: ", d_perp_cut, "    dt: ", dt_cut)
# print( h.GetMinimum())
# print("Time:", time.time() - t)
# input("wait")