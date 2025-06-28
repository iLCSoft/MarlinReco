import ROOT
import numpy as np
from utils import *

ROOT.EnableImplicitMT()

df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")\
         .Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")
df = df.Filter("tofClosest0 > 6.")
df = df.Define("beta", f"trackLengthToEcal_IKF_zedLambda/(tofClosest0*299.792458)")
df = df.Define("mass2", "harmonicMomToEcal_IKF_zedLambda*harmonicMomToEcal_IKF_zedLambda*( 1./(beta*beta) - 1.)")

h = df.Histo2D(("h", "; Momentum (GeV/c); Mass^{2} (GeV^{2}/c^{4})", 1000, 0, 10, 1000, -0.3, 1.2), "harmonicMomToEcal_IKF_zedLambda","mass2")
canvas = draw_2d_plot(h, 12000)
input("wait")