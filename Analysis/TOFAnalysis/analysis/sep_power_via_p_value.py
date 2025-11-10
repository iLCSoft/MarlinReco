#!/usr/bin/env python3

"""Produces separation power plots"""

import numpy as np
import ROOT
from utils import *

ROOT.EnableImplicitMT()

#Linear momentum bins
N_MOMENTUM_BINS, MIN_MOMENTUM, MAX_MOMENTUM = 50, 0, 18 # GeV/c
#Log momentum bins for dE/dx
N_LOG_MOMENTUM_BINS, MIN_LOG_MOMENTUM, MAX_LOG_MOMENTUM = 30, -0.3, 1.3 # (0.5 to ~20) GeV/c
LOG_MOMENTUM_BINS = np.array([ 10**(MIN_LOG_MOMENTUM + (MAX_LOG_MOMENTUM-MIN_LOG_MOMENTUM)*i/N_LOG_MOMENTUM_BINS ) for i in range(N_LOG_MOMENTUM_BINS+1) ])
N_DEDX_BINS, MIN_DEDX, MAX_DEDX = 3000, 0, 1e-6
DEDX_BINS = np.array([ MIN_DEDX + (MAX_DEDX - MIN_DEDX)*i/N_DEDX_BINS for i in range(N_DEDX_BINS+1) ])
N_MASS2_BINS, MIN_MASS2, MAX_MASS2 = 3000, -3, 3  # GeV^2/c^4
MOMENTUM_COLUMN = "harmonicMomToEcal_IKF_zedLambda"
TRACK_LENGTH_COLUMN = "trackLengthToEcal_IKF_zedLambda"
RESOLUTIONS = [0, 1, 5, 10, 17, 30, 50, 100] # ps
COLORS_RESOLUTION = [ ROOT.TColor.GetColor(c) for c in ["#00aaff", "#0091ea", "#0079d3", "#0061bd", "#004aa5", "#00348d", "#001d75", "#00045c"] ]
COLORS_DEDX = [ ROOT.TColor.GetColor(c) for c in ["#00aaff", "#cd54b9", "#c52700"] ]
MIN_SEP_POWER, MAX_SEP_POWER = -0.2, 7.5




def main():
    canvases = []
    df_init = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")

    # Get all histos first to utilise lazy RDataFrame behaviour
    histos = {"TOF" : {}, "dEdx" : {}}
    for RES in RESOLUTIONS:
        histos["TOF"][RES] = {}
        df = df_init.Filter("tofClosest0 > 6.").Define("beta", f"{TRACK_LENGTH_COLUMN}/(tofClosest{RES}*299.792458)")\
                    .Define("mass2", f"{MOMENTUM_COLUMN}*{MOMENTUM_COLUMN}*( 1./(beta*beta) - 1.)")
        for p in particles:
            histos["TOF"][RES][p] = {"lin": {}, "log" : {}}
            df_p = df.Filter(f"abs(pdg) == {p.pdg}")
            h_lin = df_p.Histo2D((get_rand_string(), "", N_MOMENTUM_BINS, MIN_MOMENTUM, MAX_MOMENTUM, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), MOMENTUM_COLUMN, "mass2")
            h_log = df_p.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_MASS2_BINS, MIN_MASS2, MAX_MASS2), MOMENTUM_COLUMN, "mass2")
            histos["TOF"][RES][p] = {"lin": h_lin, "log" : h_log}
    for p in particles:
        df_p = df_init.Filter("dEdx > 0.").Filter(f"abs(pdg) == {p.pdg}")
        # h = df_p.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_DEDX_BINS, MIN_DEDX, MAX_DEDX), MOMENTUM_COLUMN, "dEdx")
        h = df_p.Histo2D((get_rand_string(), "", N_LOG_MOMENTUM_BINS, LOG_MOMENTUM_BINS, N_DEDX_BINS, DEDX_BINS), MOMENTUM_COLUMN, "dEdx")
        histos["dEdx"][p] = h

    for tof_res in [10]:
        print(f"Calculating for dEdx for resolution {tof_res}")
        dedx_graphs = { "pik" : {}, "kp" : {} }
        _, _, _, dedx_graphs["pik"]["TOF"] = calculate_pid_graphs(histos["TOF"][tof_res][pion]["log"], histos["TOF"][tof_res][kaon]["log"])
        dedx_graphs["pik"]["TOF"].SetTitle(f"TOF {tof_res} ps;Momentum (GeV/c);#pi/K separation power")
        _, _, _, dedx_graphs["pik"]["dEdx"] = calculate_pid_graphs(histos["dEdx"][pion], histos["dEdx"][kaon])
        dedx_graphs["pik"]["dEdx"].SetTitle("dE/dx")
        c3 = draw_dedx_sep_powers(dedx_graphs["pik"]["TOF"], dedx_graphs["pik"]["dEdx"])
        c3.Modified()
        c3.Update()
        canvases.append(c3)

        _, _, _, dedx_graphs["kp"]["TOF"] = calculate_pid_graphs(histos["TOF"][tof_res][kaon]["log"], histos["TOF"][tof_res][proton]["log"])
        dedx_graphs["kp"]["TOF"].SetTitle(f"TOF {tof_res} ps;Momentum (GeV/c);K/p separation power")
        _, _, _, dedx_graphs["kp"]["dEdx"] = calculate_pid_graphs(histos["dEdx"][kaon], histos["dEdx"][proton])
        dedx_graphs["kp"]["dEdx"].SetTitle("dE/dx")
        c4 = draw_dedx_sep_powers(dedx_graphs["kp"]["TOF"], dedx_graphs["kp"]["dEdx"])
        c4.Modified()
        c4.Update()
        canvases.append(c4)
        input("wait")
    # resolution_graphs = { "pik" : {}, "kp" : {} }
    # for RES in RESOLUTIONS:
    #     print(f"Calculating for {RES} resolution")
    #     _, _, _, resolution_graphs["pik"][RES] = calculate_pid_graphs(histos["TOF"][RES][pion]["lin"], histos["TOF"][RES][kaon]["lin"])
    #     _, _, _, resolution_graphs["kp"][RES] = calculate_pid_graphs(histos["TOF"][RES][kaon]["lin"], histos["TOF"][RES][proton]["lin"])

    # c1 = draw_resolution_sep_powers(resolution_graphs["pik"])
    # resolution_graphs["pik"][0].SetTitle(";Momentum (GeV/c);#pi/K separation power")
    # c1.Modified()
    # c1.Update()
    # c2 = draw_resolution_sep_powers(resolution_graphs["kp"])
    # resolution_graphs["kp"][0].SetTitle(";Momentum (GeV/c);K/p separation power")
    # c2.Modified()
    # c2.Update()

    input("wait")

if __name__ == "__main__":
    main()