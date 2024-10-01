import ROOT
import numpy as np
from utils import *

def bias_from_fit():
    m2_true = np.array([p.mass*1000 for p in particles])
    m2_rec = np.array([139.81318714819113, 495.17769543103486, 941.0813860261262])
    m2_rec_set = np.array([141.06505909569208, 494.4509704678564, 939.0805896628319])
    # m2_rec_set = np.sqrt(m2_rec_set)

    diff = (m2_rec - m2_true)
    diff_set = (m2_rec_set - m2_true)

    canvas = create_canvas()
    frame1 = canvas.DrawFrame(0., -1.1*max(diff), 1., 1.1*max(diff), "; Mass_{true} (GeV/c^{2}); Mass_{reco} - Mass_{true} (MeV/c^{2})")
    frame1.GetYaxis().SetMaxDigits(3)
    gr = ROOT.TGraph(3, m2_true/1000., diff)
    gr_set = ROOT.TGraph(3, m2_true/1000., diff_set)
    gr.Draw("PL")
    gr_set.Draw("PLsame")
    gr_set.SetLineColor(6)
    gr_set.SetMarkerColor(6)
    legend = create_legend()
    legend.AddEntry(gr, "ECAL surface", "pl")
    legend.AddEntry(gr_set, "SET", "pl")
    legend.Draw()
    canvas.Update()
    input("wait")

bias_from_fit()

