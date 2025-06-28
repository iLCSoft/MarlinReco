import ROOT
ROOT.EnableImplicitMT()

class Algorithm:
    def __init__(self, name, momentum, track_length, tof):
        self.name = name
        self.mom = momentum
        self.len = track_length
        self.tof = tof

# colors = [ROOT.TColor.GetColor('#ff7f00') ,ROOT.TColor.GetColor('#984ea3') ,ROOT.TColor.GetColor('#4daf4a') ,ROOT.TColor.GetColor('#377eb8') ,ROOT.TColor.GetColor('#e41a1c')]
colors_3_qualitative = [ROOT.TColor.GetColor('#d95f02'), ROOT.TColor.GetColor('#7570b3'), ROOT.TColor.GetColor('#1b9e77')]
colors = colors_3_qualitative



algoClosest30 = Algorithm("Single ECAL layer", "harmonicMomToEcal", "trackLengthToEcal", "tofClosest30")
algoSET50 = Algorithm("Double SET layer", "harmonicMomToSET", "trackLengthToSET", "tofSETFront50")
algoAverage100 = Algorithm("Ten ECAL layers", "harmonicMomToEcal", "trackLengthToEcal", "tofAverage100")
algorithms = [algoClosest30, algoSET50, algoAverage100]

def draw_lines_1d(maxy):
    lines = {}
    pdgs = [211, 321, 2212]
    m_pdg = {211 : 0.13957039*1000, 321 : 0.493677*1000, 2212 : 0.938272088*1000}
    for pdg in pdgs:
        lines[pdg] = ROOT.TLine(m_pdg[pdg], 0., m_pdg[pdg], maxy)
        lines[pdg].SetLineColor(15)
        lines[pdg].SetLineWidth(2)
        lines[pdg].SetLineStyle(9)
        lines[pdg].Draw()
    return lines


df = ROOT.RDataFrame("treename", "/nfs/dust/ilc/user/dudarboh/tof/BohdanAna.root")
df = df.Filter("tofClosest0 > 6. && tofSETFront0 > 0.").Filter("abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 2212")

histos = []
for alg in algorithms:
    # NOTE beta quality cut significantly distords TOF gaussian distributions! (makes sense...)
    # df_beta = df.Define("beta", f"{alg.len}/({alg.tof}*299.792458)").Filter("beta >= 0 && beta <= 1")
    tof_true = "tofClosest0" if "SET" not in alg.tof else "tofSETFront0"
    df_res = df.Define("res", f"({alg.tof} - {tof_true})*1000")
    h = df_res.Histo1D((f"{alg.name}", f"{alg.name}; reconstructed TOF - true TOF [ps]; N entries", 1500, -300, 300), "res")
    histos.append(h)


canvas = ROOT.TCanvas()
canvas.SetLeftMargin(0.18)
canvas.SetRightMargin(0.05)
canvas.SetGridx(True)
maxy = 0
for i, h in enumerate(histos):
    h.GetYaxis().SetTitleOffset(1.65)
    h.SetLineColor(colors[i])
    h.SetMarkerColor(colors[i])
    h.SetMarkerStyle(0)
    h.SetLineWidth(3)
    maxy = max( maxy, 1.25*h.GetMaximum() )


histos[0].SetMaximum(maxy)
histos[0].Draw() # draw axis for lines, god bless...
histos[0].GetXaxis().SetRangeUser(-290, 290)

for i, h in enumerate(histos):
    h.Draw("same")

fits = []

for i, h in enumerate(histos):
    f = ROOT.TF1(f"f_{i}", "gaus", -50, 50)
    f.SetLineColor(colors[i])
    f.SetLineWidth(2)
    h.Fit(f, "RN")
    f.SetNpx(1000)
    fits.append(f)
    f.Draw("same")

legend = ROOT.TLegend(0.6, 0.44, 0.996, 0.935)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetMargin(0.1)

entry1 = legend.AddEntry(histos[0].GetPtr(), histos[0].GetTitle(), "l")
entry2 = legend.AddEntry("", "#sigma_{hit} = 30 ps", "")
entry3 = legend.AddEntry("", "#sigma_{TOF, fit} = " + f"{round(fits[0].GetParameter(2), 1)} ps", "")
gap_entry1 = legend.AddEntry("", "", "")

entry4 = legend.AddEntry(histos[1].GetPtr(), histos[1].GetTitle(), "l")
entry5 = legend.AddEntry("", "#sigma_{hit} = 50 ps", "")
entry6 = legend.AddEntry("", "#sigma_{TOF, fit} = " + f"{round(fits[1].GetParameter(2), 1)} ps", "")
gap_entry2 = legend.AddEntry("", "", "")

entry7 = legend.AddEntry(histos[2].GetPtr(), histos[2].GetTitle(), "l")
entry8 = legend.AddEntry("", "#sigma_{hit} = 100 ps", "")
entry9 = legend.AddEntry("", "#sigma_{TOF, fit} = " + f"{round(fits[2].GetParameter(2), 1)} ps", "")
gap_entry3 = legend.AddEntry("", "", "")

entry1.SetTextColor(colors[0])
entry2.SetTextColor(colors[0])
entry3.SetTextColor(colors[0])
entry4.SetTextColor(colors[1])
entry5.SetTextColor(colors[1])
entry6.SetTextColor(colors[1])
entry7.SetTextColor(colors[2])
entry8.SetTextColor(colors[2])
entry9.SetTextColor(colors[2])

legend.Draw()

latex = ROOT.TLatex()
latex.SetTextFont(52)
latex.DrawLatex(64, 28750., "ILD preliminary")

canvas.Update()
input("Finish")
