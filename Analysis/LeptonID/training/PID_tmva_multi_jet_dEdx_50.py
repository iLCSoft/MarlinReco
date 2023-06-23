#!/usr/bin/env python
# coding: utf-8

import ROOT
import argparse

parser = argparse.ArgumentParser(description="Train BDTG")
parser.add_argument("-n", "--ncores", type=int, default=2)
parser.add_argument("-t", "--ntrees", type=int, default=800)
parser.add_argument("-d", "--depth", type=int, default=3)
args = parser.parse_args()
ntrees = args.ntrees
depth = args.depth

ROOT.EnableImplicitMT(args.ncores)

# make input file an argument
dir = "../outputs/"
inFile = ROOT.TFile.Open(f"{dir}/ZH_eLpR_test_bb.root")
inTree = inFile.Get("pfos")

# make the first part of the name an argument called basename
outfile = ROOT.TFile.Open(f"PID_tmva_multi_jet_dEdx_{ntrees}t_{depth}d_cm50r50.root", "RECREATE")

factory = ROOT.TMVA.Factory(f"PID_multi_jet_dEdx_{ntrees}t_{depth}d_cm50r50", outfile, "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=multiclass")

dataloader = ROOT.TMVA.DataLoader(f"PID_multi_jet_dEdx_{ntrees}t_{depth}d_cm50r50")

# TODO: clean-up and final variable selection!
# replace this by for-loop over a list
dataloader.AddVariable("e_over_p")
dataloader.AddVariable("ecal_share")
dataloader.AddVariable("seenEcalDep")
#dataloader.AddVariable("seenHcalDep") #red2
dataloader.AddVariable("seenYokeDep")
dataloader.AddVariable("shape0")
dataloader.AddVariable("shape1")
#dataloader.AddVariable("shape2") #red3
#dataloader.AddVariable("shape3") #red1
#dataloader.AddVariable("shape4") #red1
#dataloader.AddVariable("shape6") #red1
#dataloader.AddVariable("shape7") #red1
dataloader.AddVariable("shape8")
#dataloader.AddVariable("shape9") #red3
#dataloader.AddVariable("shape11") #red1
#dataloader.AddVariable("shape12") #red4
#dataloader.AddVariable("shape13") #red1
#dataloader.AddVariable("shape14") #red5
#dataloader.AddVariable("shape15") #red3
#dataloader.AddVariable("shape16") #red4
dataloader.AddVariable("shape17") #red5 reintroduced red5.1
#dataloader.AddVariable("shape18") #red4
dataloader.AddVariable("shape19")
dataloader.AddVariable("shape20")
dataloader.AddVariable("shape21")
dataloader.AddVariable("shape22")
dataloader.AddVariable("shape23")
dataloader.AddVariable("cluEllipsoid_r1")
dataloader.AddVariable("cluEllipsoid_r2")
dataloader.AddVariable("cluEllipsoid_r3")
#dataloader.AddVariable("cluEllipsoid_vol") #red2
dataloader.AddVariable("cluEllipsoid_r_ave")
dataloader.AddVariable("cluEllipsoid_density")
dataloader.AddVariable("cluEllipsoid_eccentricity_T")
dataloader.AddVariable("cluEllipsoid_eccentricity_L")

dataloader.AddVariable("dEdxDist_e")

dataloader.AddSpectator("absTruePDG := abs(truePDG)")
dataloader.AddSpectator("absParentPDG := abs(parentPDG)")
dataloader.AddSpectator("seenP")
dataloader.AddSpectator("seenCosTheta")

dataloader.AddTree(inTree, "electrons")
dataloader.AddTree(inTree, "muons")
dataloader.AddTree(inTree, "hadrons")

dataloader.AddCut("abs(truePDG) == 11 && parentPDG != 0", "electrons")
# removes muons coming from the Z decay in ZH events
dataloader.AddCut("abs(truePDG) == 13 && parentPDG != 94 && parentPDG != 11 && parentPDG != 0", "muons")
dataloader.AddCut("(abs(truePDG) == 211 || abs(truePDG) == 321 || abs(truePDG) == 2212) && parentPDG != 0", "hadrons")

# FIXME: ideally the std::isfinite cuts should not be necessary...
dataloader.PrepareTrainingAndTestTree("std::isfinite(cluEllipsoid_density) && std::isfinite(cluEllipsoid_eccentricity_L) && std::isfinite(shape11)", "SplitMode=Random:NormMode=NumEvents:!V")

factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDTG",
    f"!H:!V:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost=True:BaggedSampleFraction=0.50:DoBoostMonitor=True:NTrees={ntrees}:MaxDepth={depth}")

#factory.OptimizeAllMethods()

factory.TrainAllMethods()

factory.TestAllMethods()

factory.EvaluateAllMethods()

outfile.Close()
