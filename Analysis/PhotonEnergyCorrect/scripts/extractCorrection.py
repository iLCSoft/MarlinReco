import ROOT
from array import array
from math import cos, log10, pi, sqrt, isnan
# ###########################################################
# this code calculates various photon corrections
#
# it takes as input the root file prepared by prepareInput.py
#   containing trees for photons of various energies
#
# as output it produces the parameters to be passed to
#     photonCorrectionProcessor via its xml steering file,
#   together with a pdf file with various plots for sanity checking:
#                   the user should verify that these look reasonable!
#
# Daniel Jeans, 2025
#
# ###########################################################

#
# the name of the input file
#
inname = "singlegammas.root"

# don't include photon samples below this energy in fits
#   (lower energies often difficult to fit)
fit_emin = 0.99

# which energy/phi/theta variables in tree to calculcate correction for?
# True  = use cluster parameters (ie not corrected)
# False = use PFO parameters (potentially already corrected,
#         if photonCorrectionProcessor was run in reconstruction
useClus = False

# derive corrections as function of MC (True) or reco (False) info?
useMC = False

# ########################################################################

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPadLeftMargin(0.15)

fin = ROOT.TFile(inname, "read")

# name of output pdf file
plname = 'extractCor_' + inname.replace('.root', '.pdf')

if useClus:
    eVar = "clEn"
    phifVar = "clPhiFold"
    thVar = "clTheta"
else:
    eVar = "pfoEn"
    phifVar = "pfoPhiFold"
    thVar = "pfoTheta"

# what trees are present in the root file?
inTrees = {}
for kk in fin.GetListOfKeys():
    obj = fin.Get(kk.GetName())
    if obj.IsA() == ROOT.TTree.Class() and 'gamTree' in obj.GetName():
        en = float(obj.GetName().split('_')[1].replace('p', '.'))
        inTrees[en] = obj

ROOT.gStyle.SetOptStat(0)

cc = ROOT.TCanvas()
cc.Print(plname + '[')

# =====================================
# gapCompensate_theta
# costh based correction
#   eCorr_costh.push_back( -0.090 ); # _costhCorr_gaus1_norm_const = pars[0];
#   eCorr_costh.push_back( 0.     ); # _costhCorr_gaus1_norm_logen = pars[1];
#   eCorr_costh.push_back( 0.235  ); # _costhCorr_gaus1_mean       = pars[2];
#   eCorr_costh.push_back( 0.0072 ); # _costhCorr_gaus1_sigm       = pars[3];
#   eCorr_costh.push_back( -0.036 ); # _costhCorr_gaus2_norm_const = pars[4];
#   eCorr_costh.push_back( 0.     ); # _costhCorr_gaus2_norm_logen = pars[5];
#   eCorr_costh.push_back( 0.588  ); # _costhCorr_gaus2_mean       = pars[6];
#   eCorr_costh.push_back( 0.0121 ); # _costhCorr_gaus2_sigm       = pars[7];
#   eCorr_costh.push_back( -0.042 ); # _costhCorr_gaus3_norm       = pars[8];
#   eCorr_costh.push_back( 0.774  ); # _costhCorr_gaus3_mean       = pars[9];
#   eCorr_costh.push_back( 0.009  ); # _costhCorr_gaus3_sigm       = pars[10]
#   eCorr_costh.push_back( 1.002  ); # _costhCorr_endcap_scale     = pars[11]
#
cosThFnBar = ROOT.TF1('cosThFnBar',
                      '[0]*(x<0.8? 1+[10]*x+gaus(1)+gaus(4)+gaus(7):[11])')
# =====================================


def gapCompensate_theta():

    print('energy correction: costh in barrel')

    cc.Clear()
    cc.Divide(5, 4)

    encthHists = {}
    ic = 1
    for kk, vv in sorted(inTrees.items()):
        cc.cd(ic)
        hn = "encth_" + str(kk)
        nsig = 3.
        if kk > 30:
            nsig = 6
        encth = ROOT.TH2F(hn, hn, 500, 0, 1,
                          200, max(0.001, 1.-nsig/sqrt(kk)), 1+nsig/sqrt(kk))
        if useMC:
            vv.Draw(eVar+"/"+str(kk)+":abs(mcCth)>>"+hn, eVar+">0")
            encth.GetXaxis().SetTitle('cosTheta [MC]')
        else:
            vv.Draw(eVar+"/"+str(kk)+":abs(cos(pfoTheta))>>"+hn, eVar+">0")
            encth.GetXaxis().SetTitle('cosTheta [PFO]')
        encth.GetYaxis().SetTitle(eVar+' / trueE')
        encthHists[kk] = encth
        encth.Draw()
        ic += 1
    cc.Print(plname)

    costhprofiles = {}

    fitPars = {}
    fitParErrs = {}

    ic = 1
    for kk, encth in sorted(encthHists.items()):
        cc.cd(ic)
        epx = encth.ProfileX()
        costhprofiles[kk] = epx
        epx.Draw()
        epx.GetYaxis().SetTitle(eVar+' / trueE')

        cosThFnBar.SetParameter(0, 1.)
        cosThFnBar.SetParameter(1, -0.1)
        cosThFnBar.SetParameter(2, 0.241)
        cosThFnBar.SetParLimits(2, 0.22, 0.25)

        cosThFnBar.SetParameter(3, 0.01)
        cosThFnBar.SetParLimits(3, 0.0, 0.05)

        cosThFnBar.SetParameter(4, -0.1)
        cosThFnBar.SetParameter(5, 0.588)
        cosThFnBar.SetParameter(6, 0.01)
        cosThFnBar.SetParameter(7, -0.1)
        cosThFnBar.SetParameter(8, 0.774)
        cosThFnBar.SetParameter(9, 0.01)

        cosThFnBar.SetParameter(10, 0.0)
        cosThFnBar.SetParameter(11, 1.)

        epx.Fit('cosThFnBar', 'q', '', 0, 0.95)

        fitPars[kk] = []
        fitParErrs[kk] = []
        for i in range(0, 12):
            fitPars[kk].append(epx.GetFunction('cosThFnBar').GetParameter(i))
            fitParErrs[kk].append(epx.GetFunction('cosThFnBar').GetParError(i))

        ic += 1
    cc.Print(plname)

    parGraphs = {}
    for iv in range(0, 12):
        x = array('d')
        y = array('d')
        dx = array('d')
        dy = array('d')
        for kk, vv in sorted(fitPars.items()):
            if kk > fit_emin:
                x.append(log10(kk))
                dx.append(0.)
                y.append(vv[iv])
                dy.append(fitParErrs[kk][iv])
        pg = ROOT.TGraphErrors(len(x), x, y, dx, dy)
        pg.SetName('barrel_theta_par'+str(iv))
        pg.SetTitle('barrel_theta_par'+str(iv))
        parGraphs[iv] = pg

# the energy dependence is not too strong
# so we assume fixed values indep of energy
# to choose value, fit to pol1 (on logE scale), and evalulate at 10 GeV.

    cc.Clear()
    cc.Divide(4, 3)
    ic = 1

    local_energyCorr_costh = {}
    for iv, pg in parGraphs.items():
        cc.cd(ic)
        pg.Draw("apl")
        pg.GetHistogram().GetXaxis().SetTitle('log10( photon energy [GeV] )')
        pg.GetHistogram().GetYaxis().SetTitle('fitted parameter')
        pg.Fit("pol1", "q")
        valAtTen = pg.GetFunction("pol1").Eval(1.)  # log10(10) GeV = 1.
        local_energyCorr_costh[iv] = valAtTen
        ic += 1
    cc.Print(plname)

    # "output" version in processor
    #   eCorr_costh.push_back( -0.090);# _costhCorr_gaus1_norm_const=pars[0];
    #   eCorr_costh.push_back( 0.    );# _costhCorr_gaus1_norm_logen=pars[1];
    #   eCorr_costh.push_back( 0.235 );# _costhCorr_gaus1_mean      =pars[2];
    #   eCorr_costh.push_back( 0.0072);# _costhCorr_gaus1_sigm      =pars[3];
    #   eCorr_costh.push_back( -0.036);# _costhCorr_gaus2_norm_const=pars[4];
    #   eCorr_costh.push_back( 0.    );# _costhCorr_gaus2_norm_logen=pars[5];
    #   eCorr_costh.push_back( 0.588 );# _costhCorr_gaus2_mean      =pars[6];
    #   eCorr_costh.push_back( 0.0121);# _costhCorr_gaus2_sigm      =pars[7];
    #   eCorr_costh.push_back( -0.042);# _costhCorr_gaus3_norm      =pars[8];
    #   eCorr_costh.push_back( 0.774 );# _costhCorr_gaus3_mean      =pars[9];
    #   eCorr_costh.push_back( 0.009 );# _costhCorr_gaus3_sigm      =pars[10]
    #   eCorr_costh.push_back( 1.002 );# _costhCorr_endcap_scale    =pars[11]

    out_energyCorr_costh = [0]*12
    out_energyCorr_costh[0] = local_energyCorr_costh[1]
    out_energyCorr_costh[1] = 0.  # no energy dep
    out_energyCorr_costh[2] = local_energyCorr_costh[2]
    out_energyCorr_costh[3] = local_energyCorr_costh[3]
    out_energyCorr_costh[4] = local_energyCorr_costh[4]
    out_energyCorr_costh[5] = 0.  # no energy dep
    out_energyCorr_costh[6] = local_energyCorr_costh[5]
    out_energyCorr_costh[7] = local_energyCorr_costh[6]
    out_energyCorr_costh[8] = local_energyCorr_costh[7]
    out_energyCorr_costh[9] = local_energyCorr_costh[8]
    out_energyCorr_costh[10] = local_energyCorr_costh[9]
    out_energyCorr_costh[11] = local_energyCorr_costh[11]

    # now set the parameters of cosThFnBar to the best pars (decided at 10 GeV)
    cosThFnBar.SetParameter(0, 1.)
    for i in range(1, 12):
        cosThFnBar.SetParameter(i, local_energyCorr_costh[i])
    cosThFnBar.SetParameter(10, 0.)

    cc.Clear()
    cosThFnBar.SetNpx(500)
    cosThFnBar.Draw()
    cosThFnBar.GetHistogram().GetXaxis().SetTitle('cos(theta)')
    cosThFnBar.GetHistogram().GetYaxis().SetTitle('corr. func. @ 10GeV')
    cc.Print(plname)

    return out_energyCorr_costh


# ################################
output_energyCorr_costh = gapCompensate_theta()
# ################################

# ========================================
# phi correction in barrel
#
# _phiBarrelCorr_pos_const
# _phiBarrelCorr_pos_logen
# _phiBarrelCorr_depth
# _phiBarrelCorr_width1
# _phiBarrelCorr_width2
#
phiFnBar = ROOT.TF1('phiFnBar',
                    '[0]*(1.0+[1]*exp(-pow(x-[2],2)/' +
                    '(2.*pow((x-[2]<0?[3]:[4]),2))))')
# ========================================


def gapCompensate_phiBarrel():

    print('energy correction: phi in barrel')

    cc.Clear()
    cc.Divide(5, 4)

    enphiHists = {}
    ic = 1
    for kk, vv in sorted(inTrees.items()):
        cc.cd(ic)
        hn = "enphi_" + str(kk)
        nsig = 3.
        if kk > 30:
            nsig = 10

        enphi = ROOT.TH2F(hn, hn, 500, 0, pi/4, 500,
                          max(0.001, 1.-nsig/sqrt(kk)), 1.+nsig/sqrt(kk))
        enphi.GetYaxis().SetTitle(eVar + ' / trueE')

        if useMC:
            vv.Draw(eVar+"/"+str(kk)+":abs(mcPhiFold)>>"+hn,
                    "abs(mcCth)<0.8 && "+eVar+">0")
            enphi.GetXaxis().SetTitle('phi(folded) [MC]')
        else:
            vv.Draw(eVar+"/"+str(kk)+":abs(pfoPhiFold)>>"+hn,
                    "abs(cos(pfoTheta))<0.8 && "+eVar+">0")
            enphi.GetXaxis().SetTitle('phi(folded) [PFO]')

        enphiHists[kk] = enphi
        enphi.Draw()
        ic += 1
    cc.Print(plname)

    cc.Clear()
    cc.Divide(5, 4)

    fitPars = {}
    fitParErrs = {}

    ic = 1
    for kk, enphi in sorted(enphiHists.items()):
        cc.cd(ic)
        epx = enphi.ProfileX()
        epx.Draw()
        epx.GetYaxis().SetTitle(eVar + ' / trueE')

        phiFnBar.SetParameter(0, 1)
        phiFnBar.SetParameter(1, -0.1)
        phiFnBar.SetParameter(2, 0.45)
        phiFnBar.SetParameter(3, 0.05)
        phiFnBar.SetParameter(4, 0.05)

        epx.Fit('phiFnBar', "q")

        fitPars[kk] = []
        fitParErrs[kk] = []
        for i in range(0, 5):
            fitPars[kk].append(epx.GetFunction('phiFnBar').GetParameter(i))
            fitParErrs[kk].append(epx.GetFunction('phiFnBar').GetParError(i))

        fitPars[kk][3] = abs(fitPars[kk][3])
        fitPars[kk][4] = abs(fitPars[kk][4])

        ic += 1
    cc.Print(plname)

    phiparGraphs = {}
    for iv in range(0, 5):
        x = array('d')
        y = array('d')
        dx = array('d')
        dy = array('d')
        for kk, vv in sorted(fitPars.items()):
            if kk > fit_emin:
                x.append(log10(kk))
                dx.append(0.)
                y.append(vv[iv])
                dy.append(fitParErrs[kk][iv])
        pg = ROOT.TGraphErrors(len(x), x, y, dx, dy)
        pg.SetName('barrel_phi_par'+str(iv))
        pg.SetTitle('barrel_phi_par'+str(iv))
        phiparGraphs[iv] = pg

    cc.Clear()
    cc.Divide(3, 2)
    ic = 1

    output_energyCorr_barrelPhi = {}

    local_barrelPhi = []

    for iv, pg in phiparGraphs.items():
        cc.cd(ic)
        pg.Draw("apl")
        pg.GetHistogram().GetXaxis().SetTitle('log10( photon energy [GeV] )')
        pg.GetHistogram().GetYaxis().SetTitle('fitted parameter')
        pg.Fit("pol1", "q")
        valAtTen = pg.GetFunction("pol1").Eval(1.)  # log10(10) GeV = 1.

        local_barrelPhi.append(valAtTen)

        # this is in the form that we will use in the
        # photonCorrectionProcessor steering
        if iv == 1:  # depth
            output_energyCorr_barrelPhi[2] = valAtTen
        elif iv == 2:  # this is the gap position, keep energy dependence
            output_energyCorr_barrelPhi[0] = \
                pg.GetFunction("pol1").GetParameter(0)
            output_energyCorr_barrelPhi[1] = \
                pg.GetFunction("pol1").GetParameter(1)
        elif iv == 3:  # width (-ve side)
            output_energyCorr_barrelPhi[3] = valAtTen
        elif iv == 4:  # width (+ve side)
            output_energyCorr_barrelPhi[4] = valAtTen
        ic += 1
    cc.Print(plname)

    for i in range(0, len(local_barrelPhi)):
        if isnan(local_barrelPhi[i]):
            local_barrelPhi[i] = 0.

    phiFnBar.SetParameter(0, local_barrelPhi[0])
    phiFnBar.SetParameter(1, local_barrelPhi[1])
    phiFnBar.SetParameter(2, local_barrelPhi[2])
    phiFnBar.SetParameter(3, local_barrelPhi[3])
    phiFnBar.SetParameter(4, local_barrelPhi[4])

    cc.Clear()
    phiFnBar.Draw()
    cc.Print(plname)

    return output_energyCorr_barrelPhi


# ################################
output_energyCorr_barrelPhi = gapCompensate_phiBarrel()
# ################################

# =====================================
# endcap gaps
#
# this corrects the gaps between ECAL endcap modules
#
# =====================================
endXEnd = ROOT.TF1('endXEnd',
                   '[0]*(1.0+[1]*exp(-pow(x-[2],2)/(2.*pow([3], 2)))+' +
                   '[4]*exp(-pow(x-[5],2)/(2.*pow([6], 2) ) ) )', 200, 1500)


def gapCompensate_endcap():

    print('energy correction: gaps in endcap')

    cc.Clear()
    cc.Divide(5, 4)

    enendxHists = {}
    ic = 1
    for kk, vv in sorted(inTrees.items()):
        cc.cd(ic)

        hn = "enendcap_" + str(kk)
        nsig = 3.
        if kk > 30:
            nsig = 8

        enphi = ROOT.TH2F(hn, hn, 500, 30, 1500, 500,
                          max(0.001, 1.-nsig/sqrt(kk)), 1.+nsig/sqrt(kk))
        enphi.GetYaxis().SetTitle(eVar + ' / trueE')

        if useMC:
            vv.Draw(eVar+"/"+str(kk)+":mcEndX>>"+hn,
                    "mcEndX>0.8 && "+eVar+">0")
            enphi.GetXaxis().SetTitle('endcap quardant dist [mm, MC]')
        else:
            vv.Draw(eVar+"/"+str(kk)+":pfoEndX>>"+hn,
                    "pfoEndX>0.8 && "+eVar+">0")
            enphi.GetXaxis().SetTitle('endcap quadrant dist [mm, PFO]')

        enendxHists[kk] = enphi
        enphi.Draw()

        ic += 1
    cc.Print(plname)

    fitPars = {}
    fitParErrs = {}

    ic = 1
    for kk, encth in sorted(enendxHists.items()):
        cc.cd(ic)
        epx = encth.ProfileX()

        epx.Draw()
        epx.GetYaxis().SetTitle(eVar + ' / trueE')

        endXEnd.SetParameter(0, 1.)
        endXEnd.SetParameter(1, -0.04)
        endXEnd.SetParameter(2, 850.)
        endXEnd.SetParameter(3, 15.)
        endXEnd.SetParameter(4, -0.04)
        endXEnd.SetParameter(5, 1250.)
        endXEnd.SetParameter(6, 15.)

        epx.Fit('endXEnd', "q")

        fitPars[kk] = []
        fitParErrs[kk] = []
        for i in range(0, 7):
            fitPars[kk].append(epx.GetFunction('endXEnd').GetParameter(i))
            fitParErrs[kk].append(epx.GetFunction('endXEnd').GetParError(i))

        fitPars[kk][3] = abs(fitPars[kk][3])
        fitPars[kk][6] = abs(fitPars[kk][6])

        if kk > fit_emin:
            epx.SetMinimum(0.9)
        else:
            epx.SetMinimum(0.7)

        ic += 1
    cc.Print(plname)

    endgapGraphs = {}
    for iv in range(0, 7):
        x = array('d')
        y = array('d')
        dx = array('d')
        dy = array('d')
        for kk, vv in sorted(fitPars.items()):
            if kk > fit_emin:
                x.append(log10(kk))
                dx.append(0.)
                y.append(vv[iv])
                dy.append(fitParErrs[kk][iv])
        pg = ROOT.TGraphErrors(len(x), x, y, dx, dy)
        pg.SetName('barrel_phi_par'+str(iv))
        pg.SetTitle('barrel_phi_par'+str(iv))
        endgapGraphs[iv] = pg

    cc.Clear()
    cc.Divide(3, 3)
    ic = 1

    local_energyCorr_endcap = {}
    for iv, pg in endgapGraphs.items():
        cc.cd(ic)
        pg.Draw("apl")
        pg.GetHistogram().GetXaxis().SetTitle('log10( photon energy [GeV] )')
        pg.GetHistogram().GetYaxis().SetTitle('fitted parameter')
        pg.Fit("pol1", "q")
        valAtTen = pg.GetFunction("pol1").Eval(1.)  # log10(10) GeV = 1.
        local_energyCorr_endcap[iv] = valAtTen
        ic += 1
    cc.Print(plname)

    endXEnd.SetParameter(0, 1.)
    for i in range(1, 7):
        endXEnd.SetParameter(i, local_energyCorr_endcap[i])

    out_energyCorr_endcap = {}
    out_energyCorr_endcap[0] = local_energyCorr_endcap[1]
    out_energyCorr_endcap[1] = local_energyCorr_endcap[2]
    out_energyCorr_endcap[2] = local_energyCorr_endcap[3]
    out_energyCorr_endcap[3] = local_energyCorr_endcap[4]
    out_energyCorr_endcap[4] = local_energyCorr_endcap[5]
    out_energyCorr_endcap[5] = local_energyCorr_endcap[6]

    cc.Clear()
    endXEnd.Draw()
    cc.Print(plname)

    return out_energyCorr_endcap


# ################################
output_energyCorr_endcap = gapCompensate_endcap()
# ################################


# =================================
# now we should try to apply all the above correcions
#  look at energy distribution after
#  then linearise the energy measurement
# =================================

# energy correction functions

def setEnThetaCorrFn():
    cosThFnBar.SetParameter(0, 1.)
    cosThFnBar.SetParameter(1, output_energyCorr_costh[0])
    cosThFnBar.SetParameter(2, output_energyCorr_costh[2])
    cosThFnBar.SetParameter(3, output_energyCorr_costh[3])
    cosThFnBar.SetParameter(4, output_energyCorr_costh[4])
    cosThFnBar.SetParameter(5, 0.)
    cosThFnBar.SetParameter(6, output_energyCorr_costh[7])
    cosThFnBar.SetParameter(7, output_energyCorr_costh[8])
    cosThFnBar.SetParameter(8, output_energyCorr_costh[9])
    cosThFnBar.SetParameter(9, output_energyCorr_costh[10])
    cosThFnBar.SetParameter(10, 0.)
    return


def getEnThetaCorr(costheta):
    setEnThetaCorrFn()
    return cosThFnBar.Eval(abs(costheta))


def setEnPhiCorrFn(en):
    phiFnBar.SetParameter(0, 1.)
    phiFnBar.SetParameter(1, output_energyCorr_barrelPhi[2])
    phiFnBar.SetParameter(2, output_energyCorr_barrelPhi[0] +
                          log10(en)*output_energyCorr_barrelPhi[1])
    phiFnBar.SetParameter(3, output_energyCorr_barrelPhi[3])
    phiFnBar.SetParameter(4, output_energyCorr_barrelPhi[4])
    return


def getEnPhiCorr(en, phiFold):
    setEnPhiCorrFn(en)
    return phiFnBar.Eval(phiFold)


def setEnEndcapCorrFn():
    endXEnd.SetParameter(0, 1.)
    endXEnd.SetParameter(1, output_energyCorr_endcap[0])
    endXEnd.SetParameter(2, output_energyCorr_endcap[1])
    endXEnd.SetParameter(3, output_energyCorr_endcap[2])
    endXEnd.SetParameter(4, output_energyCorr_endcap[3])
    endXEnd.SetParameter(5, output_energyCorr_endcap[4])
    endXEnd.SetParameter(6, output_energyCorr_endcap[5])
    return


def getEnEndcapCorr(en, endX):
    setEnEndcapCorrFn()
    return endXEnd.Eval(endX)


# get the corrected energy
def getEnCorr(en, costheta, phiFold, endX):
    corrFac = getEnThetaCorr(costheta)
    if abs(costheta) < 0.8:
        corrFac = corrFac * getEnPhiCorr(en, phiFold)
    else:
        corrFac = corrFac * getEnEndcapCorr(en, endX)
    return en / corrFac


def energy_linearise():

    print('linearise overall energy response')

    enHistos = {}

    x = array('d')
    y = array('d')
    dx = array('d')
    dy = array('d')

    cc.Clear()
    cc.Divide(5, 4)
    ic = 1
    # loop over input energies
    for kk, vv in sorted(inTrees.items()):

        cc.cd(ic)

        hen_orig = ROOT.TH1F("enOrig"+str(kk), "enOrig"+str(kk), 100,
                             max(0.001, kk-6*sqrt(kk)), kk+6*sqrt(kk))
        hen_corr = ROOT.TH1F("enCorr"+str(kk), "enCorr"+str(kk), 100,
                             max(0.001, kk-6*sqrt(kk)), kk+6*sqrt(kk))

        hen_orig.GetXaxis().SetTitle('reconstructed energy')
        hen_orig.GetYaxis().SetTitle('photons / bin')
        hen_corr.GetXaxis().SetTitle('reconstructed energy')
        hen_corr.GetYaxis().SetTitle('photons / bin')

        hen_orig.SetLineColor(1)
        hen_corr.SetLineColor(2)

        # loop over the tree entries
        for i in range(0, vv.GetEntries()):
            vv.GetEntry(i)

            if useClus:
                en = vv.clEn
            else:
                en = vv.pfoEn

            hen_orig.Fill(en)

            if useMC:
                en_corr = getEnCorr(en, vv.mcCth,
                                    vv.mcPhiFold, vv.mcEndX)
            else:
                en_corr = getEnCorr(en, cos(vv.pfoTheta),
                                    vv.pfoPhiFold, vv.pfoEndX)

            hen_corr.Fill(en_corr)

        hen_corr.Draw()
        hen_orig.Draw("same")
        ic += 1

        enHistos[kk] = (hen_orig, hen_corr)

        try:
            hen_corr.Fit('gaus', "q")
            mm = hen_corr.GetFunction('gaus').GetParameter(1)
            # gaus width, to set restricted range
            ss = hen_corr.GetFunction('gaus').GetParameter(2)
            hen_corr.Fit('gaus', 'q', '', mm-2*ss, mm+2*ss)
            mm = hen_corr.GetFunction('gaus').GetParameter(1)
            # uncertainty on mean
            ss = hen_corr.GetFunction('gaus').GetParError(2)
        except Exception:
            print('some fitting error...ignoring for now...')
            mm = 0.
            ss = 10.

        x.append(log10(kk))
        dx.append(0.)
        y.append(mm/kk)
        dy.append(ss/kk)

    cc.Print(plname)

    cc.Clear()

    pg = ROOT.TGraphErrors(len(inTrees), x, y, dx, dy)
    pg.Draw("apl")
    pg.GetHistogram().GetXaxis().SetTitle('log10( photon energy [GeV] )')
    pg.GetHistogram().GetYaxis().SetTitle('fitted parameter')
    pg.Fit("pol1", "q", "", log10(fit_emin), 3)

    energyCor_Lin = []
    energyCor_Lin.append(pg.GetFunction("pol1").GetParameter(0))
    energyCor_Lin.append(pg.GetFunction("pol1").GetParameter(1))

    cc.Print(plname)

    return energyCor_Lin


# ################################
energyCor_Linearise = energy_linearise()
# ################################


# =====================================
# now look at angle corrections
#
# reconstructed PFO direction is sometimes biased
#
# =====================================

# =====================================
# phi in the barrel
# =====================================

def angle_phiBarrel():

    print('phi correction: barrel')

    barphiphiHists = {}

    cc.Clear()
    cc.Divide(5, 4)
    ic = 1
    for kk, vv in sorted(inTrees.items()):
        cc.cd(ic)
        hn = "barphiphi_" + str(kk)
        yy = 0.03/sqrt(kk)
        phiphi = ROOT.TH2F(hn, hn, 300, 0, 3.14159/4., 200, -yy, yy)
        phiphi.GetYaxis().SetTitle(phifVar+' - truePhiFold [rad]')

        if useMC:
            vv.Draw("("+phifVar+"-mcPhiFold):mcPhiFold>>"+hn,
                    eVar+">0 && abs(mcCth)<0.8")
            phiphi.GetXaxis().SetTitle('phi [folded, MC]')

        else:
            vv.Draw("("+phifVar+"-mcPhiFold):pfoPhiFold>>"+hn,
                    eVar+">0 && abs(cos(pfoTheta))<0.8")
            phiphi.GetXaxis().SetTitle('phi [folded, PFO]')

        barphiphiHists[kk] = phiphi
        phiphi.Draw()
        ic += 1
    cc.Print(plname)

    fn_barphiphi = ROOT.TF1("fn_barphiphi",
                            "[0] + [1]*(exp(-0.5*pow( (x-[2])/[3], 2 ) ))" +
                            "+ [4]*sin(4*x) + [5]*sin(8*x)" +
                            "+ [6]*sin(12*x) + [7]*sin(16*x)")

    fitPars = {}
    fitParErrs = {}

    ic = 1
    for kk, phiphi in sorted(barphiphiHists.items()):
        cc.cd(ic)
        epx = phiphi.ProfileX()
        epx.Draw()
        epx.GetYaxis().SetTitle(phifVar+' - truePhi [rad]')

        # first fix the sharp peak pos and size
        fn_barphiphi.FixParameter(2, 0.43)
        fn_barphiphi.FixParameter(3, 0.02)

        epx.Fit("fn_barphiphi", "q")

        # then free them and refit
        fn_barphiphi.ReleaseParameter(2)
        fn_barphiphi.ReleaseParameter(3)

        epx.Fit("fn_barphiphi", "q")

        fitPars[kk] = []
        fitParErrs[kk] = []
        for i in range(0, 8):
            fitPars[kk].append(
                epx.GetFunction("fn_barphiphi").GetParameter(i))
            fitParErrs[kk].append(
                epx.GetFunction("fn_barphiphi").GetParError(i))
        ic += 1
    cc.Print(plname)

    cc.Clear()
    cc.Divide(3, 3)

    fn_barphiphi = ROOT.TF1("fn_phiphiFit",
                            "[0] + [1]/(1.0 + exp([2]*(x+[3])))")

    barphiphiGraphs = []
    output_phiphiCorr_barrel = []
    for iv in range(0, 8):
        cc.cd(iv + 1)
        x = array('d')
        y = array('d')
        dx = array('d')
        dy = array('d')
        for kk, vv in sorted(fitPars.items()):
            if kk > -fit_emin:
                x.append(log10(kk))
                dx.append(0.)
                y.append(vv[iv])
                dy.append(fitParErrs[kk][iv])
        pg = ROOT.TGraphErrors(len(x), x, y, dx, dy)
        pg.SetName('barrel_phiphi_par'+str(iv))
        pg.SetTitle('barrel_phiphi_par'+str(iv))
        barphiphiGraphs.append(pg)
        pg.Draw("apl")
        pg.GetHistogram().GetXaxis().SetTitle('log10( photon energy [GeV] )')
        pg.GetHistogram().GetYaxis().SetTitle('fitted parameter')

        if iv == 0:
            pg.Fit("fn_phiphiFit", "q")
            for i in range(0, 4):
                output_phiphiCorr_barrel.append(
                    pg.GetFunction("fn_phiphiFit").GetParameter(i))
        elif iv == 3 or iv == 4:
            pg.Fit("pol0", "q")
            output_phiphiCorr_barrel.append(
                pg.GetFunction("pol0").GetParameter(0))
        elif iv == 5 or iv == 7:
            pg.Fit("pol2", "q")
            for i in range(0, 3):
                output_phiphiCorr_barrel.append(
                    pg.GetFunction("pol2").GetParameter(i))
        else:
            pg.Fit("pol1", "q")
            for i in range(0, 2):
                output_phiphiCorr_barrel.append(
                    pg.GetFunction("pol1").GetParameter(i))

    cc.Print(plname)

    return output_phiphiCorr_barrel


# ################################
output_phiphiCorr_barrel = angle_phiBarrel()
# ################################

# =====================================
# theta in the barrel
# =====================================


def angle_thetaBarrel():

    print('theta correction: barrel')

    barthcthHists = {}

    cc.Clear()
    cc.Divide(5, 4)
    ic = 1
    for kk, vv in sorted(inTrees.items()):
        cc.cd(ic)

        hn = "barrel_thcth_"+str(kk)
        if kk < 50:
            yy = 0.01 / sqrt(kk)
        else:
            yy = 0.02 / sqrt(kk)

        cthcth = ROOT.TH2F(hn, hn, 500, -1, 1, 200, -yy, yy)
        cthcth.GetYaxis().SetTitle(thVar+' - trueTheta [rad]')
        if useMC:
            vv.Draw(thVar+"-acos(mcCth):mcCth>>"+hn, eVar+">0")
            cthcth.GetXaxis().SetTitle('Cos(Theta) [MC, rad]')
        else:
            vv.Draw(thVar+"-acos(mcCth):cos(pfoTheta)>>"+hn, eVar+">0")
            cthcth.GetXaxis().SetTitle('Cos(Theta) [PFO, rad]')

        barthcthHists[kk] = cthcth
        cthcth.Draw()
        ic += 1
    cc.Print(plname)

    fn_barthcth = ROOT.TF1("fn_barthcth",
                           "[0]*x + [1]*( 4*pow(x, 3) - 3*x )")

    fitPars = {}
    fitParErrs = {}

    ic = 1
    for kk, cthcth in sorted(barthcthHists.items()):
        cc.cd(ic)
        epx = cthcth.ProfileX()
        epx.Draw()
        epx.GetYaxis().SetTitle(cthcth.GetYaxis().GetTitle())

        fn_barthcth.SetParameter(0, -3.e-4)
        fn_barthcth.SetParameter(1, 1.e-4)

        epx.Fit("fn_barthcth", "q", "", -0.75, 0.75)

        fitPars[kk] = [epx.GetFunction("fn_barthcth").GetParameter(0),
                       epx.GetFunction("fn_barthcth").GetParameter(1)]
        fitParErrs[kk] = [epx.GetFunction("fn_barthcth").GetParError(0),
                          epx.GetFunction("fn_barthcth").GetParError(1)]

        ic += 1
    cc.Print(plname)

    cc.Clear()
    cc.Divide(2, 2)

    barthcthGraphs = []
    output_cthcthCorr_barrel = []
    for iv in range(0, 2):
        cc.cd(iv + 1)
        x = array('d')
        y = array('d')
        dx = array('d')
        dy = array('d')
        for kk, vv in sorted(fitPars.items()):
            if kk > fit_emin:
                x.append(log10(kk))
                dx.append(0.)
                y.append(vv[iv])
                dy.append(fitParErrs[kk][iv])
        pg = ROOT.TGraphErrors(len(x), x, y, dx, dy)
        pg.SetName('barrel_phi_par'+str(iv))
        pg.SetTitle('barrel_phi_par'+str(iv))
        barthcthGraphs.append(pg)
        pg.Draw("apl")
        pg.GetHistogram().GetXaxis().SetTitle('log10( photon energy [GeV] )')
        pg.GetHistogram().GetYaxis().SetTitle('fitted parameter')
        pg.Fit("pol1", "q")

        output_cthcthCorr_barrel.append(pg.GetFunction("pol1").GetParameter(0))
        output_cthcthCorr_barrel.append(pg.GetFunction("pol1").GetParameter(1))

    cc.Print(plname)

    return output_cthcthCorr_barrel


# ################################
output_cthcthCorr_barrel = angle_thetaBarrel()
# ################################

# =====================================
# theta in the endcap
# =====================================


def angle_thetaEndcap():

    print('theta correction: endcap')

    endththHists = {}

    cc.Clear()
    cc.Divide(5, 4)
    ic = 1
    for kk, vv in sorted(inTrees.items()):
        cc.cd(ic)
        hn = "endcap_thth_" + str(kk)
        yy = 0.01/sqrt(kk)
        # x-axis in radians
        thth = ROOT.TH2F(hn, hn, 50, 0.2, 0.6, 200, -yy, yy)
        thth.GetYaxis().SetTitle(thVar + ' - trueTheta')
        if useMC:
            vv.Draw("(mcCth<0 ? -1. : 1.)*( " + thVar +
                    " - acos(mcCth) ):acos(abs(mcCth))>>" + hn, eVar + ">0")
            thth.GetXaxis().SetTitle('theta [MC]')
        else:
            vv.Draw("(cos(pfoTheta)<0 ? -1. : 1.)*( " + thVar +
                    "-acos(mcCth)):acos(abs(cos(pfoTheta)))>>" + hn, eVar+">0")
            thth.GetXaxis().SetTitle('theta [PFO]')

        endththHists[kk] = thth
        thth.Draw()
        ic += 1
    cc.Print(plname)

    fn_endthth = ROOT.TF1("fn_endthth", "[0] + [1]*x + [2]*x*x")
    fitPars = {}
    fitParErrs = {}

    cc.Clear()
    cc.Divide(5, 4)
    ic = 1
    for kk, vv in sorted(endththHists.items()):
        cc.cd(ic)
        tp = vv.ProfileX()
        tp.GetYaxis().SetTitle(thth.GetYaxis().GetTitle())

        fn_endthth.SetParameter(0, 2.e-4)
        fn_endthth.SetParameter(1, -4.e-4)
        fn_endthth.SetParameter(2, 4.e-4)

        tp.Fit("fn_endthth", "q")
        tp.Draw()

        fitPars[kk] = [tp.GetFunction("fn_endthth").GetParameter(0),
                       tp.GetFunction("fn_endthth").GetParameter(1),
                       tp.GetFunction("fn_endthth").GetParameter(2)]
        fitParErrs[kk] = [tp.GetFunction("fn_endthth").GetParError(0),
                          tp.GetFunction("fn_endthth").GetParError(1),
                          tp.GetFunction("fn_endthth").GetParError(2)]

        ic += 1
    cc.Print(plname)

    cc.Clear()
    cc.Divide(2, 2)

    endththGraphs = []
    output_ththCorr_endcap = []
    for iv in range(0, 3):
        cc.cd(iv + 1)
        x = array('d')
        y = array('d')
        dx = array('d')
        dy = array('d')
        for kk, vv in sorted(fitPars.items()):
            if kk > fit_emin:
                x.append(log10(kk))
                dx.append(0.)
                y.append(vv[iv])
                dy.append(fitParErrs[kk][iv])
        pg = ROOT.TGraphErrors(len(x), x, y, dx, dy)
        pg.SetName('endcap_theta_par' + str(iv))
        pg.SetTitle('endcap_theta_par' + str(iv))
        endththGraphs.append(pg)
        pg.Draw("apl")
        pg.GetHistogram().GetXaxis().SetTitle('log10( photon energy [GeV] )')
        pg.GetHistogram().GetYaxis().SetTitle('fitted parameter')
        pg.Fit("pol0", "q")

        output_ththCorr_endcap.append(pg.GetFunction("pol0").GetParameter(0))
        # energy dependence negligible, change fn to pol0
        output_ththCorr_endcap.append(0.)  # energy dep

    cc.Print(plname)

    return output_ththCorr_endcap


# ################################
output_ththCorr_endcap = angle_thetaEndcap()
# ################################


# #################################################
# finally we write out the text which can go into the steering xml
# #################################################
print('====================================================================')
print('The following to go into steering xml for photonCorrectionProcessor:')
print('====================================================================')
print('<parameter name="energyCor_Linearise"> ', end=" ")
for i in range(0, len(energyCor_Linearise)):
    print(f"{energyCor_Linearise[i]:.4g} ", end=" ")
print('</parameter>')

print('<parameter name="energyCorr_barrelPhi"> ', end=" ")
for i in range(0, len(output_energyCorr_barrelPhi)):
    print(f"{output_energyCorr_barrelPhi[i]:.4g} ", end=" ")
print('</parameter>')

print('<parameter name="energyCorr_costh"> ', end=" ")
for i in range(0, len(output_energyCorr_costh)):
    print(f"{output_energyCorr_costh[i]:.4g} ", end=" ")
print('</parameter>')

print('<parameter name="energyCorr_endcap"> ', end=" ")
for i in range(0, len(output_energyCorr_endcap)):
    print(f"{output_energyCorr_endcap[i]:.4g} ", end=" ")
print('</parameter>')

print('<parameter name="thetaCorr_barrel"> ', end=" ")
for i in range(0, len(output_cthcthCorr_barrel)):
    print(f"{output_cthcthCorr_barrel[i]:.4g} ", end=" ")
print('</parameter>')

print('<parameter name="thetaCorr_endcap"> ', end=" ")
for i in range(0, len(output_ththCorr_endcap)):
    print(f"{output_ththCorr_endcap[i]:.4g} ", end=" ")
print('</parameter>')

print('<parameter name="phiCorr_barrel"> ', end=" ")
for i in range(0, len(output_phiphiCorr_barrel)):
    print(f"{output_phiphiCorr_barrel[i]:.4g} ", end=" ")
print('</parameter>')

print('====================================================================')
print('========== END =====================================================')
print('====================================================================')

cc.Print(plname+']')
fin.Close()
