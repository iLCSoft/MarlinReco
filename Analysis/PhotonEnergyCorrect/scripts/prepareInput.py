from pyLCIO import IOIMPL
import ROOT
import math
from array import array
import glob
#
# this prepares root trees used to derive the photon energy corrections
#
# Daniel Jeans, 2022-2025
#

# ###########################
# define input single photons samples to use
#  edit this to find your samples
#  filenames should be put into the filenames dictionary
# ###########################
filenames = {}

enlabs = ['0p2', '0p4', '0p7', '001', '002', '3p5', '005', '7p5', '010',
          '015', '030', '050', '100', '125', '175', '250', '350', '500']

# the orig files
# indir='/hsm/ilc/grid/storm/prod/ilc/mc-opt/ild/dst-merged/1-calib/single/ILD_l5_o1_v02_nobg/v02-01'
# label='rv02-01.sv02-01'

# and my samples
indir = 'runReco/output'
# label = 'noPhotonCorr'
label = 'new3PhotonCorr'

for enlab in enlabs:
    inpatt = indir + '/*' + label + '*.Ps_22_' + enlab + '*DST.slcio'
    infiles = glob.glob(inpatt)
    infiles.sort()
    filenames[enlab] = infiles

maxevt = 0    # evts per energy point. <=0 : all events

# input collections
inputPFOcollection = 'PandoraPFOs'
inputMCcollection = 'MCParticle'

# output file name
outname = "singlegammas_newSample_"+label+".root"

ROOT.gROOT.SetBatch()

# ########################################################
#
# some geometrical info for the model being considered
#  these are for the ILD_l5_v02 model (or any of the large ILD models)
#
# <constant name="TPC_Ecal_Hcal_barrel_halfZ" value="2350*mm
# <constant name="Ecal_Tpc_gap" value="35*mm"/>
# <constant name="Ecal_endcap_zmin"
#           value="TPC_Ecal_Hcal_barrel_halfZ + 61.8*mm"/>
# <constant name="Ecal_endcap_center_box_size" value="800.0*mm"/>
# <constant name="Ecal_Tpc_gap" value="35*mm"/>
# <constant name="top_TPC_outer_radius"    value="1769.8*mm"/>

# inner radius of barrel ECAL
_barrel_rin = 1769.8 + 35.     # mm
# z-pos inner face of endcap ECAL
_assumed_endZ = 2350. + 61.8   # mm
# costh of barrel-endcap interface
_barrelendcap_costhlimit = math.cos(math.atan(_barrel_rin/_assumed_endZ))
# the central box in the endcap ECAL (half-width)
_assumed_boxsize = 800./2      # mm

print('_barrel_rin', _barrel_rin)
print('_assumed_endZ', _assumed_endZ)
print('_barrelendcap_costhlimit', _barrelendcap_costhlimit)
print('_assumed_boxsize', _assumed_boxsize)

# ########################################################
#
# some geometrical transformations
#  suitable for octagonal barrel, 4-fold endcap as used in ILD
#
# ########################################################


def getBarrelFoldedPhi(phi):
    # fold the barrel phi into a single octant
    if phi < 0:
        phi += 2*math.pi
    sector = int(phi/(math.pi/4))
    foldedPhi = phi - sector*math.pi/4.
    return foldedPhi

# #######################################################


def getDistanceAcrossEndcapQuadrant(costh, phi):
    # this calculates the distance across an endcap quadrant
    #  (ie perpendicular to slab direction),
    # from the inner edge of quadrant
    if abs(costh) < _barrelendcap_costhlimit:
        # print("ERROR getDistanceAcrossEndcapQuadrant : not in endcap!")
        return -999
    assumed_endZ = _assumed_endZ
    if costh < 0:
        assumed_endZ *= -1

    endX = _assumed_endZ * math.sin(math.acos(costh)) * math.cos(phi)
    endY = _assumed_endZ * math.sin(math.acos(costh)) * math.sin(phi)
    quad = -1
    if costh > 0:
        if endX > -_assumed_boxsize and endY > _assumed_boxsize:
            quad = 0
        elif endX > _assumed_boxsize and endY < _assumed_boxsize:
            quad = 1
        elif endX < _assumed_boxsize and endY < -_assumed_boxsize:
            quad = 2
        elif endX < -_assumed_boxsize and endY > -_assumed_boxsize:
            quad = 3
    else:
        if endX < _assumed_boxsize and endY > _assumed_boxsize:
            quad = 0
        elif endX > _assumed_boxsize and endY > -_assumed_boxsize:
            quad = 1
        elif endX > -_assumed_boxsize and endY < -_assumed_boxsize:
            quad = 2
        elif endX < -_assumed_boxsize and endY < _assumed_boxsize:
            quad = 3

    foldY = 0.
    if quad >= 0:    # not in the center box
        foldPhi = phi + quad*math.pi/2.
        foldY = _assumed_endZ*math.sin(math.acos(costh))*math.sin(foldPhi)
    return foldY


# ########################################################

rfile = ROOT.TFile(outname, "recreate")

# define trees and branches
outtree = {}

mcPhi = array('f', [0])
mcCth = array('f', [0])
mcPhiFold = array('f', [0])
mcEndX = array('f', [0])

pfoEn = array('f', [0])
pfoTheta = array('f', [0])
pfoPhi = array('f', [0])
pfoPhiFold = array('f', [0])
pfoEndX = array('f', [0])

othEn = array('f', [0])

clEn = array('f', [0])
clTheta = array('f', [0])
clPhi = array('f', [0])
clPhiFold = array('f', [0])
clEndX = array('f', [0])

for enlab in filenames.keys():

    outtree[enlab] = ROOT.TTree('gamTree_'+enlab, 'gamTree_'+enlab)
    outtree[enlab].Branch('pfoEn', pfoEn, 'pfoEn/F')

    outtree[enlab].Branch('pfoTheta', pfoTheta, 'pfoTheta/F')
    outtree[enlab].Branch('pfoPhi', pfoPhi, 'pfoPhi/F')
    outtree[enlab].Branch('pfoPhiFold', pfoPhiFold, 'pfoPhiFold/F')
    outtree[enlab].Branch('pfoEndX', pfoEndX, 'pfoEndX/F')

    outtree[enlab].Branch('clEn', clEn, 'clEn/F')
    outtree[enlab].Branch('clTheta', clTheta, 'clTheta/F')
    outtree[enlab].Branch('clPhi', clPhi, 'clPhi/F')
    outtree[enlab].Branch('clPhiFold', clPhiFold, 'clPhiFold/F')
    outtree[enlab].Branch('clEndX', clEndX, 'clEndX/F')

    outtree[enlab].Branch('othEn', othEn, 'othEn/F')
    outtree[enlab].Branch('mcPhi', mcPhi, 'mcPhi/F')
    outtree[enlab].Branch('mcCth', mcCth, 'mcCth/F')
    outtree[enlab].Branch('mcPhiFold', mcPhiFold, 'mcPhiFold/F')
    outtree[enlab].Branch('mcEndX', mcEndX, 'mcEndX/F')


# #################################################################
# fill tree for a given energy point
#


def run(en):

    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    fnames = filenames[en]
    nevt = 0

    for filename in fnames:

        if maxevt > 0 and nevt > maxevt:
            break

        print(en, indir, filename)

        reader.open(filename)

        for event in reader:

            if maxevt > 0 and nevt > maxevt:
                break
            nevt += 1

            mcgamconv = False
            mcgamcosth = -99.
            mcgamphi = -99.

            try:
                mcps = event.getCollection(inputMCcollection)
                for mcp in mcps:
                    if mcp.getPDG() == 22 and mcp.getGeneratorStatus() == 1:
                        mcgamconv = mcp.isDecayedInTracker()
                        mcgamcosth = mcp.getMomentum()[2]/mcp.getEnergy()
                        mcgamphi = math.atan2(mcp.getMomentum()[1],
                                              mcp.getMomentum()[0])
                        break
            except Exception:
                print('some problem dealing with',
                      inputMCcollection, 'collection')

            if mcgamconv:
                # this photon converted in the tracker:
                #  don't use it for calo calibration
                continue

            totpfogamen = 0.0
            totpfoothen = 0.0
            totclgamen = 0.0
            ngam = 0
            try:
                pfos = event.getCollection(inputPFOcollection)
                totClMom = ROOT.TVector3()
                totPfoMom = ROOT.TVector3()
                for p in pfos:
                    if p.getType() == 22:
                        totpfogamen += p.getEnergy()
                        ngam = ngam + 1
                        totPfoMom += ROOT.TVector3(p.getMomentum())
                        for cl in p.getClusters():
                            totclgamen += cl.getEnergy()
                            pos = ROOT.TVector3(cl.getPosition())
                            totClMom += pos*(cl.getEnergy()/pos.Mag())
                    else:
                        totpfoothen += p.getEnergy()
            except Exception:
                print('some problem dealing with ',
                      inputPFOcollection, ' collection')

            mcCth[0] = mcgamcosth
            mcPhi[0] = mcgamphi
            mcPhiFold[0] = getBarrelFoldedPhi(mcgamphi)
            mcEndX[0] = getDistanceAcrossEndcapQuadrant(mcgamcosth, mcgamphi)

            pfoEn[0] = totpfogamen
            pfoTheta[0] = totPfoMom.Theta()
            pfoPhi[0] = totPfoMom.Phi()
            pfoPhiFold[0] = getBarrelFoldedPhi(totPfoMom.Phi())
            pfoEndX[0] = getDistanceAcrossEndcapQuadrant(totPfoMom.CosTheta(),
                                                         totPfoMom.Phi())

            clEn[0] = totclgamen
            clTheta[0] = totClMom.Theta()
            clPhi[0] = totClMom.Phi()
            clPhiFold[0] = getBarrelFoldedPhi(totClMom.Phi())
            clEndX[0] = getDistanceAcrossEndcapQuadrant(totClMom.CosTheta(),
                                                        totClMom.Phi())

            othEn[0] = totpfoothen
            outtree[en].Fill()

        reader.close()

    print('completed energy', en, 'nevt=', nevt)


# #########################################################
# loop over energy points
for enstr in filenames.keys():
    run(enstr)

# close output file
rfile.Write()
rfile.Close()
# #########################################################
# that's it!
