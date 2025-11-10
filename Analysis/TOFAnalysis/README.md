# TOFAnalysis

This project contains custom codes used for Bohdan Dudar's PhD studies of time-of-flight particle identification.
It contains the BohdanAna Marlin processor that dumps a lot of information to the ROOT file for further ploting with python(pyROOT) scripts in analysis folder.
It is a copy of the [BohdanAna](https://github.com/dudarboh/BohdanAna/releases/tag/phd_end_add_to_marlinreco) private repository but integrated inside MarlinReco package.

## How to run

1. Setup the key4hep/iLCSoft environment with the MarlinReco version containing this project. The easiest is to use one of the cvmfs builds:

```
# use one of these. Not all at the same time...
source /cvmfs/ilc.desy.de/key4hep/setup.sh
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh
```

2. You can check if your MarlinReco contains BohdanAna with the line below. If you do not see any output, you are likely need a newer version of the MarlinReco

```
nm -o $(echo $MARLIN_DLL | tr -t ':' '\n' | grep MarlinReco) | grep BohdanAna
```

3. Copy the steering file so you have a copy of `./xml/steer.xml` and change the path to the geometry file where your geometry is located. You can check with `echo $k4geo_DER``

4. Run ` Marlin /path/to/your/steer.xml `

Ideally it should run and produce the correposponding ROOT file.

## Running on the HTCondor

To run many instances in paralel on HTCondor, write all your input slcio files file paths into some txt file, which is used in `send_jobs.sub` and you should modify:

The txt file may contain something like:

```
/pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/2f_hadronic_eL_pR/ILD_l5_o1_v02/v02-02/00015161/000/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n000_001.d_rec_00015161_328.slcio
/pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/2f_hadronic_eL_pR/ILD_l5_o1_v02/v02-02/00015161/000/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n000_002.d_rec_00015161_221.slcio
/pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/2f_hadronic_eL_pR/ILD_l5_o1_v02/v02-02/00015161/000/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n000_003.d_rec_00015161_24.slcio
/pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/2f_hadronic_eL_pR/ILD_l5_o1_v02/v02-02/00015161/000/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n000_004.d_rec_00015161_121.slcio
/pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/2f_hadronic_eL_pR/ILD_l5_o1_v02/v02-02/00015161/000/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n000_005.d_rec_00015161_104.slcio
/pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/2f_hadronic_eL_pR/ILD_l5_o1_v02/v02-02/00015161/000/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n000_006.d_rec_00015161_69.slcio
/pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/2f_hadronic_eL_pR/ILD_l5_o1_v02/v02-02/00015161/000/rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500010.P2f_z_h.eL.pR.n000_007.d_rec_00015161_44.slcio
```

Then each file will be asined to a separate job on the cluster.
NOTE: the files must be accesible from the HTCondor nodes via shared file system.

5. To run several instances on HTCondor create your copy of `./job/job.sh`, `./job/send_jobs.sub`, and a  run `condor_submit ./job/send_jobs.sub`. You need to modify all the paths according to you.


## Content (Short)

```
.
├── src/BohdanAna.cc                Main Marlin processor.
├── src/CreateRefitPFO.cc           Marlin processor used for reffited track studies.
├── xml/steer.xml                   Marlin steering file to run.
├── src/TOF.cc                      Utility functions for various TOF calculations.
├── src/TrackLength.cc              Utility functions for various track length calculations.
├── src/BohdanUtils.cc              Utility functions (generic).
├── src/BohdanDrawing.cc            Debugging functions for drawing during Marlin execution.
├── job/job.sh                      Bash script for running jobs on HTCondor.
├── analysis                        Python scripts that make plots.
```

## Content (Detailed)

```
.
├── src/BohdanAna.cc                Main Marlin processor.
                                    It calculates a lot of info from the input slcio REC file and stores it in the ROOT file per particle.
                                    steering parameters:
                                        - eventDisplay: false (default) / true; enables event display debugging during the executing with drawDisplay()
                                        - produce_csv_output: false (default) / true; dumps some of the information into the csv file rather than ROOT file for Konrad ML studies
                                        - produce_refit_output: false (default) / true; produces a copy of output for the particle with the refitted track with the true mass hypothesis.
                                            NOTE: Rerunning vertexing with refitted PFOs slows things down significantly.
                                        - dst_mode: false (default) / true; outputs only information available from DST files, i.e. no TOF and hit info. Makes some plots possible with larger statistics of DST files.

├── src/CreateRefitPFO.cc           Marlin processor used for reffited track studies.
                                    This processor creates a new collection (updatedPandoraPFOs) of reconstructed particles from existing PandoraPFOs collection.
                                    Each PFO of the corresponding kaon/proton MCParticle is recreated with a track refitted with the true (kaon/proton) mass hypothesis adapting the PFO parameters as well.

├── xml/steer.xml                   Marlin steering file to run. Explanation of the execute section:
                                    InitDD4hep - is required to initialise the geometry. Only ILD_l5_o1_v02 is expected to work as many things are very geometry dependant.
                                    Other geometries may need serious code/algorithm adaptations
                                    MergeHitRelations - is required for some parts of the code. It makes finding a corresponding MC-Reco link easier with a single loop without the need to loop over all subdetector collections separately.
                                    Correct_dEdx - is not required, but improves the dE/dx measurement by correcting angular dependence of the 2020 ILD MC production.
                                    CreateRefitPFO - is required only if one wants outputs with refitted tracks.
                                    VertexFinderRefit - is required only if one wants outputs with vertices obtained from refitted PFO objects. All parameters (except collection names) are copied from the default ILD reconstruction.
                                    BohdanAna - is the central part of this chain. It calculates a lot of variables and stores them in the ROOT file for forther plotting/analysis with python scripts.

├── src/TOF.cc                      Utility functions for various TOF calculations.
                                    Some of them may be obsolete, but they give a general idea how to do TOF reconstruction/hit selection.
                                    As individual hit information from REC files is also written in the ROOT file, all these methods exist in analysis/tof.hpp and are tested with python scripts including more variety.
                                    E.g. TOF cut parameter scan and different time resolutions are done with analysis/tof.hpp. Thus, this file is not used so much.

├── src/TrackLength.cc              Utility functions for various track length calculations.
                                    It has all different methods for calculating the track lengths.

├── src/BohdanUtils.cc              Utility functions (generic).
                                    Some of the functions may already exist in MarlniUtil or LCIO/UTIL...

├── src/BohdanDrawing.cc            Debugging functions for drawing during Marlin execution.
                                    These functions are used purely for debugging purposes and pause the event loop after every event/particle depending on what is put in the main processor BohdanAna.
                                    These functions allow to do ROOT plots and draw event display inside Marlin event loop and plot some useful information of a particular event or particle.

├── job/job.sh                      Bash script for running jobs on HTCondor. Simply setups the environment and runs the Marlin xml/steer.xml.
                                    It is used only in job/send_jobs.sub file to submit to the HTCondor. It allows to obtaine way better statistics.
                                    After all HTCondor nodes produce their ROOT files, they cab be merged using ROOT command line tool - hadd.

├── analysis                        Python scripts that make plots.
                                    The scripts in this folder were oftenly changing, as one typically needs to make a small adjustments in the figures.
                                    Thus, these are the quick and dirty plotting scripts and do not expect high quality code here.
                                    Nevertheless, it should give a good understanding of the workflow of plotting the content of the ROOT file and analysing data with RDataFrame.
                                    Some python scripts produce purely analytical results (w/o any data).
                                    The best idea is to use these scripts as a refence when creating your own scripts that produce plots.
```










