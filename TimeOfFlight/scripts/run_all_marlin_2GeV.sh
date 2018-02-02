
#
# create root files with time of flight plots for all particles with 2 GeV
#

 Marlin mysteer.xml --global.LCIOInputFiles=/pnfs/desy.de/ilc/prod/ilc/mc-opt/ild/rec/calib/single_PDG13/ILD_l4_v02/v01-19-04_lcgeo/u014/rv01-19-04_lcgeo.sv01-19-04_lcgeo.mILD_l4_v02.Pmcparticles_PDG13_MOM002GeV.n001_1.d_rec_u014.slcio --MyAIDAProcessor.FileName=aida_file_muon_02GeV_REC

 Marlin mysteer.xml --global.LCIOInputFiles=/pnfs/desy.de/ilc/prod/ilc/mc-opt/ild/rec/calib/single_PDG321/ILD_l4_v02/v01-19-04_lcgeo/u017/rv01-19-04_lcgeo.sv01-19-04_lcgeo.mILD_l4_v02.Pmcparticles_PDG321_MOM002GeV.n001_1.d_rec_u017.slcio --MyAIDAProcessor.FileName=aida_file_kaon_02GeV_REC

 Marlin mysteer.xml --global.LCIOInputFiles=/pnfs/desy.de/ilc/prod/ilc/mc-opt/ild/rec/calib/single_PDG11/ILD_l4_v02/v01-19-04_lcgeo/u014/rv01-19-04_lcgeo.sv01-19-04_lcgeo.mILD_l4_v02.Pmcparticles_PDG11_MOM002GeV.n001_1.d_rec_u014.slcio --MyAIDAProcessor.FileName=aida_file_electron_02GeV_REC

 Marlin mysteer.xml --global.LCIOInputFiles=/pnfs/desy.de/ilc/prod/ilc/mc-opt/ild/rec/calib/single_PDG211/ILD_l4_v02/v01-19-04_lcgeo/u014/rv01-19-04_lcgeo.sv01-19-04_lcgeo.mILD_l4_v02.Pmcparticles_PDG211_MOM002GeV.n001_1.d_rec_u014.slcio --MyAIDAProcessor.FileName=aida_file_pion_02GeV_REC

 Marlin mysteer.xml --global.LCIOInputFiles=/pnfs/desy.de/ilc/prod/ilc/mc-opt/ild/rec/calib/single_PDG2212/ILD_l4_v02/v01-19-04_lcgeo/u017/rv01-19-04_lcgeo.sv01-19-04_lcgeo.mILD_l4_v02.Pmcparticles_PDG2212_MOM002GeV.n001_1.d_rec_u017.slcio --MyAIDAProcessor.FileName=aida_file_proton_02GeV_REC


root -b -q 'draw_cluster_time.C("beta_05hits_02GeV_0ps")'
root -b -q 'draw_cluster_time.C("beta_05perc_02GeV_0ps")'
root -b -q 'draw_cluster_time.C("beta_10hits_02GeV_0ps")'
root -b -q 'draw_cluster_time.C("beta_10perc_02GeV_0ps")'
root -b -q 'draw_cluster_time.C("beta_20hits_02GeV_0ps")'
root -b -q 'draw_cluster_time.C("beta_20perc_02GeV_0ps")'
