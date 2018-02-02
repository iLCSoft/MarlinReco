
#
# create root files with time of flight plots for all particles
#

Marlin mysteer.xml --global.LCIOInputFiles=../../data/mcparticles_proton_5GeV_60deg_REC.slcio --MyAIDAProcessor.FileName=aida_file_proton_5GeV_60deg_REC
Marlin mysteer.xml --global.LCIOInputFiles=../../data/mcparticles_kaon_5GeV_60deg_REC.slcio --MyAIDAProcessor.FileName=aida_file_kaon_5GeV_60deg_REC
Marlin mysteer.xml --global.LCIOInputFiles=../../data/mcparticles_muon_5GeV_60deg_REC.slcio --MyAIDAProcessor.FileName=aida_file_muon_5GeV_60deg_REC
Marlin mysteer.xml --global.LCIOInputFiles=../../data/mcparticles_pion_5GeV_60deg_REC.slcio --MyAIDAProcessor.FileName=aida_file_pion_5GeV_60deg_REC
Marlin mysteer.xml --global.LCIOInputFiles=../../data/mcparticles_electron_5GeV_60deg_REC.slcio --MyAIDAProcessor.FileName=aida_file_electron_5GeV_60deg_REC

