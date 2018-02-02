 
#
#  create root plots
#

# root -b -q 'draw_cluster_time.C("beta_05hits_50ps")'
# root -b -q 'draw_cluster_time.C("beta_05perc_50ps")'
# root -b -q 'draw_cluster_time.C("beta_10hits_50ps")'
# root -b -q 'draw_cluster_time.C("beta_10perc_50ps")'
# root -b -q 'draw_cluster_time.C("beta_20hits_50ps")'
# root -b -q 'draw_cluster_time.C("beta_20perc_50ps")'

root -b -q 'draw_cluster_time.C("clutime_0ps")'
root -b -q 'draw_cluster_time.C("time05perc_0ps")'
root -b -q 'draw_cluster_time.C("time10perc_0ps")'
root -b -q 'draw_cluster_time.C("time20perc_0ps")'
root -b -q 'draw_cluster_time.C("timeoffastesthit_0ps")'
root -b -q 'draw_cluster_time.C("time05hits_0ps")'
root -b -q 'draw_cluster_time.C("time10hits_0ps")'
root -b -q 'draw_cluster_time.C("time20hits_0ps")'
root -b -q 'draw_cluster_time.C("cor_ref_time05perc_0ps")'
root -b -q 'draw_cluster_time.C("cor_ref_time10perc_0ps")'
root -b -q 'draw_cluster_time.C("cor_ref_time20perc_0ps")'
root -b -q 'draw_cluster_time.C("cor_ref_time05hits_0ps")'
root -b -q 'draw_cluster_time.C("cor_ref_time10hits_0ps")'
root -b -q 'draw_cluster_time.C("cor_ref_time20hits_0ps")'
