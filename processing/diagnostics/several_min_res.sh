#!/bin/sh
# Plot several diagnostics

ffile="min-res-pangolin-un-CFL1.0.dat"
echo "#resolution nb_cells l2 loo" > $ffile

config=config
#for i in 20 40 80 160 320 640; do 
for i in 160 320 640 800 1000; do 
  echo '&files' > $config
  echo 'folder = "/wkdir/pae2/praga/parallel/min_res/",' >> $config
  echo 'init_file = "ratio_'$i'lat_un_CFL1.0_0.dat",' >> $config
  echo 'final_file = "ratio_'$i'lat_un_CFL1.0_T.dat",' >> $config
  echo 'out_file = "min-res-pangolin-'$i'lat-CFL1.0.dat",' >> $config
  echo 'nb_lat2 = '$i'/' >> $config

  ./pango_diagnostics --err
  file="min-res-pangolin-"$i"lat-CFL1.0.dat"
  # Skip header 
  cat $file | tail -n+2>> $ffile
done
