#!/bin/sh
# Plot several diagnostics

ffile="conv-pangolin-alternate-sp-CFL0.96.dat"
echo "#resolution nb_cells l2 loo" > $ffile

config=config
#for i in 20 40 80 160 320 640; do 
for i in 160 320 640; do 
  file="tmp-"$i"lat"
  echo '&files' > $config
  echo 'folder = "/wkdir/pae2/praga/parallel/cv_rate_alternate/",' >> $config
  echo 'init_file = "'$i'lat/ratio_1_201301010000.h5",' >> $config
  echo 'final_file = "'$i'lat/ratio_1_201301130000.h5",' >> $config
  echo 'out_file = "'$file'"' >> $config
  echo 'nb_lat2 = '$i'/' >> $config

  ./pango_diagnostics --err
  # Skip header 
  cat $file | tail -n+2>> $ffile
done
