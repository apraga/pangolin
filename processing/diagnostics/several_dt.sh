#!/bin/sh
# Plot several diagnostics

ffile="cfl-impact-un.dat"
echo "#dt l2 loo" > $ffile

config=config
for i in 9 8 6 5 4 1; do 
  file="cfl-impact-"$i"lat-undat"
  echo '&files' > $config
  echo 'folder = "/wkdir/pae2/praga/parallel/cfl_impact/gaussian_160lat/",' >> $config
  echo 'init_file = "ratio_dt_'$i'_0.dat",' >> $config
  echo 'final_file = "ratio_dt_'$i'_T.dat",' >> $config
  echo 'out_file = "'$file'",' >> $config
  echo 'nb_lat2 = 160/' >> $config

  ./pango_diagnostics --err
  # Replace number of cells (first column) by time step
  tmp=`awk -v tmp=$i '{print tmp, $3, $4, $5 }' $file`

  # Need the quote for new line
  # Skip header 
  echo "$tmp" | tail -n+2>> $ffile
done
