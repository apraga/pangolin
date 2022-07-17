#!/bin/sh
# Plot several diagnostics

# Re-generate the config file each time
config=config_extract
#folder="/wkdir/pae2/praga/parallel/cv_rate/"
folder="/wkdir/pae2/praga/parallel/min_res/"
#for i in 20 40 80 160 320 640; do 

for i in 160 320 640 800; do 
  echo '&nfile' > $config
  echo 'input = "'$folder'ratio_'$i'lat_un_CFL1.0_T.dat",' >> $config
  echo 'output = "extracted_'$i'_un_T.dat"/' >> $config
  ./extract_equator 
done

# Reference solution
i=160
echo '&nfile' > $config
echo 'input = "'$folder'ratio_'$i'lat_un_CFL1.0_0.dat",' >> $config
echo 'output = "extracted_'$i'_un_0.dat"/' >> $config
./extract_equator 

