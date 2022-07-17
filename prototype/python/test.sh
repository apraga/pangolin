#!/bin/sh
for i in `seq 1150 1200`
do
  echo $i" : " && python error_norm.py 2 "0.dat" $i".dat"
done
