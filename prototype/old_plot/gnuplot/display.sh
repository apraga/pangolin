#!/bin/bash

if [ "$#" -lt 1 ]
then
  echo "1 argument required $# provided"
  exit
fi

rm temp.gnu
echo "set pm3d" >> temp.gnu
echo "set view map" >> temp.gnu
echo "set nokey" >> temp.gnu

if [ -n "$2" ]
then
  echo "set terminal postscript color" >> temp.gnu
  echo "set output '"$2"'" >> temp.gnu
fi

echo "splot 'tmp/"$1".dat' u 1:2:3:4 palette pt 5 ps 0.6" >> temp.gnu 
#echo "splot 'tmp/"$1".dat' u 1:2:3:4 palette with lines" >> temp.gnu 

/usr/bin/gnuplot -p temp.gnu
