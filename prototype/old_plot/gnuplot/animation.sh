#!/bin/bash

if [ ! "$#" -eq 1 ]
then
  echo "1 argument required $# provided"
  exit
fi

echo "set pm3d" >> temp.gnu
echo "set view map" >> temp.gnu
echo "set terminal jpeg" >> temp.gnu
for i in {1..2499}
do
  echo "set output 'pictures/pic"$i".jpg'" >> temp.gnu
  echo "splot 'tmp/"$i".data' u 1:2:3:4 palette pt 5 ps 0.5" >> temp.gnu 
done

/usr/bin/gnuplot temp.gnu
rm temp.gnu
ffmpeg -f image2 -r 25 -i pictures/pic%d.jpg -b 600k $1

