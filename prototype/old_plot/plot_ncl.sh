#!/bin/sh
if [ $# -lt 2 ]
then
  echo "Usage : plot_ncl.sh SCRIPT [type] NUMBER_DATA"
  echo "script = contour, raster, point, polygon_oo"
  echo "type = ortho (default mercator)"
  exit 
fi

input="visu_"$1".ncl"
# Mercator by default
if [ $# -gt 2 ]
then
  file="visu_"$1"_"$3".ps"
  type="type=\""$2"\""
  data="nb="$3
  rot=""
else
  file="visu_"$1"_"$2".ps"
  type="type=\"""\""
  data="nb="$2
  rot="-rotate -90"
fi

echo "Starting NCL..."
ncl $input $type $data 
echo "NCL finished. Now converting to jpg..."

jpg=`basename "$file" .ps`
jpg=$jpg".jpg"
convert $rot $file $jpg
eog $jpg

echo "done."
