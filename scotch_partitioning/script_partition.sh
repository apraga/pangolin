#!/bin/bash

dir="/space/praga/softwares/scotch_5.1.12/bin"
name="mesh_init"
graph=$name".grf"
map=$name".map"
geom=$name".xyz"
tmp="k2.tgt"
options="-Op{e,c,l}"
out="tree.ps"
out2="tree.png"
shopt -s expand_aliases
alias gmap2=$dir"/gmap -vm"
alias gout2=$dir"/gout"
#
nb_parts=72
echo cmplt $nb_parts > $tmp
time gmap2 mesh_init.grf $tmp $map
gout2 $graph $geom  $map $options > $out
ps2raster $out -Tg
eog $out2
