#!/bin/sh

python ../python/partitioning.py
mpost partition.mp && gv partition-1.mps &
