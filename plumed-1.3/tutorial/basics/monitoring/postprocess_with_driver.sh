#!/bin/bash
CATDCD="catdcd"
TRJOCONV="trjconv_d"
DRIVER="driver_plumed_devel" 
mv COLVAR COLVAR.bak
#
# generate the pdb
#
echo 0 | $TRJCONV -f 2ala.gro -o 2ala.pdb
#
# generate the dcd 
#
$CATDCD  -o traj.dcd -trr traj.trr 
#
# run the driver
#
$DRIVER -pdb 2ala.pdb -dcd traj.dcd -ncv 2 -nopbc -plumed plumed2.dat 
