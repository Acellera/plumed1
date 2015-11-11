#!/bin/bash
#
# do the grompp 
#
grompp_d -f md.mdp -c 2ala.gro -p gromacs.top  
mdrun_d -plumed plumed 
