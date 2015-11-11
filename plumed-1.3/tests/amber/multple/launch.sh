#!/bin/bash

EXE=/Users/chicco/bin/sander_plugin

cd w1/
${EXE} -O -i amber_input -o global_info -p ../../DIALA_TOP/system.prmtop -c ../../DIALA_TOP/system.inpcrd -inf out.1.save &
cd ../w2
${EXE} -O -i amber_input -o global_info -p ../../DIALA_TOP/system.prmtop -c ../../DIALA_TOP/system.inpcrd -inf out.2.save &
cd ../
