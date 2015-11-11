#!/bin/bash

EXE=namd_plumed

cd w1/
${EXE} namd_input > namd.out.1.save &
cd ../w2
${EXE} namd_input > namd.out.2.save 
cd ../
