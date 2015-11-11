cd w1/
/Users/carlo/Codes/gromacs-4.5.5/bin/mdrun_mpi -nice 0 -s topol.tpr -g /Users/carlo/Codes/md_meta/tests/gromacs/TEST-gromacs/multiple/topol.tpr0 -plumed &>/dev/null &
cd ../w2
/Users/carlo/Codes/gromacs-4.5.5/bin/mdrun_mpi -nice 0 -s topol.tpr -g /Users/carlo/Codes/md_meta/tests/gromacs/TEST-gromacs/multiple/topol.tpr1 -plumed &>/dev/null
cd ../
