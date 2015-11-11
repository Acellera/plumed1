#!/bin/bash
nstep=25
incr=0
ks=100
exe="/Users/chicco/Programs/gronamd/md_meta/codes/namd_download/NAMD_2.6_Source/MacOSX-i686-g++/namd2"
rm -rf COLVAR.* energy CV_* metadatafile
cp init.pdb start.pdb 
rm -rf endconf_*
pi=3.14159265358979
for i in  `seq 0 $nstep`
do
val1=`echo $i\*2\* $pi/$nstep - $pi   | bc -l `
val1b=`echo $i\*2\* $pi/$nstep    | bc -l `
for j in  `seq 0 $nstep`
do
# 0 go forth 1 go back
dir=`echo ${i}%2 | bc `
if  [ "$dir" -eq "0" ]; then
val2=`echo $j\*2\* $pi /$nstep - $pi | bc -l `
val2b=`echo $j\*2\* $pi /$nstep  | bc -l `
else 
val2=`echo $pi -$j\*2\* $pi /$nstep  | bc -l `
val2b=`echo $pi -$j\*2\* $pi /$nstep + $pi | bc -l `
fi
incr=$((incr + 1)) 
echo "ANG1 $val1 ANG2 $val2 INCR $incr"
cat >metadyn.cfg <<EOF
PRINT W_STRIDE 10 
TORSION LIST 5 7  9 15  SIGMA 0.35
TORSION LIST 7 9 15 17 SIGMA 0.35
UMBRELLA CV 1 KAPPA  ${ks} AT ${val1}   
UMBRELLA CV 2 KAPPA  ${ks} AT ${val2}   
ENDMETA
EOF
$exe input.com >out
#exit
#grep "ENERGY:" out | tail -1 >>energy 
rm -rf out
awk '{printf("%12.6f %12.6f %12.6f\n",$1,$2+3.1415,$3+3.1415)}' COLVAR >CV_${incr}
echo "CV_${incr} ${val1b} ${val2b} ${ks} ${ks} ">>metadatafile 
nl=`wc -l CV_${incr} | awk '{print int($1/2)}' `
rm -rf COLVAR .l
#cp outconf.coor endconf_${incr}.pdb
mv outconf.coor start.pdb 
##exit
done
done
