# metadynamics 

 plumed            on
 plumedfile        metadyn.cfg 

# set system
coordinates start.pdb

# ffield
paratypecharmm on
parameters par_all27_prot_lipid.prm
switching      off  # Turn off the switching functions
exclude        scaled1-4
1-4scaling     0.833333  # =1/1.2, default is 1.0
scnb           2  # This is default
amber          on  # Specify this is AMBER force field
parmfile       diala.prmtop  # Input PARM file
#ambercoor      diala.inpcrd  # Input coordinate file

#approx
cutoff 22
pairlistdist 22
margin 0
stepspercycle 200000

# integrator
timestep 0.2

#output
outputenergies 1000
outputtiming 100000
binaryoutput no
dcdfreq 100000
outputname ./outconf

# protocol
temperature 300.0
#minimization on
langevin           on
langevinDamping     8
langevinTemp        300

#script 
numsteps     50000

seed		791064881

#constraint 
#constraints on  
#consref start.pdb
#conskfile start.pdb
#conskcol O
#constraintScaling 5
