#!/bin/sh

$1 pdb2gmx -ignh -ff amber99sb-ildn -f $2.pdb -o $2.gro -p $2.top -water tip3p
$1 editconf -f $2.gro -o $2-PBC.gro -bt triclinic -d 3.0
$1 grompp -f $3/mdp/em0.mdp -c $2-PBC.gro -p $2.top -o em-vac.tpr -maxwarn 2
$1 mdrun -v -deffnm em-vac

# To obtain initial energies, call this
echo 5 8 |  $1 energy -f em-vac.edr -s em-vac.tpr -o energies.xvg

#plot "./scwrl_random_2KB8-bkbone-rmsd.xvg" using 1:2 with lines

#Making video in pymol
#dss 
#show cartoon 
#viewport 640,800 
#set ray_trace_frames,1 
#mpng frame_.png
