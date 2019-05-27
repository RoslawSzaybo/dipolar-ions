#!/bin/bash

RUNFOLDER='\/net\/people\/plgofemski\/ions\/SrYb'
OUTPUTFOLDER='\/net\/people\/plgofemski\/ions\/SrYb\/05.23-omega_z-scan'
NGRANTID='coldmanybody2'
version='05.23'

n1="32"
n3="4"
n5="4"
j1="2"
j2="2"
# ~41.5k

# parameters for SrYb+
# m_{Sr} = 88 u
# m_{Yb} = 173 u
mass="261.0" # u
charge="1.0" # e
dipole="4.745" # D
# B = 0.0168 cm^{-1}
# B = 1.68 m^{-1}
B="503.6513" #MHz
#omega_z="0.16" #MHz
omega_rho="1.4" #MHz

NRAM="14000" #Mb

NQUEUEID="plgrid"
NTIME="23:59:59"

omegazs="0.01 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2"
for omega_z in $omegazs
do
	molecule="SrYb"
	NAMEX="${molecule}_w_z-${omega_z}"
	NAME="$NAMEX"

	sed -e "s/NAZWA/${NAMEX}/g" -e "s/JOB/${NAME}/g" -e "s/OUTPUTFOLDER/${OUTPUTFOLDER}/g" -e "s/RUNFOLDER/${RUNFOLDER}/g" -e "s/QUEUEID/${NQUEUEID}/g" -e "s/GRANTID/${NGRANTID}/g" -e "s/RAM/${NRAM}/g" -e "s/TIME/${NTIME}/g" -e "s/EXE/spectrum.${version}/g" /net/people/plgofemski/scripts/run_ions_new > run

	sbatch run "${n1} ${n3} ${n5} ${j1} ${j2} ${mass} ${charge} ${dipole} ${B} ${omega_rho} ${omega_z}"
	rm run
done
