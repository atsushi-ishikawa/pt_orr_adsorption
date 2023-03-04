#!/bin/sh 
#------ pjsub option --------#
#PJM -L rscgrp=n22240a
#PJM -L node=1
#PJM -L elapse=20:00:00 
#PJM -g n22240
#PJM -j
#------- Program execution -------#
NUM_NODES=${PJM_VNODES}
NUM_CORES=40
NUM_PROCS=`expr $NUM_NODES "*" $NUM_CORES`

module load intel

PRG=/lustre0/home/n22240/vasp/vasp.5.4.4/bin/vasp_std

LBL=$$

CUDIR=`pwd`
CIF=$CUDIR/"Pt.cif"
STEPS=200

vasp_script="${HOME}/ase/run_vasp.py"
echo 
echo "import os" > $vasp_script
echo "exitcode = os.system(\"mpiexec.hydra -n ${NUM_PROCS} ${PRG}\")" >> $vasp_script

python make_adsorbed_surf.py --cif=$CIF --worksubdir=$$ ${INP1} ${INP2} ${INP3}
python calc_adsorption_energy.py --worksubdir=$$ --steps=$STEPS 1> stdout_$LBL.txt 2> stderr_$LBL.txt
