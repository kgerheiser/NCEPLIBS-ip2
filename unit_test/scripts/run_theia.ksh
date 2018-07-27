#!/bin/ksh --login

#------------------------------------------------------------
# Script to run the unit test on Theia compute nodes.
#
# To run, type: 'qsub run_theia.ksh'
#
# Output is put in "unit_test.log"
#------------------------------------------------------------

#PBS -l nodes=1:ppn=1
#PBS -l vmem=5000M
#PBS -l walltime=0:15:00
#PBS -A fv3-cpu
#PBS -N ip2_unit_test
#PBS -q debug
#PBS -o unit_test.log
#PBS -e unit_test.log

set -x

export OMP_NUM_THREADS=1

rundir=$PBS_O_WORKDIR
cd $rundir

run_unit_test.ksh

exit 0
