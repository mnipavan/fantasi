#!/bin/bash
#PBS -l walltime=120:00:00
#PBS -l mpiprocs=8
#PBS -l ncpus=8
#PBS -l mem=16gb
#PBS -j oe

OUTDATE=`date`
echo "Started at $OUTDATE"
export SIMDIR=~/proj/fpe/fantasi/spintronic
export FNAME=large_field_2
export SCRIPT=${FNAME}.py
#export MPICOMMAND="mpirun -np 8 "
export MPICOMMAND=""
cd $TMPDIR
rsync -az ${SIMDIR}/${SCRIPT} ./
rsync -az ${SIMDIR}/HFields.py ./
rsync -az ${SIMDIR}/meshes ./
if [ -f ${SCRIPT} ]; then
    module load singularity
    singularity run -B $TMPDIR /opt/containers/ubuntu/18/fenics_v1_10.sif ${MPICOMMAND} python ./${SCRIPT}
fi
OUTDATE=`date`
echo "Simulation ended at $OUTDATE"
if [ -d ${FNAME}_results ]; then
    echo "Directory containing final result files found! Copying back to home directory"
    tar czf ./${FNAME}.tgz ./${FNAME}_results
    rsync -az ./${FNAME}.tgz ${SIMDIR}/.
fi
