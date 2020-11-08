#!/bin/bash
#PBS -l walltime=120:00:00
#PBS -l mpiprocs=8
#PBS -l ncpus=8
#PBS -l mem=16gb
#PBS -l vmem=16gb
#PBS -j oe

export SIMDIR=~/proj/fpe/fantasi/spintronic
export FNAME=relaxation_3
export SCRIPT=${FNAME}.py
export MPICOMMAND="mpirun -np 8"
export OMP_NUM_THREADS=1
cd $TMPDIR
rsync -az ${SIMDIR}/${SCRIPT} ./
rsync -az ${SIMDIR}/rho_series.h5 ./
rsync -az ${SIMDIR}/HFields.py ./
rsync -az ${SIMDIR}/meshes ./
if [ -f ${SCRIPT} ]; then
    module load singularity openmpi pbspro
    singularity run -B $TMPDIR -B ${PBSPROPATH}/spool /opt/containers/ubuntu/18/fenics_v1_10.sif ${MPICOMMAND} python ./${SCRIPT}
fi
if [ -d ${FNAME}_results ]; then
    echo "Directory containing final result files found! Copying back to home directory"
    tar czf ./${FNAME}.tgz ./${FNAME}_results
    rsync -az ./${FNAME}.tgz ${SIMDIR}/.
fi
