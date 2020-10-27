#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l mpiprocs=8
#PBS -l ncpus=8
#PBS -j oe

export SIMDIR=/home/finesse/xfong/proj/fpe/fantasi/spintronic
export FNAME=relaxation_1.py
cd $TMPDIR
rsync -az ${SIMDIR}/${FNAME} ./
rsync -az ${SIMDIR}/HFields.py ./
rsync -az ${SIMDIR}/meshes ./
if [ -f ${FNAME} ]; then
    module load singularity
    singularity run -B $TMPDIR /opt/containers/ubuntu/18/fenics_v1_10.sif python ./${FNAME}
fi
if [ -d result_files ]; then
    echo "Directory containing final result files found! Copying back to home directory"
    rsync -az ./result_files ${SIMDIR}/.
fi
