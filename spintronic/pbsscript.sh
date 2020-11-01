#!/bin/bash
#PBS -l walltime=120:00:00
#PBS -l mpiprocs=8
#PBS -l ncpus=8
#PBS -l mem=16gb
#PBS -l vmem=16gb
#PBS -j oe

export SIMDIR=~/proj/fpe/fantasi/spintronic
export FNAME=large_field_0.py
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
    tar czf ./${FNAME}.tgz ./result_files
    rsync -az ./${FNAME}.tgz ${SIMDIR}/.
fi
