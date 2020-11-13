#!/bin/bash
#PBS -l walltime=120:00:00
#PBS -l mpiprocs=16
#PBS -l ncpus=16
#PBS -l mem=32gb
#PBS -j oe

OUTDATE=`date`
echo "Started at $OUTDATE"
export SIMDIR=~/proj/fpe/fantasi/spintronic
export FNAME=relaxation_05
export SCRIPT=${FNAME}.py
export MPICOMMAND="mpirun -np 16"
export OMP_NUM_THREADS=1
cd $TMPDIR
rsync -az ${SIMDIR}/${SCRIPT} ./
rsync -az ${SIMDIR}/HFields.py ./
rsync -az ${SIMDIR}/rho_series.h5 ./
rsync -az ${SIMDIR}/meshes ./
if [ -f ${SCRIPT} ]; then
    module load singularity cuda/10.0 gdrcopy libnl libibverbs openucx openmpi pbspro
    singularity run -B $TMPDIR -B ${PBSPROPATH}/spool /opt/containers/ubuntu/18/fenics_v1_10.sif ${MPICOMMAND} python ./${SCRIPT} \
    1> >(tee ${PBS_O_WORKDIR}/${PBS_JOBID}.${PBS_JOBNAME}.stdout) 2> >(tee ${PBS_O_WORKDIR}/${PBS_JOBID}.${PBS_JOBNAME}.stderr)
fi
OUTDATE=`date`
echo "Simulation ended at $OUTDATE"
if [ -d ${FNAME}_results ]; then
    echo "Directory containing final result files found! Copying back to home directory"
    tar czf ./${FNAME}.tgz ./${FNAME}_results
    rsync -az ./${FNAME}.tgz ${SIMDIR}/.
fi
