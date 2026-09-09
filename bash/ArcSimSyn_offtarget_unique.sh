#!/bin/bash

#SBATCH --job-name=ArcSimSyn_offtarget_unique
#SBATCH --output=slurm/ArcSimSyn_offtarget_unique/%j.out
#SBATCH --error=slurm/ArcSimSyn_offtarget_unique/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=96:00:00

# Nextflow
pixi run nextflow run nextflow/pipelines/ArcSimSyn_offtarget_unique.nf \
    -c nextflow/configs/ArcSimSyn_offtarget_unique.config \
    -work-dir /mnt/ceph/groups_hdd/SCCGroup/computational_systems_biology_group/dkatzlberger/ArcSimulationSynthetic/work/ArcSimSyn_offtarget_unique \
    -resume