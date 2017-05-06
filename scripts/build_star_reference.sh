#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -A pyv-106-aa
#PBS -o generate_star_index_ERCC.stdout
#PBS -e generate_star_index_ERCC.stderr
#PBS -V
#PBS -N generate_star_index_ERCC

module add mugqic/star/2.5.2a

cd /gs/project/pyv-106-aa/NonCodingSequencing
mkdir -p UMD3.1_ERCC/star_index_84

STAR --runMode genomeGenerate \
    --genomeDir UMD3.1_ERCC/star_index_84 \
    --genomeFastaFiles UMD3.1_ERCC/Bos_taurus.UMD3.1_ERCC.fa \
    --sjdbGTFfile UMD3.1_ERCC/Bos_taurus.UMD3.1_ERCC.gtf \
    --sjdbOverhang 99 \
    --runThreadN 6
    