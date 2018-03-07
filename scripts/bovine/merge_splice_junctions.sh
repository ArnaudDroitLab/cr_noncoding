#!/bin/bash
#PBS -N reference_merge
#PBS -A pyv-106-aa
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=16
#PBS -o jobs/reference_merge.stdout
#PBS -e jobs/reference_merge.stderr

module load mugqic/star/2.5.1b
outputdir=output_v2

cat $outputdir/alignment_1stpass/*.out.tab |
    awk 'BEGIN {OFS="	"; strChar[0]="."; strChar[1]="+"; strChar[2]="-"} {if($5>0){print $1,$2,$3,strChar[$4]}}' |
    sort -k1,1h -k2,2n > $outputdir/alignment_1stPass/AllSamples.SJ.out.tab
mkdir -p $outputdir/reference.Merged
STAR --runMode genomeGenerate \
  --genomeDir $outputdir/reference.Merged \
  --genomeFastaFiles UMD3.1_ERCC/Bos_taurus.UMD3.1_ERCC.fa \
  --runThreadN 12 \
  --limitGenomeGenerateRAM 50000000000 \
  --sjdbFileChrStartEnd $outputdir/alignment_1stPass/AllSamples.SJ.out.tab \
  --limitIObufferSize 1000000000 \
  --sjdbOverhang 99