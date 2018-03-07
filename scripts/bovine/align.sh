outputdir=output_v2
for r1file in $outputdir/demultiplex/*.R1.fastq.gz
do
    samplename=`basename $r1file .R1.fastq.gz`
    r2file=$outputdir/demultiplex/$samplename.R2.fastq.gz
    script=jobs/$samplename.align-1.sh
    cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A pyv-106-aa
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=16
module load mugqic/python/2.7.8
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/$samplename && \
STAR --runMode alignReads \
  --genomeDir /gs/project/pyv-106-aa/NonCodingSequencing/UMD3.1_ERCC/star_index_84 \
  --readFilesIn \
    $r1file \
    $r2file \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix $outputdir/alignment_1stPass/$samplename \
  --outSAMattrRGline ID:"$samplename" 	PL:"ILLUMINA" 			SM:"$samplename" \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
EOF
    workdir=`pwd`
    #previousjob=`echo jobs/$i.cut1.sh.jobid`
    #jobid=`qsub $script -o $script.stdout -e $script.stderr -d $workdir -W afterok:$previousjob`
    jobid=`qsub $script -o $script.stdout -e $script.stderr -d $workdir`
    echo $jobid > $script.jobid

done
