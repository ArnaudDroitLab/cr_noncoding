outputdir=output_v2
for r1file in $outputdir/demultiplex/*.R1.fastq.gz
do
    samplename=`basename $r1file .R1.fastq.gz`
    r2file=$outputdir/demultiplex/$samplename.R2.fastq.gz
    script=jobs/$samplename.align-2.sh
    cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A pyv-106-aa
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=16
module load mugqic/python/2.7.8
module load mugqic/star/2.5.1b
mkdir -p $outputdir/alignment/$samplename
STAR --runMode alignReads \
  --genomeDir $outputdir/reference.Merged \
  --readFilesIn \
    $r1file \
    $r2file \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMstrandField intronMotif \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix $outputdir/alignment/$samplename \
  --outSAMattrRGline ID:"$samplename" 	PL:"ILLUMINA" 			SM:"$samplename" \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \  
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \  
  --chimSegmentMin 21  
EOF
    workdir=`pwd`
    #previousjob=`echo jobs/$i.cut1.sh.jobid`
    #jobid=`qsub $script -o $script.stdout -e $script.stderr -d $workdir -W afterok:$previousjob`
    jobid=`qsub $script -o $script.stdout -e $script.stderr -d $workdir`
    echo $jobid > $script.jobid

done
