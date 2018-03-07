for i in 1 2 3 4 5
do
    outputdir=output_v2
    script=jobs/$i.cut2.sh
    mkdir -p output_v2/cut2
    
    cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A pyv-106-aa
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=8
module load mugqic/python/2.7.8
cd `pwd`
cutadapt -m 40 -q 23 -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG -o $outputdir/cut2/$i.R1.fastq.gz -p $outputdir/cut2/$i.R2.fastq.gz $outputdir/cut1/$i.R1.fastq.gz $outputdir/cut1/$i.R2.fastq.gz
EOF
    workdir=`pwd`
    previousjob=`echo jobs/$i.cut1.sh.jobid`
    jobid=`qsub $script -o $script.stdout -e $script.stderr -d $workdir -W afterok:$previousjob`
    echo $jobid > $script.jobid
done

