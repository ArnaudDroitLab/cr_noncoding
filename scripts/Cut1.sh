for i in 1 2 3 4 5
do
    R1="HI.0588.00"$i.$i"_R1.fastq.gz"
    R2="HI.0588.00"$i.$i"_R2.fastq.gz"
    outputdir=output_v2
    script=jobs/$i.cut1.sh
    mkdir -p $outputdir/cut1
    mkdir -p jobs
    
    cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A pyv-106-aa
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=8
module load mugqic/python/2.7.8
cutadapt -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG -o $outputdir/cut1/$i.R1.fastq.gz -p $outputdir/cut1/$i.R2.fastq.gz raw/$R1 raw/$R2
EOF
    workdir=`pwd`
    jobid=`qsub $script -o $script.stdout -e $script.stderr -d $workdir`
    echo $jobid > $script.jobid
done

