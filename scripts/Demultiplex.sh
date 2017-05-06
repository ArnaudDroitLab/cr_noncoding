# Loop over all four libraries.
for i in 1 2 3 4 5
do
    # Generate index files for fastq-multx
    grep -e "^$i" input/Samples.txt | awk '{ print $2"."$3"\t"$5}' > input/$i.sample.index

    outputdir=output_v2
    mkdir -p $outputdir/demultiplexed
    
    # Generate and submit demultiplexing job.    
    script=jobs/$i.demultiplex.sh    
    cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A pyv-106-aa
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=8
module load mugqic/python/2.7.8
cd `pwd`
module load mugqic/ea-utils/1.1.2-537
fastq-multx input/$i.sample.index $outputdir/cut2/$i.R1.fastq.gz $outputdir/cut2/$i.R2.fastq.gz -o $outputdir/demultiplexed/%.R1.fastq.gz $outputdir/demultiplexed/%.R2.fastq.gz
EOF
    workdir=`pwd`
    previousjob=`echo jobs/$i.cut1.sh.jobid`
    jobid=`qsub $script -o $script.stdout -e $script.stderr -d $workdir -W afterok:$previousjob`
    echo $jobid > $script.jobid
done

