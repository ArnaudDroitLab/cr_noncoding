#!/bin/bash
# Gathers statistics about libraries/alignments. Some statistics are obtained by reading 
# input/output files, while others are obtained through the parsing of log files. Output 
# is sent to files named Output/LibraryStatistics.txt and Output/AlignmentStatistics.txt,

# Write result file headers.
echo Tissue Replicate File NumberOfReads NumberOfBases AverageQuality NumberOfReadsTrimmed PercentageOfReadsTrimmed NumberOfReadsQualityTrimmed PercentageOfReadsQualityTrimmed NumberOfBasesTrimmed PercentageOfBasesTrimmed > Output/LibraryStatistics.txt
echo Tissue Replicate LeftKeptReads RightKeptReads LeftAverageQuality RightAverageQuality MappedReads MappedPercent ProperMapped ProperMappedPercent UnmappedReads MeanFragmentSize StdDevFragmentSize > Output/AlignmentStatistics.txt


# List unique tissues
for tissue in `cut -f 2 Input/Samples.txt.full | uniq`
do
    # List reps for this tissue
    for rep in `grep -F $tissue Input/Samples.txt.full | cut -f 3 | uniq`
    do
        # List all raw library files
        for file in `grep -P "$tissue\t$rep" Input/Samples.txt.full | cut -f 1 | sed 's/.gz$//'`
        do
            # Get library average quality
            if [ ! -e logs/$file.averagequality ]
            then
                gunzip -c Input/$file.gz | perl ./Scripts/AverageQuality.pl > logs/$file.averagequality
            fi
            averageQuality=`cat logs/$file.averagequality`
            
            # Get name of cutadapt log file
            cutFile=logs/$file.cutadapt.log
            
            # Extract values from cutadapt log
            numReads=`grep -F "Processed reads:" $cutFile | awk '{print $3}'`
            numBases=`grep -F "Processed bases:" $cutFile | awk '{print $3}'`
            numReadsTrimmed=`grep -F "Trimmed reads:" $cutFile | awk '{print $3}'`
            percentTrimmed=`grep -F "Trimmed reads:" $cutFile | awk '{print $4}' | tr -d "()%"`
            numQualityTrimmed=`grep -F "Quality-trimmed:" $cutFile | awk '{print $2}'`
            percentQualityTrimmed=`grep -F "Quality-trimmed:" $cutFile | awk '{print $6}' | tr -d "()%"`
            numBasesTrimmed=`grep -F "Trimmed bases:" $cutFile | awk '{print $3}'`
            percentBasesTrimmed=`grep -F "Trimmed bases:" $cutFile | awk '{print $7}' | tr -d "()%"`

            # Output information
            echo $tissue $rep $file $numReads $numBases $averageQuality $numReadsTrimmed $percentTrimmed $numQualityTrimmed $percentQualityTrimmed $numBasesTrimmed $percentBasesTrimmed >> Output/LibraryStatistics.txt
        done
        
        # Get stats about cleaned libraries
        tophatLog=logs/$tissue.$rep.tophat2.log
        leftKeptReads=`grep -F "left reads:" $tophatLog | awk '{print $7}'`
        rightKeptReads=`grep -F "right reads:" $tophatLog | awk '{print $7}'`
        
        libraryFolder=Output/Libraries/$tissue/$rep
        if [ ! -e logs/$tissue.$rep.R1.averagequality ]
        then
            `gunzip -c $libraryFolder/R1.fastq.gz | perl ./Scripts/AverageQuality.pl` > logs/$tissue.$rep.R1.averagequality
        fi
        if [ ! -e logs/$tissue.$rep.R2.averagequality ]
        then
            `gunzip -c $libraryFolder/R2.fastq.gz | perl ./Scripts/AverageQuality.pl` > logs/$tissue.$rep.R2.averagequality
        fi
        leftAverageQuality=`cat logs/$tissue.$rep.R1.averagequality`
        rightAverageQuality=`cat logs/$tissue.$rep.R2.averagequality`
        
        # Get alignment stats
        alignmentFolder=Output/Alignments/$tissue/$rep
        if [ ! -e logs/$tissue.$rep.accepted.stats ]
        then
            ~/bamUtil_1.0.7/bamUtil/bin/bam stats --basic --in $alignmentFolder/accepted_hits.bam &> logs/$tissue.$rep.accepted.stats
        fi
        mappedReads=`grep -e "^MappedReads" logs/$tissue.$rep.accepted.stats | cut -f 2`
        properMapped=`grep -e "^ProperPair\(\s\|(e.)\)" logs/$tissue.$rep.accepted.stats | cut -f 2`
        properMappedPercent=`grep -e "^ProperPair(%)" logs/$tissue.$rep.accepted.stats | cut -f 2`

        if [ ! -e logs/$tissue.$rep.unmapped.stats ]
        then        
            ~/bamUtil_1.0.7/bamUtil/bin/bam stats --basic --in $alignmentFolder/unmapped.bam &> logs/$tissue.$rep.unmapped.stats
        fi
        
        unmappedReads=`grep -e "^TotalReads\(\s\|(e.)\)" logs/$tissue.$rep.unmapped.stats | cut -f 2`
        mappedPercent=`echo "scale=2;$mappedReads / ($mappedReads + $unmappedReads) * 100" | bc -q`
        
        # Get fragment sizes
        cufflinksLog=logs/$tissue.$rep.cufflinks.log
        meanFragmentSize=`grep -F "Estimated Mean" $cufflinksLog | awk '{print $4}'`
        stdDevFragmentSize=`grep -F "Estimated Std Dev" $cufflinksLog | awk '{print $5}'`
        
        # Output information
        echo $tissue $rep $leftKeptReads $rightKeptReads $leftAverageQuality $rightAverageQuality $mappedReads $mappedPercent $properMapped $properMappedPercent $unmappedReads $meanFragmentSize $stdDevFragmentSize >> Output/AlignmentStatistics.txt
    done
done

echo "Statistics about libraries/alignments have been gathered." | mail -s "Notification" eric.fournier.4@ulaval.ca