#!/bin/bash
# This is where the bulk of the work is accomplished. For all tissues and all replicates,
# this script performs the following operations:
#   1. Clean libraries from adapter/linker/primer sequences and remove the Ns
#      used by the sequencing center to mask R2 adapter read-through.
#   2. Aligns the the cleaned libraries to the mm10 genome.
#   3. Assemble the transcripts using cufflinks.
# This involves lots of zipping/unzipping, concatenating, etc.

rm -r tmp
mkdir logs
mkdir -p Output/Libraries
mkdir -p Output/Alignments

# List unique tissues
for tissue in `cut -f 2 Input/Samples.txt | uniq`
do
    # List reps for this tissue
    for rep in `grep -F $tissue Input/Samples.txt | cut -f 3 | uniq`
    do
        echo Processing $tissue, replicate $rep
        
        mkdir tmp
        #fileList=`grep -P "$tissue\t$rep" Input/Samples.txt | cut -f 1`
        
        # Unzip all libraries simultaneously
        for file in `grep -P "$tissue\t$rep" Input/Samples.txt | cut -f 1`
        do
            gunzip Input/$file &
        done
        wait
        
        # Cut adapters from R1 and R2 libraries, simultaneously
        for R1File in `grep -P "$tissue\t$rep" Input/Samples.txt | cut -f 1 | sed 's/.gz//' | grep -F "R1"`
        do
            echo $R1File
            libName=`basename Input/$R1File`
            # gnl|uv|NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0
            cutadapt -f fastq -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG -e 0.15 -q 23 -o tmp/$libName.cut Input/$R1File > logs/$libName.cutadapt.log &
        done
        
        for R2File in `grep -P "$tissue\t$rep" Input/Samples.txt | cut -f 1 | sed 's/.gz//' | grep -F "R2"`
        do
            echo $R2File
            libName=`basename Input/$R2File`
            # gnl|uv|NGB00361.1:1-92 Illumina PCR Primer
            cutadapt -f fastq -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG -e 0.15 -q 23 -o tmp/$libName.cut Input/$R2File > logs/$libName.cutadapt.log &
        done
        wait

        # Rezip original files
        for file in `grep -P "$tissue\t$rep" Input/Samples.txt | cut -f 1 | sed 's/.gz//'`
        do
             gzip Input/$file &
        done
        wait
        
        # fastq-mcf will barf on empty sequences
        for f in tmp/*.cut
        do
            sed -i 's/^$/N/' $f &
        done
        wait
        
        
        # Remove Ns, small reads, in pairs
        for pair in `grep -P "$tissue\t$rep" Input/Samples.txt | cut -f 1 | sed 's/.R..fastq.gz//' | sort | uniq`
        do
            fastq-mcf -l 25 /dev/null \
                      tmp/$pair.R1.fastq.cut \
                      tmp/$pair.R2.fastq.cut \
                      -o tmp/$pair.R1.fastq.cut.mcf \
                      -o tmp/$pair.R2.fastq.cut.mcf > logs/$pair.mcf.log &
        done
        wait
        
        # Immediately remove temporary files to reduce hard disk load (Example: GV library creates 49Gb of "cut" libraries)
        rm tmp/*.cut
        
        # Remove adapters at the end of R2 reads
        for lib in `grep -P "$tissue\t$rep" Input/Samples.txt | cut -f 1 | sed 's/.gz//' | grep -F "R2"`
        do
            file=tmp/$lib.cut.mcf
            adapter=`grep -F $lib Input/Samples.txt | cut -f 4`
            revAdapter=`echo $adapter | tr "ACGT" "TGCA" | rev`
            perl ./Scripts/CutTag.pl $revAdapter < $file > $file.notag &
        done
        wait
        
        # Once again, immediately remove temporary files we are done with.
        rm tmp/*R2*.cut
        
        # Concatenate files and move them to the library folder.
        mkdir -p Output/Libraries/$tissue/$rep
        cat tmp/*R1.fastq.cut.mcf > Output/Libraries/$tissue/$rep/R1.fastq
        cat tmp/*R2.fastq.cut.mcf.notag > Output/Libraries/$tissue/$rep/R2.fastq
        
        # Clean up temporary files.
        rm -r tmp
      
        # Perform genomic alignment
        mkdir -p Output/Alignments/$tissue/$rep
        tophat2 --num-threads 6 --segment-mismatches 3  \
                --GTF Input/genes_ERCC_rDNA.gff \
                --transcriptome-index mm10Index2/transcriptome \
                -o Output/Alignments/$tissue/$rep \
                mm10Index2/mm10_rDNA_Index \
                Output/Libraries/$tissue/$rep/R1.fastq \
                Output/Libraries/$tissue/$rep/R2.fastq &> logs/$tissue.$rep.tophat2.log
                
        # Zip libraries 
        gzip Output/Libraries/$tissue/$rep/R1.fastq &
        gzip Output/Libraries/$tissue/$rep/R2.fastq &
        wait
        
        # Assemble transcripts
#                  --max-bundle-frags 60000000 \
        cufflinks -q -p 6 -L $tissue.$rep \
                  --GTF-guide Input/genes_ERCC_rDNA.gff \
                  -o Output/Alignments/$tissue/$rep/cufflinks \
                  Output/Alignments/$tissue/$rep/accepted_hits.bam  &> logs/$tissue.$rep.cufflinks.log
    done
done

# Send e-mail to notify that this process has completed.
echo "All libraries have been cleaned/aligned/assembled." | mail -s "Notification" eric.fournier.4@ulaval.ca