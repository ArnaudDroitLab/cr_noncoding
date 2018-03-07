#!/bin/bash

# Rename input files so that the R1/2 part is preceded by a period, as it makes name concatenation much easier later on.
# Shouldn't be needed anymore since I did not keep a copy of the original-name files.
# rename 's/_R1/.R1/' Input/*
# rename 's/_R2/.R2/' Input/*

# Download mouse genome if not present
if [ ! -d Mus_musculus ]
do
    wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
    tar -xzvf Mus_musculus_UCSC_mm10.tar.gz
done

# Build bowtie2 index using all chromosomes, ERCC sequences and rDNA sequence.
pushd Mus_musculus/UCSC/mm10/Sequence/Chromosomes/
bowtie2-build chr1.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chr10.fa,chr11.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17.fa,chr18.fa,chr19.fa,chrX.fa,chrM.fa,chrY.fa,/home/efournier/Projects/NonCodingMouse/Input/ERCC.fa,/home/efournier/Projects/NonCodingMouse/Input/mm-rDNA.fasta mm10Index
popd

mkdir mm10Index
mv Mus_musculus/UCSC/mm10/Sequence/Chromosomes/mm10Index* mm10Index/