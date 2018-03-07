#! /bin/bash
# Produce a merged reference assembly using all replicate-specific assemblies.

# List all assembly GTFs in a text file.
ls Output/Alignments/*/*/cufflinks/transcripts.gtf > AllTranscripts.txt

# Run cuffmerge
cuffmerge -g Input/genes_ERCC_rDNA.gff -o Output/MergedAssembly -p 6 AllTranscripts.txt &> logs/cuffmerge

# Clean temporary files.
rm AllTranscripts.txt