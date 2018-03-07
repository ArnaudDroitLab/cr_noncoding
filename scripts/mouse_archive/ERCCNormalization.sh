grep -e "^ERCC" Output/MergedAssembly/merged.gtf | awk '{print $1 "\t" $12}' | tr -d '";' > Output/ERCC-Transcript.map

Rscript Scripts/ERCCNormalization.R