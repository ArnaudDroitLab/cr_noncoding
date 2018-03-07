#!/bin/bash
# Re-calculate FPKMs on the reference assembly so that all transcripts can be directly
# compared across stages/replicates.

# Loop over all tissues, replicates.
for tissue in `cut -f 2 Input/Samples.txt | uniq`
do
    # List reps for this tissue
    for rep in `grep -F $tissue Input/Samples.txt | cut -f 3 | uniq`
    do
        echo Processing $tissue, replicate $rep

        # Run cufflinks using merged assembly.
        cufflinks -q -p 6 -L $tissue.$rep \
                  --GTF Output/MergedAssembly/merged.gtf \
                  --no-effective-length-correction \
                  -o Output/Alignments/$tissue/$rep/cufflinks_merged_nolength \
                  Output/Alignments/$tissue/$rep/accepted_hits.bam  &> logs/$tissue.$rep.cufflinks.on.merged.logs.nolength
    done
done

# Send e-mail to notify that this process has completed.
echo "Cufflinks is done re-calculating FPKMs against the reference assembly" | mail -s "Notification" eric.fournier.4@ulaval.ca