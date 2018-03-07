#! /bin/bash

# Create directory to hold all temporary files
mkdir tmpSummary

for tissue in `cut -f 2 Input/Samples.txt | uniq`
do
    # List reps for this tissue
    for rep in `grep -F $tissue Input/Samples.txt | cut -f 3 | uniq`
    do
        echo Processing $tissue, replicate $rep
        
        # Add suffix to FPKM columns in header
        head -n 1 Output/Alignments/$tissue/$rep/cufflinks_merged_nolength/isoforms.fpkm_tracking > tmp_header
        awk '{OFS="\t"} {$9=$9 ".'$tissue.$rep'";$10=$10 ".'$tissue.$rep'";$11=$11 ".'$tissue.$rep'";$12=$12 ".'$tissue.$rep'";$13=$13 ".'$tissue.$rep'";} (1)' tmp_header > suffix_header
        
        # Sort file (without header)
        sed '1d' Output/Alignments/$tissue/$rep/cufflinks_merged_nolength/isoforms.fpkm_tracking | sort > sort_tracking
        
        # Concatenate header and sorted file
        cat suffix_header sort_tracking > tmpSummary/$tissue.$rep
        
        # Remove temporary files
        rm tmp_header suffix_header sort_tracking
    done
done

# Concatenate all columns from all files, and keep only the first header columns + the data columns
paste tmpSummary/* | cut -f 1-8,9-13,22-26,35-39,48-52,61-65,74-78,87-91,100-104,113-117,126-130,139-143,152-156,165-169,178-182,191-195,204-208 > Output/All.FPKM

# Clean up temporary directory
rm -r tmpSummary

 
                      
