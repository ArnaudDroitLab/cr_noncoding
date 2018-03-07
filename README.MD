## Synopsis

This is the code for the analysis of Claude Robert's RNA-Seq libraries of ten different embryonic bovine tissues.

## Libraries (Bovine)

The libraries were generated using a mix of random and dT oligos, which should provide a fair coverage of the full length of transcripts, rather than just their 3' ends. Genes lacking a poly-A tail should also be represented. All libraries were ribo-depleted, and ERCC spike-ins were added in the extraction buffer, insuring that they would go through all subsequent processing steps. The libraries were subjected to a linear amplification, then were barcoded and pooled into 5 sequencing lanes.

| Library | Stage               | Rep | Tag_Name | Tag_Seq |
|---------|---------------------|-----|----------|---------|
| 1       | GV                  | a   | L2-BC1   | ACCC    |
| 1       | GV                  | b   | L2-BC2   | CGTA    |
| 1       | GV                  | c   | L2-BC3   | GAGT    |
| 1       | MII                 | a   | L2-BC4   | TTAG    |
| 1       | MII                 | b   | L2-BC5   | AGGG    |
| 1       | MII                 | c   | L2-BC6   | CCAT    |
| 2       | 1-cells             | a   | L2-BC7   | GTCA    |
| 2       | 1-cells             | b   | L2-BC8   | TATC    |
| 2       | 1-cells             | c   | L2-BC1   | ACCC    |
| 2       | 2-cells             | a   | L2-BC2   | CGTA    |
| 2       | 2-cells             | b   | L2-BC3   | GAGT    |
| 2       | 2-cells             | c   | L2-BC4   | TTAG    |
| 3       | 4-cells             | a   | L2-BC5   | AGGG    |
| 3       | 4-cells             | b   | L2-BC6   | CCAT    |
| 3       | 4-cells             | c   | L2-BC7   | GTCA    |
| 3       | 8-cells-early       | a   | L2-BC8   | TATC    |
| 3       | 8-cells-early       | b   | L2-BC1   | ACCC    |
| 3       | 8-cells-early       | c   | L2-BC2   | CGTA    |
| 4       | 8-cells-late        | a   | L2-BC3   | GAGT    |
| 4       | 8-cells-late        | b   | L2-BC4   | TTAG    |
| 4       | 8-cells-late        | c   | L2-BC5   | AGGG    |
| 4       | Morula              | a   | L2-BC6   | CCAT    |
| 4       | Morula              | b   | L2-BC7   | GTCA    |
| 4       | Morula              | c   | L2-BC8   | TATC    |
| 5       | Young-blastocyst    | a   | L2-BC1   | ACCC    |
| 5       | Young-blastocyst    | b   | L2-BC2   | CGTA    |
| 5       | Young-blastocyst    | c   | L2-BC3   | GAGT    |
| 5       | Expanded-blastocyst | a   | L2-BC4   | TTAG    |
| 5       | Expanded-blastocyst | b   | L2-BC5   | AGGG    |
| 5       | Expanded-blastocyst | c   | L2-BC6   | CCAT    |

## Processing (Bovine)

In orders, the following steps must be carried out:

1. Adapters need to be [trimmed](scripts/cut1.sh) [twice](scripts/cut1.sh) to remove the Illumina primers:

```
>gnl|uv|NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0 (Oligonucleotide sequence copyright 2007-2012 Illumina, Inc. All rights reserved.)
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>gnl|uv|NGB00360.1:1-58 Illumina PCR Primer (Oligonucleotide sequence copyright 2007-2012 Illumina, Inc. All rights reserved.)
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
```

2. The libraries must be [demultiplexed](scripts/demultiplex.sh)
3. Readthrough multiplexing barcodes should be removed. (TODO)
4. The libraries must be aligned to a reference including both the bovine genome (bosTau8/UMD3.1.1) and the [ERCC sequences](input/ERCC.fa) using HiSAT2. (TODO)
5. The novel transcripts from each library must be identified using StringTIE. (TODO)
6. The StringTIE transcripts must be merged. (TODO)
7. The merged transcripts must be requantified. (HTseq-count?)
8. The raw quantifications must be normalized using the ERCC abundances.
9. The transcripts must be annotated. (TODO)

## Libraries (Mouse)

| File                               | Tissue       |Rep| Tag  |
|------------------------------------|--------------|---|------|
| 1_cell-A_ACCC_L002.R1.fastq.gz     | 1c           | A | ACCC |
| 1_cell-A_ACCC_L002.R2.fastq.gz     | 1c           | A | ACCC |
| 1_cell-A_ACCC_L003.R1.fastq.gz     | 1c           | A | ACCC |
| 1_cell-A_ACCC_L003.R2.fastq.gz     | 1c           | A | ACCC |
| 1_cell-B_CGTA_L002.R1.fastq.gz     | 1c           | B | CGTA |
| 1_cell-B_CGTA_L002.R2.fastq.gz     | 1c           | B | CGTA |
| 1_cell-B_CGTA_L003.R1.fastq.gz     | 1c           | B | CGTA |
| 1_cell-B_CGTA_L003.R2.fastq.gz     | 1c           | B | CGTA |
| 1_cell-C_GAGT_L002.R1.fastq.gz     | 1c           | C | GAGT |
| 1_cell-C_GAGT_L002.R2.fastq.gz     | 1c           | C | GAGT |
| 1_cell-C_GAGT_L003.R1.fastq.gz     | 1c           | C | GAGT |
| 1_cell-C_GAGT_L003.R2.fastq.gz     | 1c           | C | GAGT |
| 2_cells-A_TTAG_L004.R1.fastq.gz    | 2c           | A | TTAG |
| 2_cells-A_TTAG_L004.R2.fastq.gz    | 2c           | A | TTAG |
| 2_cells-A_TTAG_L005.R1.fastq.gz    | 2c           | A | TTAG |
| 2_cells-A_TTAG_L005.R2.fastq.gz    | 2c           | A | TTAG |
| 2_cells-B_AGGG_L004.R1.fastq.gz    | 2c           | B | AGGG |
| 2_cells-B_AGGG_L004.R2.fastq.gz    | 2c           | B | AGGG |
| 2_cells-B_AGGG_L005.R1.fastq.gz    | 2c           | B | AGGG |
| 2_cells-B_AGGG_L005.R2.fastq.gz    | 2c           | B | AGGG |
| 2_cells-C_CCAT_L004.R1.fastq.gz    | 2c           | C | CCAT |
| 2_cells-C_CCAT_L004.R2.fastq.gz    | 2c           | C | CCAT |
| 2_cells-C_CCAT_L005.R1.fastq.gz    | 2c           | C | CCAT |
| 2_cells-C_CCAT_L005.R2.fastq.gz    | 2c           | C | CCAT |
| 4_cells_GAGT_L006.R1.fastq.gz      | 4c           | A | GAGT |
| 4_cells_GAGT_L006.R2.fastq.gz      | 4c           | A | GAGT |
| 4_cells_GAGT_L007.R1.fastq.gz      | 4c           | A | GAGT |
| 4_cells_GAGT_L007.R2.fastq.gz      | 4c           | A | GAGT |
| 4_cells_GAGT_L008.R1.fastq.gz      | 4c           | A | GAGT |
| 4_cells_GAGT_L008.R2.fastq.gz      | 4c           | A | GAGT |
| 8_cells-A_TTAG_L002.R1.fastq.gz    | 8c           | A | TTAG |
| 8_cells-A_TTAG_L002.R2.fastq.gz    | 8c           | A | TTAG |
| 8_cells-A_TTAG_L003.R1.fastq.gz    | 8c           | A | TTAG |
| 8_cells-A_TTAG_L003.R2.fastq.gz    | 8c           | A | TTAG |
| 8_cells-B_AGGG_L002.R1.fastq.gz    | 8c           | B | AGGG |
| 8_cells-B_AGGG_L002.R2.fastq.gz    | 8c           | B | AGGG |
| 8_cells-B_AGGG_L003.R1.fastq.gz    | 8c           | B | AGGG |
| 8_cells-B_AGGG_L003.R2.fastq.gz    | 8c           | B | AGGG |
| Blast_GTCA_L006.R1.fastq.gz        | Blast        | A | GTCA |
| Blast_GTCA_L006.R2.fastq.gz        | Blast        | A | GTCA |
| Blast_GTCA_L007.R1.fastq.gz        | Blast        | A | GTCA |
| Blast_GTCA_L007.R2.fastq.gz        | Blast        | A | GTCA |
| Blast_GTCA_L008.R1.fastq.gz        | Blast        | A | GTCA |
| Blast_GTCA_L008.R2.fastq.gz        | Blast        | A | GTCA |
| Early_blast_CCAT_L006.R1.fastq.gz  | Early_blast  | A | CCAT |
| Early_blast_CCAT_L006.R2.fastq.gz  | Early_blast  | A | CCAT |
| Early_blast_CCAT_L007.R1.fastq.gz  | Early_blast  | A | CCAT |
| Early_blast_CCAT_L007.R2.fastq.gz  | Early_blast  | A | CCAT |
| Early_blast_CCAT_L008.R1.fastq.gz  | Early_blast  | A | CCAT |
| Early_blast_CCAT_L008.R2.fastq.gz  | Early_blast  | A | CCAT |
| Expand_blast_TATC_L006.R1.fastq.gz | Expand_blast | A | TATC |
| Expand_blast_TATC_L006.R2.fastq.gz | Expand_blast | A | TATC |
| Expand_blast_TATC_L007.R1.fastq.gz | Expand_blast | A | TATC |
| Expand_blast_TATC_L007.R2.fastq.gz | Expand_blast | A | TATC |
| Expand_blast_TATC_L008.R1.fastq.gz | Expand_blast | A | TATC |
| Expand_blast_TATC_L008.R2.fastq.gz | Expand_blast | A | TATC |
| GV_CGTA_L006.R1.fastq.gz           | GV           | A | CGTA |
| GV_CGTA_L006.R2.fastq.gz           | GV           | A | CGTA |
| GV_CGTA_L007.R1.fastq.gz           | GV           | A | CGTA |
| GV_CGTA_L007.R2.fastq.gz           | GV           | A | CGTA |
| GV_CGTA_L008.R1.fastq.gz           | GV           | A | CGTA |
| GV_CGTA_L008.R2.fastq.gz           | GV           | A | CGTA |
| Morula-A_GTCA_L004.R1.fastq.gz     | Morula       | A | GTCA |
| Morula-A_GTCA_L004.R2.fastq.gz     | Morula       | A | GTCA |
| Morula-A_GTCA_L005.R1.fastq.gz     | Morula       | A | GTCA |
| Morula-A_GTCA_L005.R2.fastq.gz     | Morula       | A | GTCA |
| Morula-B_TATC_L004.R1.fastq.gz     | Morula       | B | TATC |
| Morula-B_TATC_L004.R2.fastq.gz     | Morula       | B | TATC |
| Morula-B_TATC_L005.R1.fastq.gz     | Morula       | B | TATC |
| Morula-B_TATC_L005.R2.fastq.gz     | Morula       | B | TATC |
| UF_ACCC_L006.R1.fastq.gz           | UF           | A | ACCC |
| UF_ACCC_L006.R2.fastq.gz           | UF           | A | ACCC |
| UF_ACCC_L007.R1.fastq.gz           | UF           | A | ACCC |
| UF_ACCC_L007.R2.fastq.gz           | UF           | A | ACCC |
| UF_ACCC_L008.R1.fastq.gz           | UF           | A | ACCC |
| UF_ACCC_L008.R2.fastq.gz           | UF           | A | ACCC |

## Processing (Mouse)