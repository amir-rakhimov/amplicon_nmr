#!/usr/bin/env bash
conda activate picrust2
../picrust2-2.5.1/scripts/picrust2_pipeline.py \
-s pooled-rep-seqs-trimmed-dada2-234-no-mit-no-chlor-NMR.fasta \
-i pooled-table-trimmed-dada2-234-no-mit-no-chlor-NMR.biom \
-o picrust2_out_NMR