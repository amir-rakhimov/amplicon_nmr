#!/usr/bin/env bash
conda activate picrust2
trunc_len=234
custom_classes=("NMR" "SPFmouse" "FukomysDamarensis" "hare" "rabbit" "spalax"  "pvo")
custom_classes_filename=$(IFS=-; echo "${custom_classes[*]}")
custom_classes_picrust=$(IFS=_; echo "${custom_classes[*]}")
../picrust2-2.5.2/scripts/picrust2_pipeline.py \
-s pooled-single-filtered-rep-seqs-trimmed-dada2-${trunc_len}-$custom_classes_filename.fasta \
-i pooled-single-filtered-table-trimmed-dada2-${trunc_len}-$custom_classes_filename.biom \
-o picrust2_out_${custom_classes_picrust}