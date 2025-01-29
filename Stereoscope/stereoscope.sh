stereoscope run --sc_cnt /Users/tili/Desktop/CIIR/data/Wu_etal_2021_BRCA_scRNASeq/stereoscope/tnbc/BCSA1TumA2_tnbc_mtx.tsv \
--sc_labels /Users/tili/Desktop/CIIR/data/Wu_etal_2021_BRCA_scRNASeq/stereoscope/tnbc/BCSA1TumA2_tnbc_meta.tsv \
--st_cnt /Users/tili/Desktop/CIIR/data/ST_data/stereoscope/BCSA1TumA2.mtx.tsv \
-gl /Users/tili/Desktop/CIIR/data/ST_data/stereoscope/features.txt \
-stb 2048 -scb 2048 -lr 0.01 --gpu -ste 50000 -sce 50000 \
-o /Users/tili/Desktop/CIIR/results/stereoscope