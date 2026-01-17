# SDAS

sample=$1

SDAS cellAnnotation scimilarity -i $sample.h5ad -o outdir --bin_size 20 \
--model_dir ./model_resolution_0.2 \   # resolution = 0.1  0.2 0.3 ....
--reference_database scimilarity_ref
