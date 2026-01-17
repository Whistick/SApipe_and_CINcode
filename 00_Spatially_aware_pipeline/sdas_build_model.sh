# SDAS

sample=$1

# Running spatial cell annotation with single cell RNA-seq reference
SDAS cellAnnotation scimilarityMakeRef -o ./scimilarity_ref --reference $sample.h5ad --label_key annotation \
--model_dir ./model_resolution_0.1 \  # resolution = 0.1 0.2 0.3 ...
--remove_tmp
