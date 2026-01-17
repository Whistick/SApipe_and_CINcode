# SDAS

samp=$1

./sdas-1.0.0/SDAS cellAnnotation tangram -i $sample.h5ad -o outdir --reference sc.h5ad --bin_size 20 --label_key annotation2 --filter_rare_cell 0
