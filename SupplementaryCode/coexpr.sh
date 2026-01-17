# SDAS

sample=$1

./sdas-1.0.0/SDAS coexpress nest -i ./$sample.h5ad -o ./out --bin_size cellbin  --selected_genes top5000 --n_cpus 16 --seed 42 --hotspot_min_size 30 --hotspot_min_samples 4 --min_cells 100