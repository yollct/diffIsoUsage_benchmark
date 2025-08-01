[![DOI](https://zenodo.org/badge/572552561.svg)](https://doi.org/10.5281/zenodo.16692043)

# Differential usage isoform benchmarking analysis

## Simulation
### --dir - path where simulation file will be written to
### --type - "pair" or "single" ended
### --rep - number of replicates
### --depth - sequencing depth
### --noise - background noise (the ratio of genes that remain unchanged across conditions)

```bash ./simulation/mknewdir.sh --dir mame_of_dir --type pair --rep 4 --depth 5000000 --noise 0```

## Run analysis 
### --config - path to the config file
### --name - the name of directory where results will be written to

```bash main_script.sh --config /nfs/proj/is_benchmark/paired4_config.sh --name name_of_dir ```


