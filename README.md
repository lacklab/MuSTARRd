Analysis of Saturation Mutagenesis STARR-seq data
===

## Running the pipeline
1. Clone the repo
```sh
git clone https://github.com/pauldrinn/MuSTARRd.git
cd MuSTARRd
```
2. Change `config/config.yaml`, `config/samples.tsv`, `config/units.tsv` and `profile/config.yaml` (documentation WIP)
3. Run snakemake with the command:
```py
snakemake --profile profile/ 
```

Documentation is a WIP
