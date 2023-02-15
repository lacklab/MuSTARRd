Analysis of Saturation Mutagenesis STARR-seq data
===

## Running the pipeline
1. Install Snakemake
```sh
mamba create -c bioconda -c conda-forge -n snakemake snakemake
conda activate snakemake
```
2. Clone the repo
```sh
git clone https://github.com/lacklab/MuSTARRd.git
cd MuSTARRd
```
3. Modify `config/config.yaml`, `config/samples.tsv` and `config/units.tsv` (see [here](config/README.md))
4. Run the snakemake workflow with the following:
```py
snakemake --profile profile/ # use your own profile - modify profile/config.yaml
```

Documentation is a WIP
