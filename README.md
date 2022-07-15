Analysis of Saturation Mutagenesis STARR-seq data
===

## Running the pipeline
1. Install Snakemake
```sh
mamba create -c bioconda -c conda-forge --name snakemake snakemake
conda activate snakemake
```
3. Clone the repo
```sh
git clone https://github.com/pauldrinn/MuSTARRd.git
cd MuSTARRd
```
3. Modify `config/config.yaml`, `config/samples.tsv`, `config/units.tsv` and `profile/config.yaml` (documentation WIP)
4. Run the snakemake workflow with the following:
```py
snakemake --profile profile/ # use your own profile
```

Documentation is a WIP
