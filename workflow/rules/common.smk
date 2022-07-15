import pandas as pd

samples = pd.read_csv(config["SAMPLES"], sep="\t", dtype=str, comment="#").set_index(
    "sample", drop=False
)
samples.index.names = ["sample_id"]

units = pd.read_csv(config["UNITS"], sep="\t", comment="#").set_index(
    ["sample", "unit"], drop=False
)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])


def get_raw_data(wildcards):
    return list(units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]])


def get_technical_replicates(wildcards):
    treps = list(units[units["sample"] == wildcards.sample]["unit"])
    libtype = samples.loc[wildcards.sample, "type"].lower()
    return expand(
        f"results/aligned_{libtype}/{wildcards.sample}_{{trep}}.bam", trep=treps
    )


outputs = ["results/DNA.RNA.counts.tsv.gz"]