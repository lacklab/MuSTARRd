rule combine_RNA_counts:
    input:
        expand("results/{sample}.count.tsv", sample=samples[samples["type"] == "RNA"]["sample"]),
    output:
        "results/RNA.count.tsv.gz",
    shell:
        "bash workflow/scripts/multijoin.sh {input} | pigz > {output}"

rule combine_DNA_associations:
    input:
        expand("results/{sample}.association.tsv", sample=samples[samples["type"] == "DNA"]["sample"]),
    output:
        "results/DNA.association.tsv.gz",
    shell:
        """
        awk '
            FNR==1 && NR!=1 {{ while (/^BC/) getline; }}
            1 {{print}}
        ' {input} | pigz > {output}
        """

rule combine_RNA_and_DNA_counts:
    input:
        DNA="results/DNA.association.tsv.gz",
        RNA="results/RNA.count.tsv.gz",
    output:
        "results/DNA.RNA.counts.tsv.gz"
    shell:
        "join --header -t $'\t' <(pigz -dc {input.RNA}) <(pigz -dc {input.DNA}) > {output}"