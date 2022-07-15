rule map_RNA_reads:
    input:
        get_raw_data,
    output:
        coordsorted="results/aligned_rna/{sample}_{unit}.bam",
        qc="results/qc/total-reads.{sample}_{unit}",
    params:
        ref_bwa=config["REF_BWA_IDX"],
    threads: 64
    conda:
        "../envs/count.yaml"
    shell:
        """
        bc -l <<< $(zcat {input[0]} | wc -l)/4 > {output.qc}
        bwa mem -t {threads} {params.ref_bwa} {input} | samtools view -Sb | samtools sort -@{threads} -m 5G > {output.coordsorted} && samtools index -@{threads} {output.coordsorted}
        """


rule merge_RNA_treplicates:
    input:
        get_technical_replicates,
    output:
        "results/aligned_rna/{sample}.bam",
    conda:
        "../envs/count.yaml"
    shell:
        """
        if [[ "{input}" = *" "* ]]; then
            samtools merge -o {output} {input}
            samtools index {output}
        else
            mv {input} {output}
            samtools index {output}
        fi
        """


rule remove_RNA_indels:
    input:
        "results/aligned_rna/{sample}.bam",
    output:
        "results/aligned_rna/{sample}.noindels.bam",
    threads: 64
    conda:
        "../envs/count.yaml"
    shell:
        """
        samtools view -@{threads} -h {input} \
        | LC_ALL=C awk '$1 ~ "^@" || $6 !~ "I|D"' \
        | samtools sort -@{threads} -m 5G -nu \
        | samtools fixmate -m - - \
        | samtools sort -@{threads} -m 5G > {output}

        samtools index -@{threads} {output}
        """


rule filter_RNA_reads:
    input:
        "results/aligned_rna/{sample}.noindels.bam",
    output:
        "results/filtered_rna/{sample}.filtered.bam",
    params:
        cap=config["PRIMERS"],
        umi1=config["UMI1_LEN"],
        umi2=config["UMI2_LEN"],
        insert_size=config["EXPECTED_TLEN"],
    threads: 64
    conda:
        "../envs/count.yaml"
    shell:
        """
        bash workflow/scripts/filter.sh {input} {output} {params} RNA {threads}
        """


rule fix_and_sort_RNA_reads:
    input:
        "results/filtered_rna/{sample}.filtered.bam",
    output:
        "results/filtered_rna/{sample}.filtered.qname.bam",
    threads: 64
    conda:
        "../envs/count.yaml"
    shell:
        """
        samtools view -h -f2 -F3852 {input} \
        | samtools sort -@ {threads} -m 5G -nu \
        | samtools fixmate -m - - > {output}
        """


rule count_RNA_UMIs:
    input:
        "results/filtered_rna/{sample}.filtered.qname.bam",
    output:
        counts=temp("results/{sample}.count.tsv"),
        tmp=temp("results/{sample}.count.tmp"),
    threads: 64
    conda:
        "../envs/count.yaml"
    shell:
        """
        paste -d "" \
        <(samtools view -@ {threads} {input} | awk '$9 > 0' | cut -f10 | cut -c -12 | sed 's/./N/2') \
        <(samtools view -@ {threads} {input} | awk '$9 < 0' | cut -f10 | rev | cut -c -12 | sed 's/./N/2' | rev) > {output.tmp}

        parsort -S5G {output.tmp} | uniq -c | awk '{{printf("%s\t%d\\n", $2, $1);}}' | parsort -k1,1 -S5G > {output.counts}
        """
