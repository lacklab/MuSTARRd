rule merge_DNA_reads:
    input:
        get_raw_data,
    output:
        "results/aligned_dna/{sample}_{unit}_merged.fq.gz",
    conda:
        "../envs/association.yaml"
    shell:
        """
        bbmerge.sh in1={input[0]} in2={input[1]} out={output}
        """


rule map_DNA_reads:
    input:
        "results/aligned_dna/{sample}_{unit}_merged.fq.gz",
    output:
        coordsorted="results/aligned_dna/{sample}_{unit}.bam",
        qc="results/qc/total-reads.{sample}_{unit}",
    params:
        ref_bwa=config["REF_BWA_IDX"],
    threads: 64
    conda:
        "../envs/association.yaml"
    shell:
        """
        bc -l <<< $(zcat {input[0]} | wc -l)/4 > {output.qc};
        bwa mem -t {threads} {params.ref_bwa} {input} | samtools view -Sb | samtools sort -@{threads} -m 5G > {output.coordsorted} && samtools index -@{threads} {output.coordsorted};
        """


rule merge_DNA_treplicates:
    input:
        get_technical_replicates,
    output:
        "results/aligned_dna/{sample}.bam",
    conda:
        "../envs/association.yaml"
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


rule remove_DNA_indels:
    input:
        "results/aligned_dna/{sample}.bam",
    output:
        "results/aligned_dna/{sample}.noindels.bam",
    threads: 64
    conda:
        "../envs/association.yaml"
    shell:
        """
        samtools view -@{threads} -h {input} \
        | LC_ALL=C awk '$1 ~ "^@" || $6 !~ "I|D"' \
        | samtools view -@{threads} -Sb - > {output}

        samtools index -@{threads} {output}
        """


rule filter_DNA_reads:
    input:
        "results/aligned_dna/{sample}.noindels.bam",
    output:
        "results/filtered_dna/{sample}.filtered.bam",
    params:
        cap=config["PRIMERS"],
        umi1=config["UMI1_LEN"],
        umi2=config["UMI2_LEN"],
        insert_size=config["EXPECTED_TLEN"],
    threads: 64
    conda:
        "../envs/association.yaml"
    shell:
        """
        bash workflow/scripts/filter.sh {input} {output} {params} DNA {threads}
        """


rule reformat_DNA_reads_to_fastq:
    input:
        "results/filtered_dna/{sample}.filtered.bam",
    output:
        "results/filtered_dna/{sample}.filtered.fq.gz",
    conda:
        "../envs/association.yaml"
    shell:
        """
        reformat.sh in={input} out={output} primaryonly
        """


rule move_UMIs_seq2qname:
    input:
        "results/filtered_dna/{sample}.filtered.fq.gz",
    output:
        "results/umi_dna/{sample}.UMI.fq.gz",
    conda:
        "../envs/association.yaml"
    shell:
        """
        umi_tools extract --extract-method=regex --bc-pattern='(?P<umi_1>.{{6}}).+(?P<umi_2>.{{6}})$' --stdin={input} --stdout={output}
        """


rule map_qname_UMI:
    input:
        "results/umi_dna/{sample}.UMI.fq.gz",
    output:
        "results/umi_dna/{sample}.UMI.bam",
    threads: 64
    params:
        fixed_ref=config["REF_FIXED"],
    conda:
        "../envs/association.yaml"
    shell:
        """
        bwa mem -t {threads} {params.fixed_ref} {input} | samtools view -Sb > {output}
        """


rule get_UMIs:
    input:
        "results/umi_dna/{sample}.UMI.bam",
    output:
        "results/umi_dna/{sample}.umis.tsv.gz",
    conda:
        "../envs/association.yaml"
    shell:
        """        
        paste <(samtools view -F3852 {input} | cut -f1 | cut -f2 -d _ | cut -c-6) \
              <(samtools view -F3852 {input} | cut -f10 | cut -c -6) \
              <(samtools view -F3852 {input} | cut -f10 | rev | cut -c -6 | rev) \
              <(samtools view -F3852 {input} | cut -f1 | cut -f2 -d _ | cut -c7-) \
              <(samtools view -F3852 {input} | cut -f10) \
              <(samtools view -F3852 {input} | cut -f12) \
              <(samtools view -F3852 {input} | cut -f13) \
              | awk -v OFS="" '{{print $1,$2,$3,$4,"\t",$1,$4,"\t",$2,$3,"\t",$5,"\t",$6,"\t",$7}}' \
              | pigz > {output}
        """


rule mask_anchors:
    input:
        "results/umi_dna/{sample}.umis.tsv.gz",
    output:
        "results/umi_dna/{sample}.umis.masked.tsv.gz",
    conda:
        "../envs/association.yaml"
    shell:
        """
        paste <(pigz -dc {input} | cut -f1 | sed 's/./N/2; s/./N/23') \
              <(pigz -dc {input} | cut -f2 | sed 's/./N/2; s/./N/11') \
              <(pigz -dc {input} | cut -f3-6) \
              | parsort -k1,1 -S10G \
              | pigz > {output}
        """


rule count_UMIs:
    input:
        "results/umi_dna/{sample}.umis.masked.tsv.gz",
    output:
        "results/umi_dna/{sample}.umicounts.masked.tsv.gz",
    threads: 64
    conda:
        "../envs/association.yaml"
    shell:
        """
        pigz -dc {input} | cut -f1 | LC_ALL=C parsort -S10G --parallel={threads} | uniq -c | mawk -v OFS="\t" '{{print $2,$1}}' | pigz > {output}
        """


rule count_UMIs_per_region:
    input:
        "results/umi_dna/{sample}.umicounts.masked.tsv.gz",
    output:
        "results/umi_dna/{sample}.regioncounts.masked.tsv.gz",
    threads: 64
    conda:
        "../envs/association.yaml"
    shell:
        """
        pigz -dc {input} | LC_ALL=C mawk '{{print $2}}' | cut -c 7-18 | LC_ALL=C parsort -S10G --parallel={threads} | uniq -c | mawk -v OFS="\t" '{{print $2,$1}}' | pigz > {output}
        """


rule get_read_per_BC:
    input:
        umis="results/umi_dna/{sample}.umis.masked.tsv.gz",
        umicounts="results/umi_dna/{sample}.umicounts.masked.tsv.gz",
    output:
        temp("results/umi_dna/{sample}.RperBC.7.tsv.gz"),
    conda:
        "../envs/association.yaml"
    shell:
        """
        join -1 1 -2 1 <(pigz -dc {input.umis}) <(pigz -dc {input.umicounts}) | tr ' ' '\t' | pigz > {output}
        """


rule get_read_per_mut:
    input:
        "results/umi_dna/{sample}.RperBC.7.tsv.gz",
    output:
        out="results/umi_dna/{sample}.table.tsv.gz",
        tmp=temp("results/umi_dna/{sample}.tmp1.tsv"),
    conda:
        "../envs/association.yaml"
    shell:
        """
        pigz -dc {input} | awk '{{print $0,"\t",$1"-"$4}}' | parsort -k8,8 -S10G > {output.tmp}

        join -1 8 -2 1 {output.tmp} \
        <(pigz -dc {input} | awk '{{print $0,"\t",$1"-"$4}}' | cut -f8 | parsort -S10G | uniq -c | awk '{{printf("%s\t%d\\n", $2, $1);}}' | parsort -k1,1 -S10G) \
        | tr ' ' '\t' \
        | cut -f2- | pigz > {output.out}
        """


rule apply_UMI_filters:
    input:
        "results/umi_dna/{sample}.table.tsv.gz",
    output:
        "results/umi_dna/{sample}.umifiltered.tsv.gz",
    conda:
        "../envs/association.yaml"
    shell:
        "pigz -dc {input} | awk '($7 > 5) && ($8 / $7 > 0.75)' | pigz > {output}"


rule parse_mutations:
    input:
        "results/umi_dna/{sample}.umifiltered.tsv.gz",
    output:
        "results/{sample}.association.tsv.gz",
    conda:
        "../envs/association.yaml"
    shell:
        "pigz -dc {input} | uniq | awk -f workflow/scripts/parse_mutations.awk > {output}"
