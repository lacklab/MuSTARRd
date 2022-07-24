In order to configure your analysis, make changes to `config.yaml`.

1. `SAMPLES`: A tab-separated file with the following example should be provided to specify the samples:

| sample  | condition | type |
|---------|-----------|------|
| ETOH.R1 | EtOH      | RNA  |
| ETOH.R2 | EtOH      | RNA  |
| ETOH.R3 | EtOH      | RNA  |
| DHT.R1  | DHT       | RNA  |
| DHT.R2  | DHT       | RNA  |
| DHT.R3  | DHT       | RNA  |
| INPUT   |           | DNA  |

sample: Sample name

condition: Treatment condition (blank if none; blank for DNA)

type: Sample type (RNA or DNA)

2. `UNITS`: A tab-separated file with the following example should be provided to specify the units (replicates or lanes):

| sample  | unit | fq1                               | fq2                               |
|---------|------|-----------------------------------|-----------------------------------|
| ETOH.R1 | 1    | reads/EM_LNCaP_R1_EtOH_L1_1.fq.gz | reads/EM_LNCaP_R1_EtOH_L1_2.fq.gz |
| ETOH.R2 | 1    | reads/EM_LNCaP_R2_EtOH_L1_1.fq.gz | reads/EM_LNCaP_R2_EtOH_L1_2.fq.gz |
| ETOH.R3 | 1    | reads/EM_LNCaP_R3_EtOH_L1_1.fq.gz | reads/EM_LNCaP_R3_EtOH_L1_2.fq.gz |
| ETOH.R3 | 2    | reads/EM_LNCaP_R3_EtOH_L2_1.fq.gz | reads/EM_LNCaP_R3_EtOH_L2_2.fq.gz |
| DHT.R1  | 1    | reads/EM_LNCaP_R1_DHT_L1_1.fq.gz  | reads/EM_LNCaP_R1_DHT_L1_2.fq.gz  |
| DHT.R2  | 1    | reads/EM_LNCaP_R2_DHT_L1_1.fq.gz  | reads/EM_LNCaP_R2_DHT_L1_2.fq.gz  |
| DHT.R3  | 1    | reads/EM_LNCaP_R3_DHT_L1_1.fq.gz  | reads/EM_LNCaP_R3_DHT_L1_2.fq.gz  |
| INPUT   | 1    | reads/NL18_L2_1.fq.gz             | reads/NL18_L2_2.fq.gz             |

sample: Sample name (same as in samples.tsv)

unit: Unit no

fq1: Path to the 1st FASTQ file

fq2: Path to the 2nd FASTQ file

3. `PRIMERS`: A tab-separated file with the following specificiations (and example row) should be provided to specify the regions of analysis:

| Region name         | Forward primer  | Reverse primer         | Chromosome | Start   | End     |
|---------------------|-----------------|------------------------|------------|---------|---------|
| overlapped_read_114 | AGCGCGGCTTAGTGA | TACCAGGAGACTATTTCCAACA | chr8       | 6456903 | 6457213 |

This file shouldn't have a header. The primers should be 5' to 3' and the positions should be BED-like (0-based).

4. `REF`
   1. `FA`: Path to the FASTA file of the reference genome
   2. `BWA_IDX`: Path to the BWA index files (prefix)
   3. `FIXED`: Path to the FASTA file matching the wild type plasmids (or the same as `FA`)

5. `SEQ`
   1. `UMI1_LEN`: Length of the 5' UMI
   2. `UMI2_LEN`: Length of the 3' UMI
   3. `EXPECTED_TLEN`: Template (insert) length of the plasmids (excluding UMIs but including primers)