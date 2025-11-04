rule oarfish_quant_allelic:
    input:
        index = "index/spliced_index",
        reads = "allelic_reads/{sample}_{allele}.fastq.gz"
    output:
        meta_json = "oarfish_output/{sample}_{allele}.meta_info.json",
        bootstrap = "oarfish_output/{sample}_{allele}.infreps.pq",
        quant = "oarfish_output/{sample}_{allele}.quant"
    params:
        prefix = "oarfish_output/{sample}_{allele}",
        woutdir = "oarfish_output",
        numBootstraps = config["num_bootstraps"]
    envmodules: "oarfish/0.9.0"
    threads: 16
    shell: """
          oarfish -j {threads} --reads {input.reads} --index {input.index} --seq-tech ont-cdna -o {params.prefix} --num-bootstraps {params.numBootstraps} --filter-group no-filters --model-coverage
          """ 


rule mod_json_allelic:
    input:
        meta_json_list = expand("oarfish_output/{sample}_{allele}.meta_info.json",sample=samples,allele=alleles),
        exons_fasta = "index/spliced.fa"
    output:
        expand("oarfish_output/{sample}_{allele}/aux_info/meta_info.json",sample=samples,allele=alleles)
    params:
        woutdir = "oarfish_output"
    conda: "envs/R.yaml"
    script: "../rscripts/mod_json.R"


rule pq_to_bootstrap_allelic:
    input:
        pq_list = expand("oarfish_output/{sample}_{allele}.infreps.pq",sample=samples,allele=alleles)
    output:
        temp(expand("oarfish_output/{sample}_{allele}/aux_info/bootstrap/bootstraps.tsv",sample=samples,allele=alleles)),
        expand("oarfish_output/{sample}_{allele}/aux_info/bootstrap/bootstraps.gz",sample=samples,allele=alleles)
    params:
        woutdir = "oarfish_output"
    conda: "envs/R.yaml"
    script: "../rscripts/pq_to_bootstrap.R"


rule rename_cols_allelic:
    input:
        quant_list = expand("oarfish_output/{sample}_{allele}.quant",sample=samples,allele=alleles)
    output:
        expand("oarfish_output/{sample}_{allele}/quant.sf",sample=samples,allele=alleles)
    params:
        woutdir = "oarfish_output"
    conda: "envs/R.yaml"
    script: "../rscripts/rename_cols.R"
