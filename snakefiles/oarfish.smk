rule oarfish_quant:
    input:
        index = "index/spliced_index",
        reads = os.path.join(config["input_dir"], "{sample}_fastq.gz")
    output:
        meta_json = "oarfish_output/{sample}.meta_info.json",
        bootstrap = "oarfish_output/{sample}.infreps.pq",
        quant = "oarfish_output/{sample}.quant"
    params:
        prefix = "oarfish_output/{sample}",
        woutdir = "oarfish_output",
        numBootstraps = config["num_bootstraps"]
    envmodules: "oarfish/0.9.0"
    threads: 16
    shell: """
          oarfish -j {threads} --reads {input.reads} --index {input.index} --seq-tech ont-cdna -o {params.prefix} --num-bootstraps {params.numBootstraps} --filter-group no-filters --model-coverage
          """ 


rule mod_json:
    input:
        meta_json_list = expand("oarfish_output/{sample}.meta_info.json",sample=samples),
        exons_fasta = "index/spliced.fa"
    output:
        expand("oarfish_output/{sample}/aux_info/meta_info.json",sample=samples)
    params:
        woutdir = "oarfish_output"
    conda: "envs/R.yaml"
    script: "../rscripts/mod_json.R"


rule pq_to_bootstrap:
    input:
        pq_list = expand("oarfish_output/{sample}.infreps.pq",sample=samples)
    output:
        temp(expand("oarfish_output/{sample}/aux_info/bootstrap/bootstraps.tsv",sample=samples)),
        expand("oarfish_output/{sample}/aux_info/bootstrap/bootstraps.gz",sample=samples)
    params:
        woutdir = "oarfish_output"
    conda: "envs/R.yaml"
    script: "../rscripts/pq_to_bootstrap.R"


rule rename_cols:
    input:
        quant_list = expand("oarfish_output/{sample}.quant",sample=samples)
    output:
        expand("oarfish_output/{sample}/quant.sf",sample=samples)
    params:
        woutdir = "oarfish_output"
    conda: "envs/R.yaml"
    script: "../rscripts/rename_cols.R"
