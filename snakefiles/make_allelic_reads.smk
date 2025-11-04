rule haplotag_bam:
    input:
        bam = os.path.join(config["input_dir"], "{sample}.bam"),
        vcf = config["phasedVcf"],
        genome_fasta = config["genome_fasta"]
    output:
        haplotagged_bam = "haplotagged_bams/{sample}.bam",
        hap_list = "haplotagged_bams/{sample}_haplotype_list.tsv"
    envmodules: "whatshap/2.3"
    shell: """
        whatshap haplotag --ignore-read-groups -o {output.haplotagged_bam} --reference {input.genome_fasta} --output-threads=4 --output-haplotag-list={output.hap_list} {input.vcf} {input.bam}
        """


rule split_bam:
    input:
        haplotagged_bam = "haplotagged_bams/{sample}.bam",
        hap_list = "haplotagged_bams/{sample}_haplotype_list.tsv"
    output:
        allelic_bams = expand("allelic_bams/{{sample}}_{allele}.bam",allele=alleles)
    params:
        h1 = "allelic_bams/{sample}_h1.bam",
        h2 = "allelic_bams/{sample}_h2.bam"
    envmodules: "whatshap/2.3"
    shell: """
        whatshap split  --output-h1 {params.h1} --output-h2 {params.h2} {input.haplotagged_bam} {input.hap_list}
        """

rule namesort_bam:
    input:
        allelic_bam = "allelic_bams/{sample}_{allele}.bam"
    output:
        namesorted_bam = "namesorted_bams/{sample}_{allele}.bam"
    threads: 4
    envmodules: "samtools/1.21"
    shell: """
        samtools sort -n {input.allelic_bam} -o {output.namesorted_bam} -@ {threads} -m 2G
        """

rule make_reads:
    input:
        namesorted_bam = "namesorted_bams/{sample}_{allele}.bam",
    output:
        reads = "allelic_reads/{sample}_{allele}.fastq.gz"
    threads: 4
    envmodules: "samtools/1.21"
    shell: """
        samtools fastq -F 0x100  -n -@ {threads} {input.namesorted_bam} | gzip -c - > {output.reads}
        """
