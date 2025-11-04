This is a snakemake workflow covering index generation as well as quantitation of ONT cDNA data with oarfish, as well as differential expression analysis with dowstream tools.
Inputs to index generation are a genome fasta file and an annotation gtf file. Inputs to oarfish quantification are the transcriptome index file, as well as fastq files.
For the allele-specific mode, input is a bam file and a phased vcf file.

To create the conda environment for running the workflow, run:
```
conda env create -f env.yaml
```

To run the workflow, run:
```
conda activate oarfish_snakemake_env
snakemake -j 1 --use-envmodules --use-conda

```

The workflow is currently configured to run locally.



Citations:
oarfish:  Oarfish: Enhanced probabilistic modeling leads to improved accuracy in long read transcriptome quantification
Zahra Zare Jousheghani, Rob Patro
bioRxiv 2024.02.28.582591; doi: https://doi.org/10.1101/2024.02.28.582591 

edgeR/catchSalmon: Pedro L Baldoni, Yunshun Chen, Soroor Hediyeh-zadeh, Yang Liao, Xueyi Dong, Matthew E Ritchie, Wei Shi, Gordon K Smyth, Dividing out quantification uncertainty allows efficient assessment of differential transcript expression with edgeR, Nucleic Acids Research, Volume 52, Issue 3, 9 February 2024, Page e13, https://doi.org/10.1093/nar/gkad1167
