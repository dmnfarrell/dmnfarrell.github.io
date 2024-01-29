---
layout: post
title:  "kSNP for alignment free phylogenetics"
date:   2024-01-22 11:40:00
categories: bioinformatics
tags: [python,genomics]
thumbnail: /img/tree_branches.png
---

## Background

<div style="width: 200px; float:right;">
 <a href="/img/tree_branches.png"> <img src="/img/tree_branches.png" width="180px"></a>
</div>

Performing phylogenetic analysis with whole or core genome sequences maximizes the information used to estimate phylogenies and the resolution of closely related species. Usually sequences are aligned with a reference species or strain. However genome alignment is a process that does not scale well computationally. kSNP4 is a program that identifies SNPs without doing alignments and a reference genome. This permits the inclusion of hundreds of microbial genomes that can be processed in a realistic time scale. Such an alignment free technique comes about through the insight that SNPs can be detected in small odd length chunks of sequence, **kmers**. So you split up the genomes into kmers and compare them. If the kmers are otherwise identical and are long enough not to be random, they can be compared for SNPs between many samples. 

kSNP is mainly used for the analysis of viral and prokaryotic genomes. The input data are genome sequences in FASTA format. It can also annotate the SNPs if you include at least one annotation file.

## Usage

```
MakeKSNP4infile -indir entrez_assemblies/ -outfile entrez_genomes.txt
kSNP4 -core -NJ -k 31 -outdir ksnp_assembled -in assembled_genomes.txt
```

## Speed



## Links

* [kSNP3.0: SNP detection and phylogenetic analysis of genomes without genome alignment or reference genome](https://academic.oup.com/bioinformatics/article/31/17/2877/183216)
* [Building Phylogenetic Trees From Genome Sequences With kSNP4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10640685/)