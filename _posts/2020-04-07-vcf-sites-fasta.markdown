---
layout: post
title:  "Create a fasta alignment from a multi sample VCF"
date:   2020-04-07 10:20:00
categories: bioinformatics
tags: [python,genomics,ngs]
thumbnail: /img/fasta_from_vcf.svg
---

**UPDATE**: This method is not recommended, use the function in the [follow up post](/bioinformatics/vcf-sites-fasta2) instead.

## Background

<div style="width: 400px; float: right;">
 <a href="/img/fasta_from_vcf.svg"> <img src="/img/fasta_from_vcf.svg" width="350px"></a>
</div>

A typical endpoint of microbial whole genome sequencing analysis is to construct a MSA (multiple sequence alignment) of the variable sites, most commonly the SNVs (ignoring indels). This is then used to infer a phylogeny. The input is a vcf with all sample sites. To produce a multi sample vcf, you can either call the variants for each sample merge all the single vcfs together or call all samples at once. (Merging many vcfs appears to be memory intensive so I would prefer the latter). The resulting file is then converted to a fasta alignment. There are scripts available to do this but I was not able to find one written in Python that I could easily use. The following details a Python method to do this.

Note: vcf files should always be filtered to get high quality sites for this to produce a reliable phylogeny. False negatives will give erroneous results. This is not detailed here.

Requires these Python packages: pyvcf, pyfaidx, biopython and pandas. The programs tabix (htslib) is used to index the vcf.

## Imports

```python
import sys,os
import pandas as pd
from Bio import SeqIO
from pyfaidx import Fasta
from pyfaidx import FastaVariant
import vcf
```

## Code

This particular function returns a list of `SeqRecord` objects and a matrix of the sites in the form of a pandas DataFrame. It iterates over each sample once to produce a set of all unique sites in all samples, using the `FastaVariant` object from pyfadix. The nucleotide at each sites is recorded for the reference genome. In a second loop over each sample we extract the variant nucleotide at each of the sites and store these in a sequence. The results is a list of sequences per sample. This can then be saved to a fasta file and used as an MSA for input into RAXml or other tree constuction programs. This code has not been fully compared to other algorithms of this kind so use it at your discretion.

```python
def fasta_alignment_from_vcf(vcf_file, ref):
    """
    Get a fasta alignment for all snp sites in a multi
    sample vcf file, including the reference sequence.
    """

    #index vcf
    cmd = 'tabix -p vcf -f {i}'.format(i=vcf_file)
    tmp = subprocess.check_output(cmd,shell=True)
    #get samples from vcf
    vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
    samples = vcf_reader.samples
    print ('%s samples' %len(samples))
    result = []

    #reference sequence
    reference = Fasta(ref)
    chrom = list(reference.keys())[0]

    #get the set of all sites first
    sites=[]
    for sample in samples:      
        variant = FastaVariant(ref, vcf_file,
                                 sample=sample, het=True, hom=True)
        pos = list(variant[chrom].variant_sites)
        sites.extend(pos)

    sites = sorted(set(sites))
    print ('using %s sites' %len(sites))
    #get reference sequence for site positions
    refseq=[]
    for p in sites:
        refseq.append(reference[chrom][p-1].seq)
    refseq = ''.join(refseq)
    #seqrecord for reference
    refrec = SeqRecord(Seq(refseq),id='ref')
    result.append(refrec)

    sites_matrix = {}
    #iterate over variants in each sample
    for sample in samples:        
        seq=[]
        variant = FastaVariant(ref, vcf_file,
                                 sample=sample, het=True, hom=True)     
        for p in sites:        
            rec = variant[chrom][p-1:p]    
            seq.append(rec.seq)
        seq = ''.join(seq)
        #make seqrecord for the samples sites  
        seqrec = SeqRecord(Seq(seq),id=sample)
        result.append(seqrec)
        sites_matrix[sample] = list(seqrec)
    df = pd.DataFrame(sites_matrix)
    df.index=sites    
    return result, df
```

## Links

* [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)
* [pyfaidx](https://pythonhosted.org/pyfaidx/)
* [Processing and Merging SNVs from multiple isolates to get an MSA](https://mtbgenomicsworkshop.readthedocs.io/en/latest/material/day4/merge_snvs_from_vcf.html)
* [biostars thread](https://www.biostars.org/p/94573/)
