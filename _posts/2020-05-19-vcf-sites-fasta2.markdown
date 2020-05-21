---
layout: post
title:  "Fasta alignment from a multi sample VCF - a less terrible method"
date:   2020-05-19 15:20:00
categories: bioinformatics
tags: [python,genomics,ngs]
thumbnail: /img/fasta_from_vcf2.png
---

## Background

I [previously posted](/bioinformatics/vcf-sites-fasta) a method to extract all the snps from a vcf with many samples to a fasta file. It turns out to have been pretty rubbish. Use this method instead. This was mainly because I didn't properly read the documentation for `pyvcf` which does this perfectly easily and rapidly. It also shows the disadvantages of using the first solution you come across on biostars or stackoverflow. Though this was meant as a placeholder solution at the time.
I have left the previous post up as it serves as a comparison to this page. That solution did not scale well since it used `pyfaidx` to randomly access the reference fasta file at every site for each sample. This isn't necessary at all since the reference allele is (obviously) stored in the vcf file anyway.

## Imports

```python
import sys,os
import pandas as pd
from Bio import SeqIO
import vcf
```

## Code

This function returns a list of `SeqRecord` objects and a matrix of the sites in the form of a pandas DataFrame. It iterates over all records, each of which has a `samples` object. Then it simply adds the `sample.gt_bases` value to a list, one for each sample. This can be made into a SeqRecord object and then saved to a fasta file. In only return the matrix for my own purposes and it may not be required.

```python
def fasta_alignment_from_vcf(vcf_file):
    """Get snp site alt bases as sequences from all samples in a vcf file"""

    import vcf
    from collections import defaultdict
    vcf_reader = vcf.Reader(open(vcf_file, 'rb'))  
    def default():
        return []
    result = defaultdict(default)
    sites = []
    for record in vcf_reader:
        ref = record.REF
        result['ref'].append(record.REF)
        sites.append(record.POS)
        for sample in record.samples:
            name = sample.sample
            if sample.gt_bases != None:
                result[name].append(sample.gt_bases)
            else:
                result[name].append(record.REF)
    print ('found %s sites' %len(sites))
    recs = []
    for sample in result:
        seq = ''.join(result[sample])
        seqrec = SeqRecord(Seq(seq),id=sample)
        recs.append(seqrec)

    smat = pd.DataFrame(result)
    smat.index = sites
    return recs, smat
```

## Links

* [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)
* [pyVCF](https://pyvcf.readthedocs.io/en/latest/index.html)
