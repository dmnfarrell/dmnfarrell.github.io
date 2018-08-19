---
layout: post
title:  "Create a bacterial GFF from a genbank file for BCFtools/csq"
date:   2018-08-14 13:45:00
categories: bioinformatics
tags: [bcftools/csq]
---

## Consequence calling

Consequence calling is the computational prediction of functional consequences from nucleotide sequence variants such as SNPs and indels. A typical effect is a single amino acid being altered inside a coding sequence. There are a number of programs for making predictions such as [snpEff](http://snpeff.sourceforge.net/) and [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) from ensembl. These two programs are somewhat difficult to get started with. bcftools/csq is relatively easy to use with a single command. The inputs requires are the vcf file obtained from variant calling, a reference sequence and annotation file in GFF3 format.

## GFF format problems

However bcftools/csq has quite specific requirements for the gff file format. For bacterial genomes, gff files can have varying formats due to how they are generated both by ensembl and the NCBI. See [this issue](https://github.com/samtools/bcftools/issues/530) on the bcftools github page for more details.

In short the format required looks like this:

```
NC_002945.4	feature	gene	1	1524	.	+	.	ID=gene:BQ2027_MB0001;Name=dnaA;biotype=protein_coding;
NC_002945.4	feature	mRNA	1	1524	.	+	.	ID=transcript:BQ2027_MB0001;Parent=gene:BQ2027_MB0001;biotype=protein_coding
NC_002945.4	feature	CDS	1	1524	.	+	0	ID=CDS:BQ2027_MB0001;Parent=transcript:BQ2027_MB0001;biotype=protein_coding
```

In attempting to use the [M.Bovis](https://www.ncbi.nlm.nih.gov/genome/161?genome_assembly_id=354758) annotation, I found that the gff file provided on genbank would not work as it is nothing like the format above. The easiest solution was to make a new gff file from the genbank file using a Python script. The code is given below and may be of use to others using non-standard bacterial genomes. Though it is possible this issue will only affect a small number of older annotations. Also the code might need to be tweaked to suit your genbank input. For example the M.bovis file has very long 'note' qualifiers that are probably not in most annotations.

## Running bcftools/csq

With the proper gff file, you can then run the calling using this command:

```bcftools csq -f Mbovis_AF212297.fa -g  Mbovis_csq_format.gff sample.vcf -Ot -o sample.csq.tsv```

## Code

To run this you need the biopython and bcbio-gff packages.

```python
def GFF_bcftools_format(in_file, out_file):
    """Convert a bacterial genbank file from NCBI to a GFF3 format that can be used in bcftools csq.
    see https://github.com/samtools/bcftools/blob/develop/doc/bcftools.txt#L1066-L1098.
    Args:
        in_file: genbank file
        out_file: name of GFF file
    """

    from BCBio import GFF
    in_handle = open(in_file)
    out_handle = open(out_file, "w")
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import FeatureLocation
    from copy import copy

    for record in SeqIO.parse(in_handle, "genbank"):
        #make a copy of the record as we will be changing it during the loop
        new = copy(record)
        new.features = []
        #loop over all features
        for i in range(0,len(record.features)):          
            feat = record.features[i]
            q = feat.qualifiers
            #remove some unecessary qualifiers
            for label in ['note','translation','product','experiment']:
                if label in q:
                    del q[label]
            if(feat.type == "CDS"):
                #use the CDS feature to create the new lines
                tag = q['locus_tag'][0]
                q['ID'] = 'CDS:%s' %tag
                q['Parent'] = 'transcript:%s' %tag
                q['biotype'] = 'protein_coding'

                #create mRNA feature
                m = SeqFeature(feat.location,type='mRNA',strand=feat.strand)
                q = m.qualifiers
                q['ID'] = 'transcript:%s' %tag
                q['Parent'] = 'gene:%s' %tag
                q['biotype'] = 'protein_coding'
                new.features.append(m)

            elif(record.features[i].type == "gene"):
                #edit the gene feature
                q=feat.qualifiers
                q['ID'] = 'gene:%s' %q['locus_tag'][0]
                q['biotype'] = 'protein_coding'
                if 'gene' in q:
                    q['Name'] = q['gene']
            new.features.append(feat)
        #write the new features to a GFF                                      
        GFF.write([new], out_handle)
        return
```

## Links

* [2011 blog post on NCBI's gff files](https://blastedbio.blogspot.com/2011/08/why-are-ncbi-gff3-files-still-broken.html)
* [bcftools csq paper](https://www.ncbi.nlm.nih.gov/pubmed/28205675)
* [man page](http://samtools.github.io/bcftools/bcftools-man.html#csq)
* [github page for bcftools](https://github.com/samtools/bcftools)
