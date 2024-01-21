---
layout: post
title:  "Finding genes in a genome or assembly with Python"
date:   2024-01-21 11:40:00
categories: bioinformatics
tags: [python,genomics]
thumbnail: /img/python_dna.png
---

## Background

<div style="width: 200px; float:right;">
 <a href="/img/python_dna.png"> <img src="/img/python_dna.png" width="180px"></a>
</div>

There are lots of tools for finding specific genes inside genome sequences. The well established technique is to blast the query gene sequence to your own (the target). You generally need to use some threshold of percentage identity and coverage of the sequence to filter results. Here is a method in Python that uses blast to search a genome sequence in a fasta file. The sequence can be anything like full genomes, assembly contigs or a short segment of sequence.

## Method

This method takes the target sequence and makes a local blast database from it using `makeblastdb`. The function `blast_sequences` uses the blastn command to run the search. The code for these functions isn't given here for brevity, they are part of the tools module of the snipgenie package. You can copy them from the source [here](https://github.com/dmnfarrell/snipgenie/blob/master/snipgenie/tools.py). Hits are returned as a dataframe with one hit per row. From this coverage can be calculated. Some other fields are derived here because I was blasting contigs from Spades. Those lines can be removed if not needed. `ident` and `coverage` are used to filter the hits according to how specific and sensitive you want the search. It's determined by your goal i.e. finding paralogs in a complete genome you might use high coverage and moderate identity. **It's often easier to set a threshold low and then see what you get.**

The function assigns the gene name from the query sequence id. So make sure you query sequences have names that are meanginful to you. Something like this for example:

```
>vspL
ATGAAAAAATCAAAGTTTTTACTACTTGGATCAGTAGCATCTTTAGCTTCAATTCCCTTTGTAGCAGCTA
AATGTGGTGAAACCAAAGAAGAAAAGAAACCTGAAGCTGATAAACCAAAGCTTAGCGAAACATTAAAATC
TATTACTGGTAATGATTTAGGAAAAGTACAAGTTGCTGA
>vspK
ATGAAAAAATCAAAGTTTTTACTACTTGGATCAGTAGCTTCATTAGTTTCAATTCCCTTTGTAGCAGCTA
AATGTGGTGAGACCAAAGAAGAAAAGAAACCTGAGCCCGACAAAAATCCAGGTGGAGATAAAAACCCTGG
AGGAGAAAAGA
```

## Code

First install the ncbi blast command line tool on your system. You can so it like this on Ubuntu for example:

```bash
sudo apt install ncbi-blast+
```

```python
import pandas as pd
#you don't need this if you just copy the required tools functions from the module
from snipgenie import tools

def find_genes(target, query, ident=90, coverage=75, duplicates=False, threads=2, **kwds):
    """Find ref genes by blasting the target sequences"""

    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import SeqIO

    queryseqs = list(SeqIO.parse(query,'fasta'))
    print ('blasting %s sequences' %len(queryseqs))
    bl = blast_sequences(target, queryseqs, maxseqs=1, evalue=1e-4,
                               cmd='blastn', show_cmd=True, threads=int(threads))

    #print (bl.iloc[0])
    bl['qlength'] = bl.qseq.str.len()
    bl['coverage'] = bl.length/bl.qlength*100
    bl = bl[bl.coverage>coverage]
    bl = bl[bl.pident>ident]
    #used when we are blasting to contigs from spades. remove if not needed.
    bl['filename'] = bl.sseqid.apply(lambda x: x.split('~')[0],1)
    bl['id'] = bl.filename.apply(lambda x: os.path.basename(x),1)
    bl['contig'] = bl.sseqid.apply(lambda x: x.split('_')[1],1)
    #try to extract gene name from query sequence id, otherwise use it directly
    try:
        bl['gene'] = bl['qseqid'].apply(lambda x: x.split('~~~')[1],1)
    except:
        bl['gene'] = bl.qseqid

    #remove exact and close duplicates
    print (len(bl))
    bl = bl.sort_values(['bitscore'], ascending=False).drop_duplicates(['contig','sstart','send'])
    print (len(bl))
    #this is also optional
    if duplicates == False:
        dist = 20
        x=bl.sort_values(by=["contig","sstart"],ascending=False)
        #print (x[:15][x.columns[:5]])
        unique = x.sstart.diff().abs().fillna(dist)
        bl = bl[unique>=dist]
    cols = ['gene','id','qseqid','pident','coverage','sstart','send','contig','filename','bitscore']
    #print (bl)
    bl = bl[cols]
    return bl
```

The code is called like this. Where `querygenes` contains the sequences you want to search for. `bl` is a DataFrame.

```python
target = 'mygenes.fa'
make_blast_database(target,dbtype='nucl')
bl = find_genes(target,'querygenes.fa',ident=95,coverage=95)
```

The results will look like this:

```
             gene         qseqid   pident  coverage  sstart  send contig  bitscore sample
66  MBOVPG45_0817  MBOVPG45_0817  100.000     100.0    1072   101     67    1754.0    221
62  MBOVPG45_0815  MBOVPG45_0815  100.000     100.0    2406  1534     47    1575.0    221
7            vspK           vspK   99.739     100.0    4425  5189     47    1371.0    221
63           vspJ           vspJ  100.000     100.0     853   128     70    1310.0    221
```

## Links

* [snipgenie source](https://github.com/dmnfarrell/snipgenie/blob/master/snipgenie/tools.py)