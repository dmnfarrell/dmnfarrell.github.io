---
layout: post
title:  "Find PFAM domains in protein sequences with Python"
date:   2020-11-10 12:06:00
categories: bioinformatics
tags: [annotation,mtb]
thumbnail: /img/pfam_scr1.png
---

## Background

We have previously updated the genome annotation for the M.bovis AF2122/97 genome, [detailed here](https://www.microbiologyresearch.org/content/journal/acmi/10.1099/acmi.0.000129), to add additional functional annotation. This used data from Doerks et al. and other sources to add product information to hypothetical proteins. There are still approximately 467 hypotheticals remaining in the genome. Some of these are still identifiable using domain information from PFAM or InterPro. Here Python is used to fetch PFAM annotation for these proteins. You can also use a batch request on the PFAM site to achieve the same goal.

## Code

```python
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from snipgenie import tools
from prody import *
```

First we take the latest version of the genome and load it into a pandas DataFrame for easier processing. The method `genbank_to_dataframe` is part of the *snipgenie* package. But you can just use the method on it's own without the library by copying it from [here](https://github.com/dmnfarrell/snipgenie/blob/master/snipgenie/tools.py). It uses BioPython to parse the genbank file.

```python
#put the genome annotation into a dataframe (one row per feature)
g = tools.genbank_to_dataframe('LT708304_latest.gb')
#remove 'gene' features
g = g[g.feat_type!='gene']

#find the hypotheticals in the genome using the product field
labels=['hypothetical protein','conserved protein','conserved hypothetical',
        'unknown protein','uncharacterized','hypothetical alanine']
def find_unknown(x):
    for l in labels:
        if l in str(x).lower():           
            return True
    return False

g['hypothetical'] = g['product'].apply(find_unknown,1)
#filter
hypo = g[g.hypothetical==True]
```

We can then use the following function to fetch the pfam domain for each sequence, if available. The function uses the [prody](http://prody.csb.pitt.edu/) library to make the call to PFAM which returns a dictionary. Each dict is a domain and there can be several in one protein. The dict is parsed and added to the list `new`, which is then converted to a DataFrame at the end. Note that this function can take a previous table of results so that you can

```python
def find_pfam_domains(g, res=None):
    """Find pfam domains"""

    new=[]
    for i,r in g.iterrows():
        if res!=None:
            if r.locus_tag in list(res.locus_tag):
                continue
        print (r.locus_tag)
        try:
            d = searchPfam(r.translation)
        except Exception as e:
            print (e)
            continue
        found = []
        for k in d:
            found = d[k]       
            found['locus_tag'] = r.locus_tag
            found['length'] = len(r.translation)
            found['product'] = r['product']
            new.append(found)    
    df=pd.DataFrame(new)
    if len(df)>0 and g is not None:
        locs = df['locations'].apply(pd.Series)
        df=pd.concat([df, locs], axis=1)
        df=df.drop(columns=['locations'])
    df = pd.concat([res,df])
    return df
```

We can then run as follows:

```python
res = find_pfam_domains(hypo)
res.to_csv('hypothetical_pfam.csv',index=False)

#get only the non DUF domains
known = res[~res.id.str.contains('DUF')]
#we can merge with the original dataframe if we want
cols=['locus_tag','translation','protein_id','gene']
known = known.merge(g[cols],on='locus_tag')
#get only the non DUF domains
```

## Result

The search produces domains for 188 of 467 of the hypotheticals. Some of these are marked as DUF, meaning domains of unknown function. That leaves 104 with named domains. The results are returned as a dataframe so we can just look at the table:

|      locus_tag |                        product |            id |  accession | start | end | length |
|---------------:|-------------------------------:|--------------:|-----------:|------:|----:|-------:|
|  BQ2027_MB0028 | conserved hypothetical protein | T7SS_ESX_EspC |  PF10824.8 |     1 |  99 |    105 |
| BQ2027_MB0037c |              conserved protein |       MDMPI_N |  PF11716.8 |    11 | 147 |    257 |
| BQ2027_MB0037c |              conserved protein |  Wyosine_form | PF08608.12 |   183 | 234 |    257 |
| BQ2027_MB0037c |              conserved protein |        DinB_2 |  PF12867.7 |    10 | 165 |    257 |
|  BQ2027_MB0060 |           HYPOTHETICAL PROTEIN |          DarT |  PF14487.6 |    29 | 230 |    230 |
| BQ2027_MB0081c |           HYPOTHETICAL PROTEIN |        AbiEii | PF08843.11 |     9 | 196 |    197 |

Looking at the genbank file we can see that many of these proteins already have `/db_xref` links to the relevant InterPro family already. This information was added by GenBank subsequent to the original annotation but the prodict fields were never changed. For example the first row in the table above, MB0028, is in fact an EspC like protein:

```
CDS             31173..31490
               /locus_tag="BQ2027_MB0028"
               /codon_start=1
               /product="conserved hypothetical protein"
               /transl_table=11
               /note="Mb0028, -, len: 105 aa. Equivalent to Rv0027, len:
               105 aa, from Mycobacterium tuberculosis strain H37Rv,
               (100.0% identity in 105 aa overlap). Hypothetical unknown
               protein. Protein product from Mb0028 detected using SWATH
               mass spectrometry. Mb0028 found to be expressed during
               exponential growth in Sauton's minimal media by
               RNA-sequencing."
               /db_xref="GOA:P64668"
               /db_xref="InterPro:IPR022536"
               /db_xref="UniProtKB/Swiss-Prot:P64668"
               /translation="MTDRIHVQPAHLRQAAAHHQQTADYLRTVPSSHDAIRESLDSLGP
               IFSELRDTGRELLELRKQCYQQQADNHADIAQNLRTSAAMWEQHERAASRSLGNIIDGS
               R"
```

We could now use this table to add product information to these 104 proteins.

## Links

* [M. bovis genome on GenBank: LT708304.1](https://www.ncbi.nlm.nih.gov/nuccore/LT708304.1)
* [M. bovis annotation github page](https://github.com/dmnfarrell/gordon-group/tree/master/mbovis_annotation)
* [PFAM](https://pfam.xfam.org/)
