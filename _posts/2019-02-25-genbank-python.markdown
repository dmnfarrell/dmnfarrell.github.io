---
layout: post
title:  "Reading and writing genbank/embl files with Python"
date:   2019-02-25 19:00:00
categories: bioinformatics
tags: [annotation,genbank,python]
thumbnail: https://upload.wikimedia.org/wikipedia/commons/thumb/2/24/Genbank100CD.jpg/450px-Genbank100CD.jpg
---

<div style="width: 300px; float:right;">
<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/2/24/Genbank100CD.jpg/450px-Genbank100CD.jpg" width="250px">
</div>

## Background

The GenBank and Embl formats go back to the early days of sequence and genome databases when annotations were first being created. They are a (kind of) human readable format but rather impractical for programmatic manipulation. These formats were designed for annotation and store locations of gene features and often the nucleotide sequence. Though they are not practical for tasks like variant calling, they are still very much used within the main INSDC databases. Thus programming languages with bio libraries like Python have functionality for using them. The Biopython package contains the __SeqIO__ module for parsing and writing these formats which we use below. You could also use the sckit-bio library which I have not tried. Note this method is useful if you want to bulk edit features automatically. For small edits it's much easier to do it manually in a text editor or interactively in Artemis, for example.

## Genbank features

We have recently had the task of updating annotations for protein sequences and saving them back to embl format. Such files contain one or more records with a __feature__ for each coding sequence (or other genetic element). Each feature attribute is called a __qualifier__ e.g. the protein_id (see below). It was useful to be able to write the features to a pandas dataframe, edit this and then rewrite the features using this dataframe to a new embl file.

```
CDS             complement(695671..696906)
                /locus_tag="Y980_RS03120"
                /old_locus_tag="Y980_0597c"
                /inference="COORDINATES: similar to AA
                sequence:RefSeq:NP_215111.1"
                /note="Derived by automated computational analysis using
                gene prediction method: Protein Homology."
                /codon_start=1
                /transl_table=11
                /product="ATP-binding protein"
                /protein_id="WP_003403128.1"
                /translation="MGVVERAIAPSVLAALADTPVVVVNGARQVGKTTLVARLDYPGS
                SEVVSLDDVANRDAARDDPRAFVSRPVDTLVIDEAQLEPGLFRAIKAEVDRDRRPGRF
                LLTGSARLLSAPDMADALVGRVEIIELWPFSQGERAGIADGFVDALFTAPRELIHGSD
                MRRADLVDRIATGGFPDIVARSPSRRRAWFDNYLTTATQSVIREISPIERLAEMPRVL
                RLCAARTGAELNVSALANDLSIPARTTAGYLALLEAAFLIHRVPAWSTNLSRKVIRRP
                KLVVSDSGLACHLLGVTGATLDRPGRPLGPLLETFVANEIRKQLTWSTERPSLWHFRD
                RGGAEVDLVLEHPDGRVCGIEVKATSTPRAEDLRGLRYLAERLDDRFQFGVLLTAAPE
                ATPFGPTLAALPVSTLWAG"
```

## Code

This code requires pandas and biopython to run. We first make a function converting to a dataframe where the features are rows and columns are qualifier values:

```python
def features_to_dataframe(recs, cds=False):
    """Get genome records from a biopython features object into a dataframe
      returns a dataframe with a row for each cds/entry"""

    genome = recs[0]
    #preprocess features
    allfeat = []
    for (item, f) in enumerate(genome.features):
        x = f.__dict__
        q = f.qualifiers
        x.update(q)
        d = {}
        d['start'] = f.location.start
        d['end'] = f.location.end
        d['strand'] = f.location.strand
        for i in featurekeys:
            if i in x:
                if type(x[i]) is list:
                    d[i] = x[i][0]
                else:
                    d[i] = x[i]
        allfeat.append(d)

    import pandas as pd
    df = pd.DataFrame(allfeat,columns=featurekeys)
    df['length'] = df.translation.astype('str').str.len()
    #print (df)
    df = check_tags(df)
    if cds == True:
        df = get_cds(df)
        df['order'] = range(1,len(df)+1)
    #print (df)
    if len(df) == 0:
        print ('ERROR: genbank file return empty data, check that the file contains protein sequences '\
               'in the translation qualifier of each protein feature.' )
    return df

```

Then we can wrap this in a function to easily read in files and return a dataframe:

```python
def embl_to_dataframe(infile, cds=False):
    recs = list(SeqIO.parse(infile,'embl'))
    df = features_to_dataframe(recs, cds)
    return df
```

Say we edit the dataframe table in python (or even in a spreadsheet). We then want to update the feature records and write a new file. A convenient way to handle the features is to scan through them and build up a mapping (a python dictionary) the locus tag to the feature index (from code by Peter Cock). The key used should be unique so locus_tag is best.

```python
def index_genbank_features(gb_record, feature_type, qualifier):
    """Index features by qualifier value for easy access"""

    answer = dict()
    for (index, feature) in enumerate(gb_record.features):
        #print (index, feature)
        if feature.type==feature_type:
            if qualifier in feature.qualifiers:
                values = feature.qualifiers[qualifier]
                if not type(values) is list:
                    values = [values]
                for value in values:
                    if value in answer:
                        print ("WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index))
                    else:
                        answer[value] = index
    return answer
```

This index is then used to find the appropriate feature for updating. We need to use the same key as used in the index, the locus_tag in this case.

```python
def update_features(genome, df, field, key='locus_tag'):
    """Use a dataframe to update a genbank file with new or existing qualifier
    Returns a seqrecord object.
    """

    index = sequtils.index_genbank_features(genome, "CDS","locus_tag")    
    c=0
    for i,r in df.iterrows():        
        t=r[key]
        if t not in index:
            continue
        #print (r)
        #print (t,index[t])
        new = r[field]
        cds = genome.features[index[t]]
        if field not in cds.qualifiers:
            cds.qualifiers[field] = new
            c+=1
        else:
            curr = cds.qualifiers[field][0]
            #print (curr,new)
            if new != curr:
                cds.qualifiers[field] = new
                c+=1
    print ('updated %s features' %c)
    return genome
```

## Usage

Here is how we use all that code together to make new embl files. Here we have edited the product field. The new values will replace the old ones.

```python
from Bio import SeqIO
df = embl_to_dataframe('file.embl','embl')
#edit the dataframe in some way
feats = SeqIO.read('file.embl','embl')
new = update_features(feats, df, 'product')
SeqIO.write(new, 'new.embl', "embl")
```

## Links

* [Dealing with GenBank files in Biopython](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank/)
* [Artemis](https://www.sanger.ac.uk/science/tools/artemis)
* [sckit-bio page on genbank format](http://scikit-bio.org/docs/0.5.2/generated/skbio.io.format.genbank.html)
* [The DDBJ/ENA/GenBank Feature Table Definition](http://www.insdc.org/files/feature_table.html)
