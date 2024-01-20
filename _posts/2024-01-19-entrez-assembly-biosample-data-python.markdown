---
layout: post
title:  "Fetch assembly and associated biosample data using Entrez tools in Python"
date:   2024-01-19 16:10:00
categories: bioinformatics
tags: [python,genomics,microbiology]
thumbnail: /img/entrez_dna.png
---

## Background

This is a somewhat altered version of code from an [old post](/bioinformatics/assemblies-genbank-python) for downloading assemblies. GenBank provides access to information on all it's assembled genomes via the **assembly** database. You have several options: download assemblies individually or in bulk via the website or the Entrez database search system using command line tools without the website. The later is good for large automated searches/downloads. The BioPython package provides an interface to the entrez tools too. 

* esearch - searches an NCBI database for a query and finds the unique identifiers (UIDs) for all records that match.
* esummary - returns document summaries (DocSums) for a list of input UIDs.

## Method

In this case we wish to get meta data on set of assemblies and also retrieve extra fields from the corresponding entry in the **biosample** database. We use the `BioSampleId` to perform a linked search on the biosample db and add it to our result. The final table will also contain links to download the assembly fasta file. Obviously this can be customised to use any other linked db like bioproject. In this case I wanted the `geo_loc_name` field for the biosample.

First the esearch is performed and the returned ids are used to get an esummary for each in turn. This contains assembly accession and biosample id which we put in a dictionary (called `row`). We then run another esearch using that biosampleid. This returns sample data, some of which must be parsed from the XML using BeautifulSoup. All the attributes are added to the dict. (Some of these attributes will be empty and many you likely don't need.) The dict is added to a list which is used to create a pandas DataFrame at the end. Note that I use the `FtpPath_GenBank` to get the link to the assembly file as the RefSeq one is sometimes not present. At the end there is a short function to download the links in the table.

## Code

You should first install these packages (on Ubuntu based systems):

```
sudo apt install ncbi-entrez-direct
```

```python
import os,sys,glob,re
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import Entrez

def get_assembly_summary(id, db="assembly"):
    """Get esummary for an entrez id"""
    
    esummary_handle = Entrez.esummary(db=db, id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def get_assemblies(term):
    """Download genbank assembly meta data for a given search term.
    Args:
        term: search term, usually organism name        
    """

    #provide your own mail here
    Entrez.email = "A.N.Other@example.com"
    handle = Entrez.esearch(db="assembly", term=term, retmax='5000')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    links = []
    result = []
    for id in tqdm(ids[:4]):
        row = {'id':id}
        #get summary
        rec = get_assembly_summary(id)
        
        asm_summ = rec['DocumentSummarySet']['DocumentSummary'][0]       
        fields = ['AssemblyAccession','BioSampleAccn','BioSampleId','SubmitterOrganization']
        for key in fields:
            row[key] = asm_summ[key]
        row['GenbankAccession'] = asm_summ['Synonym']['Genbank']
        
        #biosample info is a separate request using the BioSampleId
        handle = Entrez.esummary(db="biosample", id=asm_summ['BioSampleId'], report="full")        
        rec2 = Entrez.read(handle)
        sampledata = rec2['DocumentSummarySet']['DocumentSummary'][0]['SampleData']
        #parse xml in sampledata
        from bs4 import BeautifulSoup
        soup = BeautifulSoup(sampledata)
        all_attr = soup.findAll('attribute')        
        for attr in all_attr:
            #print (attr,attr['attribute_name'],attr.text)
            row[attr['attribute_name']] = attr.text
        #get url
        #url = asm_summ['FtpPath_RefSeq']
        url = asm_summ['FtpPath_GenBank']
        #print (url)
        if url != '':           
            label = os.path.basename(url)
            #get the fasta link - change this to get other formats
            link = os.path.join(url,label+'.fna.gz')        
            row['link'] = link
        result.append(row)
    result = pd.DataFrame(result)
    return result

def download_links(df, path):
    for i,r in df.iterrows():
        label = r.AssemblyAccession
        urllib.request.urlretrieve(link, os.path.join(path, f'{label}.fna.gz'))
        
```

We call the method like this: 

```python
res = get_assemblies('Mycoplasmopsis bovis')
```

The resulting table will look like this:

```
  id         AssemblyAccession BioSampleAccn BioSampleId  ..
0 19098541   GCF_032463445.1  SAMN37573123    37573123 ..
1 18700801   GCF_016452265.2  SAMN15246619    15246619 ..
2 18700791   GCF_016452245.2  SAMN15246620    15246620 ..
3 18486841   GCF_016452225.2  SAMN15246621    15246621 ..
```

## Links

* [Biopython: Searching-the-Entrez-databases](https://biopython-tutorial.readthedocs.io/en/latest/notebooks/09%20-%20Accessing%20NCBIs%20Entrez%20databases.html#ESearch:-Searching-the-Entrez-databases)
* [The E-utilities In-Depth: Parameters, Syntax and More](https://www.ncbi.nlm.nih.gov/books/NBK25499/)