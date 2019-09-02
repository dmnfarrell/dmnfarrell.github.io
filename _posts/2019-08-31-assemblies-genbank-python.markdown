---
layout: post
title:  "Retrieving genome assemblies via Entrez with Python"
date:   2019-08-31 11:40:00
categories: bioinformatics
tags: [genbank,python]
thumbnail: https://upload.wikimedia.org/wikipedia/commons/8/8f/Entrez.svg
---

## Background

GenBank provides access to information on all it's assembled genomes via the Assembly database. It's not that hard to download assemblies individually or in bulk via the website. Though you have to know what you are looking for when it comes to micro-organisms since there can be many assembled genomes for a species. A typical example of a search result would be this page for all [*M. tuberculosis* assemblies](https://www.ncbi.nlm.nih.gov/assembly/?term=mycobacterium+tuberculosis). This page returns 6622 results and it shows a big button allowing you to download the assemblies all at once in multiple formats such as the fasta file or gff and so on.

The above is fine but you may want to do batch downloads without using the web pages. For this the NCBI provides programmatic access via the Entrez query and database system. These utilities can be used via the command line (esearch) but for assemblies I found Python was more flexible. The BioPython package is used to access the Entrez utilities. For the case of assemblies it seems the only way to download the fasta file is to first get the assembly ids and then find the ftp link to the RefSeq or GenBank sequence using `Entrez.esummary`. Then a url request can be used to download the fasta file. The code is presented below and may be adapted to download any of the other formats.

## Code

This code requires pandas and biopython to run.

```python
def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    from Bio import Entrez
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record

def get_assemblies(term, download=True, path='assemblies'):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """

    from Bio import Entrez
    #provide your own mail here
    Entrez.email = "A.N.Other@example.com"
    handle = Entrez.esearch(db="assembly", term=term, retmax='200')
    record = Entrez.read(handle)
    ids = record['IdList']
    print (f'found {len(ids)} ids')
    links = []
    for id in ids:
        #get summary
        summary = get_assembly_summary(id)
        #get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        #get the fasta link - change this to get other formats
        link = os.path.join(url,label+'_genomic.fna.gz')
        print (link)
        links.append(link)
        if download == True:
            #download link
            urllib.request.urlretrieve(link, f'{label}.fna.gz')
    return links
  ```

We can then call the function as follows:

```python
links = get_assemblies("mycobacterium tuberculosis", download=True)
```

## Using the command line

To show how you can accomplish the same task as above using the command line you can try the code below. Unless you are really familiar with using bash commands this might be very hard to read! (This code is from a [biostars post](https://www.biostars.org/p/344959/))

```bash
esearch -db assembly -query 'mycobacterium tuberculosis' \
    | esummary \
    | xtract -pattern DocumentSummary -element FtpPath_GenBank \
    | while read -r line ;
    do
        fname=$(echo $line | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
        wget "$line/$fname" ;
    done
```

## The ncbi-genome-download tool

There are other tools to do this from the command line. These are very well covered in [this blog](https://astrobiomike.github.io/bash/ncbi_eutils). In particular the **ncbi-genome-download** tool is very convenient and flexible to use. However I am not sure if it can be used to accomplish the above for genome assemblies.

## Links

* [GenBank Assembly](https://www.ncbi.nlm.nih.gov/assembly)
* [Parsing GenBank records - BioPython tutorial page](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc55)
* [Fetching metadata of NCBI/GenBank/RefSeq assembly identifiers](https://www.biostars.org/p/345510/)
* [Downloading from NCBI](https://astrobiomike.github.io/bash/ncbi_eutils)
* [ncbi-genome-download tool](https://github.com/kblin/ncbi-genome-download)
