---
layout: post
title:  "Creating a local RefSeq protein blast database"
date:   2018-10-09 15:02:00
categories: bioinformatics
tags: [refseq]
---

## Background

Blasting online sequence databases is a way to retrieve orthologs for a protein of interest. These results can be parsed further. However using the remote blast service can be slow. Using a local version is much faster. The example here is for creating a refseq protein db for bacterial genomes. You can get makeblastdb on Ubuntu by installing the ncbi-blast+ package. For windows see link at bottom.

## Fetching the protein fasta files

see <http://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/> for further info

Download the /refseq/bacteria/assembly_summary.txt file (<ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt>)

List the FTP path (column 20) for the assemblies of interest, in this case those that have "Complete Genome" assembly_level (column 12) and "latest" version_status (column 11). One way to do this would be using the following awk command:

```
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths
```

Append the filename of interest, in this case "*_protein.faa.gz" to the FTP directory names. One way to do this would be using the following awk command:

```
awk 'BEGIN{FS=OFS="/";filesuffix="protein.faa.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
```

Use curl or wget to download the data file for each FTP path in the list, e.g:

```wget -i ftpfilepaths```


Put all the files together and make the database
```
gunzip *.gz
cat *.faa > bac_refseq.fa
#Make the local blast database:
makeblastdb -in bac_refseq.fa -out bacterial_refseq -dbtype prot
```

All commands together
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary_refseq.txt > ftpdirpaths
awk 'BEGIN{FS=OFS="/";filesuffix="protein.faa.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
wget -i ftpfilepaths
gunzip *.gz
cat *.faa > bac_refseq.fa
makeblastdb -in bac_refseq.fa -out bacterial_refseq -dbtype prot
```

### For viral genomes

The same steps as above except we use the following file: <ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt>


## Links

* [About RefSeq](https://www.ncbi.nlm.nih.gov/refseq/about/)
* [NCBI blast+ tools](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
