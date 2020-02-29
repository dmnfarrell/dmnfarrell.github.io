---
layout: post
title:  "Updates to a genome annotation on the ENA via Webin-CLI"
date:   2020-02-18 11:30:00
categories: bioinformatics
tags: [genomics]
thumbnail: img/ena-scr.jpg
---

## Background

<div style="width: 350px; float:right;">
 <a href="/img/ena-src.jpg"> <img src="/img/ena-scr.jpg" width="300px"></a>
</div>

The European Nucleotide Archive (ENA) stores and manages genome and sequence data and shares these with DDBJ and NCBI as part of its role in the INSDC. If all these acronyms are already confusing see this [previous post](sequence-databases). The ENA allows submission of new data by different routes depending on the source type. For raw reads they can use the interactive web form, by command line using their **Webin-CLI** program or programmatically. However it seems that for some things you now must use the Webin-CLI. According to the documentation: "Genome and transcriptome assemblies can only be submitted using the Webin-CLI submission interface." NCBI have their own [submission portal](https://submit.ncbi.nlm.nih.gov/) which some might find easier to use. _It is only necessary to submit to one of the INSDC databases since they are synced._

## Updating (older) genome annotations is a bit of a pain

Using the Webin command line tool seems tricky at first but it is at least [documented quite well](https://ena-docs.readthedocs.io/). One particular task that I found difficult though is genome annotation updates. This wouldn't be a very common problem because I would guess that once the vast majority of genomes are submitted, the annotation is not changed. That's because they are nearly all automatically generated. On Genbank when submitting microbial assemblies you can use their Prokaryotic Annotation Pipeline (PGAP) to do the annotation. For the ENA you can use another program like RAST and upload your file. However if there is already an annotation and you want to update it, it's a little more confusing. Some existing reference genomes were annotated in the days before automation and were manually curated. Some still are. I am referring here to mainly bacterial genomes. Updating them used to require sending an e-mail to the ENA. So using this tool is certainly better. Here is how it's done.

## Preparing annotations

The content of this article assumes you are already working with an existing flat file in `embl` or `genbank` format. As I said, these are largely made automatically nowadays. The problem comes when you want to edit them. You can't do it online interactively, hence the need to go through the submission process. Artemis is a Sanger program traditionally used to edit these files. In my case I used Python to do it programmatically as I wished to do batch updates to multiple features.

## Updating a genome annotation

 First you download the Webin-CLI tool from the [github page](https://github.com/enasequence/webin-cli/releases). This requires java and is run using `java -jar`. The next thing is to get an account on Webin [here](https://www.ebi.ac.uk/ena/submit/sra/#home). Then you make a manifest file with the following key details about your genome sequence. An already existing genome will have a STUDY and SAMPLE id. The following example is the one I used for updating the Mbovis AF212297 annotation:

```
STUDY   ERP016893
SAMPLE   ERS1465382
ASSEMBLYNAME     ASM19583v2
ASSEMBLY_TYPE   isolate
COVERAGE   86
PROGRAM   mummer/soapdenovo
PLATFORM     MiSeq/PacBio RSII
FLATFILE	 LT708304_updated_aug19.embl.gz
CHROMOSOME_LIST	chromosome.list.gz
```

For genomes and assemblies you also may have to add a chromosome list file or you'll get an error like this: `ERROR: Invalid number of sequences : 1, Minimum number of sequences for CONTIG is: 2`. The chromosome list file just contains 3 columns with sequence name, chromosome name (specified in the embl file) and location. It has to be gzipped also.

```
LT708304	Mycobacterium_bovis_AF212297	Chromosome
```

You then want to validate your file. Note the context is `genome`:

```bash
java -jar /local/webin-cli-2.2.0.jar -validate -context genome -manifest=submission.manifest -userName yourname -password ******
```

However you may find that your annotation file won't pass the validation step because it had less rigorous checking applied when first submitted. These erors will be put into a report file created in the same folder under `genome/<assembly>/validate/`. Typical reasons are that a `locus_tag` is missing or non-unique. These have to be corrected before you can submit. You might have to edit them in a text file or in Artemis.

Once the file is fixed and you don't get any errors you can submit by replacing the `-validate` option in the command with `-submit`.

## Links

* [ENA submissions](https://ena-docs.readthedocs.io/)
* [NCBI Submission Portal](https://submit.ncbi.nlm.nih.gov/)
* [Artemis](https://www.sanger.ac.uk/science/tools/artemis)
* [Reading and writing genbank/embl files with Python](http://dmnfarrell.github.io/bioinformatics/genbank-python))
