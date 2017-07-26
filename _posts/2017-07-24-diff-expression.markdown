---
layout: post
title:  "Differential gene expression plugin"
date:   2017-07-24 18:23:00
categories: dataexplore
tags: [plugins]
---

## About

This plugin, added in version 0.8.1, will be of interest to researchers in biology. This kind of analysis is used for finding diffrences in genes, for example represented by the amount of RNA measured by sequencing. It involves taking the normalised read count data and performing statistical analysis to discover quantitative changes in expression levels between experimental groups or conditions. This can be a somewhat involved process and is often done by specialists. This plugin allows any user to do the analysis assuming their data is prepared correctly.

## Requirements

This plugin requires that you install the R language and some R packages. You will not need to use R directly.

#### Linux

Installation is via the package managers so on Ubuntu:

```
sudo apt install r-base
```

#### Windows/Mac

Go to https://cran.r-project.org/ and download the installers.

Install edgeR and limma packages. These are provided as part of the bioconductor project. You can install from the command line as follows:

```
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("edgeR")
```

## Example: miRNA expression

As an example, consider an experiment using miRNA data from the prefrontal cortex of six human individuals suffering Alzheimers disease and from six control subjects (Lau et al. 2013). Assuming the genes have been counted from the raw files you should have a set of gene expression values for each sample as read counts. These are usually in a text file with a column for each sample (the counts file). We also need a file that has labels for each sample indicating which experimental conditions they correspond to, e.g. a time point, infection state. This allows us to map the column names in the counts file to group by the conditions we want to compare.

These files have already been imported and are available as a dataexplore project file <a href="/other/de_example.dexpl">here</a>. If you open this you will see two sheets. One is the counts file ('mirbase_mirna_counts') and the other are the sample labels ('phenodata'). Notice that the counts file has a column for each sample. To perform the analysis:

* go to the *mirbase_mirna_counts* sheet and open the plugin
* a dialog appears below the table, select *phenodata* as the sample_labels
* select *sampleFile* as sample_col and *diseaseStatus* as factors_col
* finally conditions are entered as a comma separated pair, in this case: *CT,AD*
* then click 'Run DE'

The dialog should look like this:

<div style="width: 350px;">
<a href="/img/de_dialog.png"><img src="/img/de_dialog.png" width="500px"></a>
</div>

The result is a set of genes significantly different between conditions. You can open this as a table and plot the results comparing the sets of samples:

<div style="width: 400px;">
<a href="/img/de_genes_plot.png"><img src="/img/de_genes_plot.png" width="500px"></a>
</div>

Changing the fold change cutoff and reads cutoff are useful filters. You can also make an MD plot and gene cluster heatmap of the filtered (significant) genes.

<div style="width: 400px;">
<a href="/img/de_clustermap.png"><img src="/img/de_clustermap.png" width="500px"></a>
</div>


## Limitations

This plugin currently implements pairwise comparisons only.

## Links

* [Plugins](https://github.com/dmnfarrell/pandastable/wiki/Plugins)
* [Explanation of DE at the EBI](https://www.ebi.ac.uk/training/online/course/functional-genomics-ii-common-technologies-and-data-analysis-methods/differential-gene)
* [Information on example dataset](https://www.ncbi.nlm.nih.gov//geo/query/acc.cgi?acc=GSE48552)

## References

* Lau P et al. (2013) Alteration of the microRNA network during the progression of Alzheimer's disease EMBO molecular medicine 5:1613-1634 doi:10.1002/emmm.201201974
* D. Farrell, “smallrnaseq : short non coding RNA-seq analysis with Python,” Bioarxiv, 2017.
