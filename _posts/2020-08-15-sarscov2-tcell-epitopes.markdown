---
layout: post
title:  "Predicting cross-reactive T cell epitopes in Sars-CoV-2"
date:   2020-08-15 11:12:00
categories: bioinformatics
tags: [epitope,sarscov2,immunology]
thumbnail: /img/scov2_tcell_mhc2.png
---

## Background

<div style="width: 350px; float:right;">
<a href="https://science.sciencemag.org/content/sci/early/2020/08/04/science.abd3871/F4.large.jpg?width=800&height=600&carousel=1"> <img src="https://science.sciencemag.org/content/sci/early/2020/08/04/science.abd3871/F4.medium.gif" width="300px"></a>
</div>

Eight months after the initial outbreak, puzzles remain about the human immune response to the SARS-CoV-2 virus. Despite it's potential lethality to susceptible individuals, the disease is asymptomatic or mild in most people. However antibody tests have often revealed lower than expected rates of seropositivity in populations where the virus has spread. It is almost certain that other components of the immune system are important in protecting individuals or making them non-susceptible to infection. How do we explain the seeming resistance of children to infection? Robust innate immune responses are one candidate. Another possibility is T cells. These modulate the so-called [cell-mediated immune response](https://www.ecdc.europa.eu/en/covid-19/latest-evidence/immune-responses). SARS-CoV-2 reactive CD4+ T cells have been reported in unexposed individuals, suggesting pre-existing cross-reactive T cell memory in 20-50% of people. It's possible these are immune cells that would have been kept around (as 'memory' cells) after one or more previous exposures to the common cold coronaviruses which circulate widely.

A [paper in Science by Mateus et al.](https://science.sciencemag.org/content/early/2020/08/04/science.abd3871) [1] has identified such cross-reactive CD4+ epitopes by generating 42 short term T cell lines from previously identified epitopes in PBMCs from unexposed donors. Then homologs to these peptides in the HCoVs were tested for a response. These tests were done in unexposed and convalescent COVID19 patients. Responding cells in unexposed donors were predominantly found in the effector memory CD4+ T cell population. Cross reactivity was found in 10/42 of the T cell lines. In three cases, HCoV analogs were better antigens than the SARS-CoV-2 peptide. Though this does not prove immunity, it is an important step towards finding the answer [4].

## The human coronaviruses

Human coronaviruses were first identified in the mid-1960s. The seven coronaviruses that can infect people are below. The first four generally cause common cold symptoms. These have been circulating long enough for humans to develop some immunity, even if initially they may have caused more harmful infections. The other three can cause severe acute respitory disease and are newer. SARS-CoV is no longer in circulation.

| name       | description     |
| :--------- | :-------------- |
|[HCoV-OC43](https://en.wikipedia.org/wiki/Human_coronavirus_OC43)| Betacoronavirus 1 which infects humans and cattle. Four HCoV-OC43 genotypes (A to D), have been identified  |
|[HCoV-HKU1](https://en.wikipedia.org/wiki/Human_coronavirus_HKU1)| Betacoronavirus first discovered in January 2005, but dates back further. Found primarily in young children, the elderly, and immunocompromised patients. |
|[HCoV-NL63](https://en.wikipedia.org/wiki/Human_coronavirus_NL63)| Alphacoronavirus identified in late 2004. Enters the host cell by binding to ACE2. Likely circulated in humans for centuries. |
|[HCoV-229E](https://en.wikipedia.org/wiki/Human_coronavirus_229E)| An Alphacoronavirus. Most frequently codetected with other respiratory viruses. |
|MERS-CoV| Causes Middle East Respiratory Syndrome, or MERS. Extremely acute infection but has not spread widely. |
|SARS-CoV | A beta coronavirus that causes severe acute respiratory syndrome, or SARS. |
|SARS-CoV-2 | A beta coronavirus that causes COVID-19. Uses ACE2 receptor and highly transmissable. |

## Predicting epitopes

Using computational methods it's possible to predict such potential cross-reactive CD4+ eptiopes just using the sequences. The method we use here is as follows:
* Predict MHC-binders in each SCoV2 protein sequence and selected the top scoring candidates. Here we use [epitopepredict](https://github.com/dmnfarrell/epitopepredict) to predict the most 'promiscuous' binders across the 8 most representative human MHC-II alleles. Each protein sequence is split into 15-mer peptides and scored. We can really use any length from 9-20 but 15 matches the lengths used in the experimental study.
* Select the top n scoring peptides in each protein. Since the proteins are different sizes, various strategies can be used. In this case we select the peptides above a 95% percentile threshold for each allele and sort them by the numnber of alleles they are present in. We also limit the total for each protein to prevent a very long protein like ORF1ab from dominating the selection.
* Then we can check how well conserved each peptide is with it's closest homologous sequence in each of the other HCoVs by finding it's closest match in the genome and the percentage identity. Then rank them by conservation.
* We can compare our results to those of Mateus et al. who identified 10/42 SCov2 epitopes cross-reactive to human coronaviruses.

This analysis was done in a Jupyter notebook available [here](https://github.com/dmnfarrell/teaching/blob/master/sarscov2/epitopes.ipynb).

## Results

### Most conserved predicted peptides

Using a limit of 70 peptides per protein we found 282 predicted peptides. The top most conserved are shown below. Out of these, 162 were conserved with >67% identity in at least one HCoV (most commonly with Sars). Note that for a peptide to be cross-reactive it does not necessarily have to share all residues in common with it's homolog. The 9-mer core binding sequence could be conserved with perhaps similar residues at ends. 67% is the arbitrary identity used by Mateus et al. The full list can be found [here](https://github.com/dmnfarrell/teaching/blob/master/sarscov2/scov2_netmhciipan_conserved.csv).

|   gene |         peptide |  pos | alleles | sars | 229E | NL63 | OC43 | HKU1 |
|:-------|:----------------|-----:|--------:|-----:|-----:|-----:|-----:|-----:|
| ORF1ab | VNRFNVAITRAKVGI | 5881 |       8 | 0.85 | 0.85 | 0.85 | 0.85 | 0.85 |
| ORF1ab | YLRKHFSMMILSDDA | 5139 |       8 | 0.85 | 0.85 | 0.85 | 0.85 | 0.85 |
| ORF1ab | NFKSVLYYQNNVFMS | 5172 |       8 | 0.85 | 0.69 | 0.77 | 0.85 | 0.85 |
| ORF1ab | MNLKYAISAKNRART | 4933 |       8 | 0.85 | 0.69 | 0.69 | 0.85 | 0.85 |
| ORF1ab | NLKYAISAKNRARTV | 4934 |       8 | 0.85 | 0.69 | 0.69 | 0.85 | 0.85 |
| ORF1ab | NEFYAYLRKHFSMMI | 5134 |       8 | 0.85 | 0.69 | 0.77 | 0.77 | 0.77 |
|      S | AEVQIDRLITGRLQS |  988 |       8 | 0.85 | 0.69 | 0.69 | 0.77 | 0.77 |
| ORF1ab | ATNYDLSVVNARLRA | 5702 |       8 | 0.85 | 0.69 | 0.69 | 0.69 | 0.69 |
| ORF1ab | YDYYRYNLPTMCDIR | 4844 |       8 | 0.85 | 0.69 | 0.00 | 0.69 | 0.69 |
| ORF1ab | SIKNFKSVLYYQNNV | 5169 |       8 | 0.77 | 0.00 | 0.69 | 0.69 | 0.69 |
|      N | RWYFYYLGTGPEAGL |  106 |       4 | 0.85 | 0.00 | 0.00 | 0.85 | 0.85 |
|      M | RSMWSFNPETNILLN |  106 |       5 | 0.85 | 0.00 | 0.00 | 0.77 | 0.77 |
| ORF1ab | RIMASLVLARKHTTC | 5022 |       8 | 0.85 | 0.00 | 0.00 | 0.77 | 0.69 |
| ORF1ab | WTAFVTNVNASSSEA | 6987 |       8 | 0.85 | 0.00 | 0.00 | 0.69 | 0.77 |
|      S | DRLITGRLQSLQTYV |  993 |       7 | 0.85 | 0.69 | 0.69 | 0.00 | 0.00 |
| ORF1ab | AFVNLKQLPFFYYSD | 6359 |       8 | 0.85 | 0.00 | 0.00 | 0.69 | 0.69 |
|      S | FIEDLLFNKVTLADA |  816 |       6 | 0.85 | 0.00 | 0.00 | 0.69 | 0.69 |
| ORF1ab | GDAVVYRGTTTYKLN | 5529 |       8 | 0.85 | 0.00 | 0.00 | 0.69 | 0.69 |
|      M | SFRLFARTRSMWSFN |   98 |       8 | 0.85 | 0.69 | 0.00 | 0.00 | 0.00 |
|      M | LSYFIASFRLFARTR |   92 |       8 | 0.77 | 0.77 | 0.00 | 0.00 | 0.00 |
| ORF1ab | YFVLTSHTVMPLSAP | 5547 |       8 | 0.85 | 0.00 | 0.69 | 0.00 | 0.00 |
|      S | ALNTLVKQLSSNFGA |  957 |       8 | 0.85 | 0.00 | 0.00 | 0.69 | 0.00 |
|      M | ASFRLFARTRSMWSF |   97 |       8 | 0.85 | 0.69 | 0.00 | 0.00 | 0.00 |

### Hits to experimental data

We finally check our peptides against the 10 epitopes identified by Mateus et al. These were found by testing responses in short term T cell lines generated using a primary screen with SCoV2 peptide pools. See Figure 4 in the [paper](https://science.sciencemag.org/content/early/2020/08/04/science.abd3871). The hits in our 162 peptides are shown below along with the experimental results from Mateus et al. We find a hit in 6/10 cases. Some hits are two peptides in our set overlapping which probably indicates the same core epitope.

| Sequence        | Protein | Start | "+"/tested | SFC   | CD4R-30 | CD4S-31 | hit                                |
|-----------------|---------|-------|------------|-------|---------|---------|------------------------------------|
| PSGTWLTYTGAIKLD | N       | 326   | 1/15       | 1067  | Yes     | No      | GTWLTYTGAIKLDDK                   |
| SFIEDLLFNKVTLAD | S       | 816   | 7/15       | 30487 | No      | Yes     | FIEDLLFNKVTLADA, DLLFNKVTLADAGFI |
| YEQYIKWPWYIWLGF | S       | 1206  | 1/17       | 200   | No      | Yes     | None                               |
| VLKKLKKSLNVAKSE | nsp8    | 3976  | 1/16       | 5660  | Yes     | No      | VVLKKLKKSLNVAKS, EVVLKKLKKSLNVAK |
| KLLKSIAATRGATVV | nsp12   | 4966  | 1/17       | 187   | Yes     | No      | RQFHQKLLKSIAATR                  |
| EFYAYLRKHFSMMIL | nsp12   | 5136  | 2/18       | 787   | Yes     | No      | NEFYAYLRKHFSMMI, YLRKHFSMMILSDDA |
| LMIERFVSLAIDAYP | nsp12   | 5246  | 2/17       | 3870  | Yes     | No      | None                               |
| TSHKLVLSVNPYVCN | nsp13   | 5361  | 1/17       | 160   | Yes     | No      | None                               |
| NVNRFNVAITRAKVG | nsp13   | 5881  | 1/18       | 760   | Yes     | No      | VNRFNVAITRAKVGI                  |

## Summary

Mateus et al. [1] have shown that CD4+ T cell memory in some donors could be a contributing factor to immunity. T cells receptors are much less specific to their epitope than antibodies. HCoV binding antibodies may be cross reactive with SARS-CoV-2 but many of these may not be neutralizing antibodies [5]. T cell epitopes are more readily conserved. This exercise shows that it is possible to computationally predict the same peptides with some degree of success. Note that the authors did also use predicted epitopes for their pooling of non-spike epitopes, so this enhances the chances of overlap to our set.

## References

1. J. Mateus et al., “Selective and cross-reactive SARS-CoV-2 T cell epitopes in unexposed humans,” Science (80-. )., vol. 3871, no. August, p. eabd3871, Aug. 2020.
2. A. Grifoni et al., “A sequence homology and bioinformatic approach can predict candidate targets for immune responses to SARS-CoV-2,” Cell Host Microbe, pp. 1–10, 2020.
3. V. Baruah and S. Bose, “Immunoinformatics-aided identification of T cell and B cell epitopes in the surface glycoprotein of 2019-nCoV,” J. Med. Virol., no. February, pp. 495–500, 2020.
4. Sekine, Takuya, et al. ‘Robust T Cell Immunity in Convalescent Individuals with Asymptomatic or Mild COVID-19’. BioRxiv, June 2020, p. 2020.06.29.174888. www.biorxiv.org, doi:10.1101/2020.06.29.174888.
5. Lv H, Wu NC, Tsang OTY, Yuan M, Perera RAPM, Leung WS, et al.. Cross-reactive Antibody Response between SARS-CoV-2 and SARS-CoV Infections. Cell Rep. 2020; doi: 10.1016/j.celrep.2020.107725.

## links

* [notebook on github](https://github.com/dmnfarrell/teaching/blob/master/sarscov2/epitopes.ipynb)
* [Common coronaviruses](https://www.cdc.gov/coronavirus/types.html)
* [Recently discovered human coronaviruses](https://pubmed.ncbi.nlm.nih.gov/19892230/)
* [Immune responses and immunity to SARS-CoV-2](https://www.ecdc.europa.eu/en/covid-19/latest-evidence/immune-responses)
* [Table S1, Mateus et al.](https://science.sciencemag.org/content/sci/suppl/2020/08/03/science.abd3871.DC1/abd3871-Mateus-SM.pdf)
