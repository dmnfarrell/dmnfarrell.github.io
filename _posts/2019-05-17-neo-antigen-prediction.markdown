---
layout: post
title:  "Predicting neoantigens"
date:   2019-05-17 11:04:00
categories: bioinformatics
tags: [epitope,neoantigen,cancer,immunotherapy]
thumbnail: /img/neoantigen_tcr.svg
---

## Background

<figure>
<img src="/img/neoantigen_tcr.svg" alt="T cell recognition" class="scaled">
<figcaption>Tumor cells can present mutated peptides to T cells which respond by killing them if they induce the correct response. If inhibitory signals are produced by PD-L1, the response can be curtailed.</figcaption>
</figure>

A recent development in oncology has been the use of cancer immunotherapy. This leverages methods of enhancing the bodies own immune defences to target cancer cells. Adoptive T cell therapy is one approach. White blood cells (T cells) are taken out of the patients blood and grown in the laboratory. They are augmented in some way and put back into the body where they persist and help the immune system fight the cancer. Such therapy can sometimes induce tumor regression in patients, though response rates are still inadequate. One reason for this is because the tumor microenvironment inhibits the activity of T cells. Drugs called [Checkpoint Blockade Inhibitors](https://en.wikipedia.org/wiki/Checkpoint_inhibitor) are very successful in some cancers in preventing this problem. Personalized cancer vaccines are another approach that boost the number of tumor specific T cells present by administering the right antigens. These then encourage the proliferation of T cells that target cells with that antigen presented at the surface. How then is this antigen selection to be made?

## What are neoantigens?

Cancer antigens are those proteins more likely to be present (or more expressed) in the tumor compared to normal tissue and induce an immune response. Many of these are 'self' antigens because they are also present in normal cells. It is hard to target them because the immune system is likely to ignore them as it does all normal tissue. If this was not the case the immune system would attack your own body as can happen in auto-immunity. Even if it works there can be 'off target' effects on normal tissue also. There are another type of cancer antigen called neoantigens. These are mutated proteins created from somatically mutated genes in the tumor. This creates 'non-self' peptide sequences that are processed by the tumor cell in the normal way and presented to a T cell as shown above. The peptides are called **neoepitopes** if they can be recognised by certain T cells with the right receptor. Neoantigens can safely be targeted by T cells. They are specific to the patient and type of cancer and some cancers make more mutations than others. Clinical response to checkpoint inhibitor treatment does appear to correlate with neoantigen load. If the neoepitope sequences can be predicted efficiently, then cancer vaccines might be improved greatly in the the context of checkpoint inhibitors.

## How to predict them

<div style="width: 470px; float:right;">
<img src="/img/neopeptide_prediction.svg" width="450px">
</div>

Early work in computationally predicting these peptides has already shown some promise. This is achieved by applying whole [exome sequencing](https://en.wikipedia.org/wiki/Exome_sequencing) of matched cancer and normal tissues. Variant calling is then performed and the potential new mutated peptides are predicted. This should include non-synonymous substitutions, indels and frameshift mutations. Usually RNA-seq quantification of tumour gene expression is needed also so non-expressed genes can be filtered out and the right isoforms determined. In theory potential epitopes can then be predicted. The primary determining factors are still not well understood though clearly MHC-binding affinity is important. This is calculated using the standard prediction methods such as netMHC and also requires we HLA type the patient to know their MHC alleles.

If only the binding affinity is considered accurate prediction of actual HLA presented ligands is poor. In fact only a small proportion of binders are found to be neoantigens in practice. Clearly a more nuanced model of ranking epitopes is required. Some studies have noted that the difference between wild type and mutant peptides is significant. This is probably because neantigens too similar to the closest wild type form of the peptide are also subject to immunological tolerance. The opposite rational applies to known viral antigenic peptides in that neoepitopes similar to these might indicate enhanced immunogencity. This combination of multiple factors points to a prediction model that could be fit using known positive and negative peptides using the forementioned features as variables.

Other factors might be included. At the moment making such a model generally applicable is difficult with minimal training data is available. It should also be noted that current methods rely solely on MHC-I binding affinity. This ignores the possible importance of MHC-II binding for presentation to CD4+ T cells which are also involved in tumor responses. Improvements in binding prediction are ongoing. This includes the inclusion of newer data from mass spectrometry identifying peptides eluted from MHC molecules and their length distribution (the 'peptidome'). This implicitly integrates information on the whole peptide processing pathway including MHC binding. Such data has already been utilised in netMHC version 4.

In summary, algorithms effectively able to rank the most likely candidates from the many somatic mutations are still some way off. A more comprehensive modeling of antigen processing is required to improve computation tools. There remains much research to be done to solve this problem.

## Links

* [The problem with neoantigen prediction](https://www.nature.com/articles/nbt.3800)

## References

* Aurisicchio, L. et al. (2018) ‘The perfect personalized cancer therapy: Cancer vaccines against neoantigens’, Journal of Experimental and Clinical Cancer Research. Journal of Experimental & Clinical Cancer Research, 37(1), pp. 1–10. doi: 10.1186/s13046-018-0751-1.
* T. A. Tokuyasu and J.-D. Huang, “A primer on recent developments in cancer immunotherapy, with a focus on neoantigen vaccines,” J. Cancer Metastasis Treat., vol. 4, no. 1, p. 2, 2018.
* Lee CH et al., "Update on Tumor Neoantigens and Their Utility: Why It Is Good to Be Different" Trends Immunol. 2018 Jul;39(7):536-548.
* Abelin, Jennifer G., et al. "Mass spectrometry profiling of HLA-associated peptidomes in mono-allelic cells enables more accurate epitope prediction." Immunity 46.2 (2017): 315-326.
* S. Kim et al., “Neopepsee: Accurate genome-level prediction of neoantigens by harnessing sequence and amino acid immunogenicity information,” Ann. Oncol., vol. 29, no. 4, pp. 1030–1036, 2018.
* Balachandran, Vinod P., et al. "Identification of unique neoantigen qualities in long-term survivors of pancreatic cancer." Nature 551.7681 (2017): 512
* A. M. Bjerregaard et al., “An analysis of natural T cell responses to predicted tumor neoepitopes,” Front. Immunol., vol. 8, no. NOV, pp. 1–9, 2017.
* M. A. Wood et al., “Population-level distribution and putative immunogenicity of cancer neoepitopes.,” BMC Cancer, vol. 18, no. 1, p. 414, 2018.
