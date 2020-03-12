---
layout: post
title:  "Model of the SARS-CoV-2 spike protein in Blender"
date:   2020-03-02 11:30:00
categories: bioinformatics
tags: [genomics,SARS-CoV-2,blender]
thumbnail: /img/spike_thumb.png
---

## Background

SARS-CoV-2 makes use of a glycosylated spike (S) protein to gain entry into host cells. This is a trimeric heterodimer that contains two functional domains: a receptor binding domain (RBD), and a second domain which contains sequences that mediate fusion of the viral and cell membranes. One of the three heterodimers binds the cellular receptor angiotensin-converting enzyme 2 (ACE2). In order to engage a host-cell receptor, the receptor-binding domain of S1 (blue) undergoes hinge-like conformational movements that transiently hide or expose the determinants of receptor binding (this is called the  "up" conformation). Binding takes place and this destabilizes the trimer, resulting in shedding of the S1 subunit, allowing fusion to take place. Due to the indispensable function of the S protein, it represents a target for antibody-mediated neutralization (Wrapp et al.).

This model was created by SWISS-model from template structures, mainly [6ACG](https://www.rcsb.org/structure/6acg). It was then loaded into Pymol, displayed as surface representation and each domain exported as a .wrl image. These can then be imported into a scene in Blender 2.8 as detailed in a [previous article](/bioinformatics/proteins-blender). The image shows the three heterodimers with the blue colored one contacting the ACE2 protein in red. The ACE2 receptor is attached to the human cell and has a transmembrane component that is not shown.

<div style="width: auto; float:left;">
 <a href="/img/sarscov2_spike_blender.png"> <img class="scaled" src="/img/sarscov2_spike_blender.png"></a>
</div>


## References

* Wan, Y., Shang, J., Graham, R., Baric, R. S. & Li, F. Receptor recognition by novel coronavirus from Wuhan: An analysis based on decade-long structural studies of SARS. J. Virol. (2020) doi:10.1128/JVI.00127-20.
* Hoffmann, M. et al., The novel coronavirus 2019 (2019-nCoV) uses the SARS-coronavirus receptor ACE2 and the cellular protease TMPRSS2 for entry into target cells. bioRxiv 2020.01.31.929042 (2020) doi:10.1101/2020.01.31.929042.
* Song W, Gui M, Wang X, Xiang Y (2018) Cryo-EM structure of the SARS coronavirus spike glycoprotein in complex with its host cell receptor ACE2. PLoS Pathog 14(8): e1007236. https://doi.org/10.1371/journal.ppat.1007236
* Wrapp D, Wang N, Corbett KS, et al. Cryo-EM structure of the 2019-nCoV spike in the prefusion conformation [published online ahead of print, 2020 Feb 19]. Science. 2020;eabb2507. doi:10.1126/science.abb2507

## Links

* [SWISS-MODEL SARS-CoV-2](https://swissmodel.expasy.org/repository/species/2697049)
