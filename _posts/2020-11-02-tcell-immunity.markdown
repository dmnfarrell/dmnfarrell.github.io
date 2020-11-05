---
layout: post
title:  "Covid-19 and T cell immunity"
date:   2020-11-02 11:06:00
categories: general
tags: [covid19]
thumbnail: /img/dendritic_cell.png
---

## Background

<div style="width: 250px; float:right;">
<img src="/img/sarscov2.png" width="230px">
</div>

Though spread of SARS-Cov-2 has been widespread in some locations, antibody tests have often revealed lower than expected rates of seropositivity (presence of virus specific antibodies) in these populations. It is probable that many people do not generate antibodies against the virus despite being exposed and avoiding serious disease. In children this could be because the virus simply cannot enter host cells efficiently due to lack of available ACE2 receptors. In some adults or children it may be that their mucosal (IgA) antibodies are sufficient to neutralize the virus and these aren't typically tested for. Even in those who do develop an antibody response, levels in serum will decline over a period of months, depending on the antibody type. This is to be expected, but longer term immunity can still be stored in the form of memory T and B cells. Such immunity can last for many years. This ability to retain a memory response will also be relevant for effective vaccines.

For adults, levels of resistance to disease might also be due to the presence of T cells that have been developed from a previous coronavirus infection and stored in the form of memory cells. There are [four circulating species](https://www.cdc.gov/coronavirus/types.html) of human coronavirus (229E, NL63, OC43 and HKU1) that normally cause common colds. They are similar enough to SARS-Cov-2 that some of the T cells developed in response to them can also be reactive to Sars-Cov-2 antigen.  Such 'cross reactive' T cells have been shown to be present in the blood of perhaps 20-50% individuals not previously exposed to SARS-Cov-2. Whether these are protective or not is not yet fully known.

So what are T cells? T and B cells are classes of white blood cells that are central to the so-called cell-mediated immune response. T cells act to recognise infected cells and either kill them or assist in doing so. B cells produce antibodies that kill the virus before it can spread to cells.

## Cell mediated immunity and viruses

Viruses can enter the body by many routes. (SARS-Cov-2 appears to enter epithelial cells in the nasal passage and throat by binding to a receptor called ACE2 on the surface of the cells). At this stage of infection, innate immune mechanisms are initiated in response to the binding of pathogens by [pattern-recognition receptors (PRRs)](https://www.immunology.org/public-information/bitesized-immunology/receptors-and-molecules/pattern-recognition-receptor-prrs). These recognise generic molecules specific to viruses or bacteria but not normally present in the body. This promotes the migration of certain cells that express these receptors (such as macrophages) to secondary lymphoid organs. Here, they present virus-derived peptides on [MHC](https://immunobites.com/2018/07/23/what-is-mhc-and-why-does-it-matter/) class II molecules to 'naive' helper (CD4+) T cells. If these have the right virus/peptide specific receptor they will be activated. The activated CD4+ T cells undergo cell division and differentiate into different types. Many become 'effector' cells which migrate to the infected tissue and make molecules called cytokines that stimulate the other cells. The image below shows some of the essential components in this 'cell mediated' response. Many details are excluded but the essential role of CD4+ T Cells is clear. They help other cells specific to the virus to go into action. CD8+ T cells are another type, important for directly attacking and killing virus-infected cells, whereas helper T cells are crucial to prime both CD8+ T cells and B cells. During the infection some T cells differentiate into the memory type. Clones of memory T cells can persist for decades.

<div style="width: auto; float:center;">
 <a href="/img/tcells_infection.png"> <img class="scaled" src="/img/tcells_infection.png"></a>
 <p class="caption">A simplified picture of T Cell involvement in the viral immune response.</p>
</div>
<br>

## Memory T Cells

All T cells have their own unique receptor (TCR) so they will be specific only to a subset of antigenic peptides. As mentioned above, when a naive T cell finds it antigen, it makes many clones. Some of the effector cells will become memory cells. The exact process is still not fully understood. Dendritic cells are a type of cell that can take up a pathogen and can present antigen to T CD4 cells by MHC class II, shown below and above. If it's a memory type CD4 cell it will quickly expand to produce many more cells and boost the immune reponse of the corresponding antigen/virus specific B cells which can make antibodies. (To do this the B cell also presents antigen to the T cells). This process is much faster than the initial response, hence you will be wholly or partially immune to infection. The cross reactivity comes in because these memory cells could have been formed during another infection with a similar virus.


<div class ="image-gallery">
<div class="box">
 <a href="/img/dendritic_cell.png"> <img class="scaled" src="/img/dendritic_cell.png"></a>
 <p class="caption">Dendritic cells can activate a memory T cell that recognises the right virus derived peptide. This memory cell could have been formed during a previous infection with the same virus or be cross reactive from another related virus.</p>
 </div>
 <div class="box">
 <a href="/img/tcr_mhc2.png"> <img class="scaled" src="/img/tcr_mhc2.png"></a>
  <p class="caption">Close up of the TCR/peptide-MHC interface. The TCR has to 'fit' the shape of the p-MHC for a stable binding interaction to occur. This is an MHC-II so it also has to bind the CD4 receptor on the T cell surface. </p>
 </div>
</div>

## T cell cross-reactivity

How is it possible that T cells specific to antigen in one virus can also recognise one from a different species? This is because of the way T cells are activated. Cells present bits of viral protein as peptides in the MHC molecule to the T cell receptor. Strong binding initates the response (and the peptide is called an epitope). If the peptide is very close in sequence, say one or two amino acids, it's structural differences could be minor. This is called a [conservative replacement](https://en.wikipedia.org/wiki/Conservative_replacement). Also the TCR has flexible loops that can accomodate variations in the peptide, as shown below. Finally, most TCRs focus on two to four upward-facing peptide residues and multiple amino acid substitutions at other positions can be tolorated without changing specificity. Therefore the peptide-MHC may still bind to the TCR of a memory cell previously generated.

<div style="width: auto; float:center;">
 <a href="/img/tcr_cross_reactive.png"> <img class="scaled" src="/img/tcr_cross_reactive.png"></a>
 <p class="caption">Two related peptides can be bound to the same TCR. The loops are flexible and can accomodate certain changes to the peptide. </p>
</div>

We can see this illustrated in part of the sequence alignment of the Spike protein of four HCoVs and Sarscov2. The highlighted part, `SFIEDLLFNKVTLADAG`, is a partially conserved peptide. Some of the differences are conservative changes (e.g. A->V) and may not alter the peptide binding significantly. In fact some HCoV homologs of this peptide were found by Mateus et al. [2] to elicit a response in Sars-Cov-2 specific T cells.

<pre><code>
scov2           NFSQILPDPSK---------PSKR<b>SFIEDLLFNKVTLADAG</b>-FIKQYGDCLGDIAARDLI
sars            NFSQILPDPLK---------PTKR<b>SFIEDLLFNKVTLADAG</b>-FMKQYGECLGDINARDLI
OC43            NVDDINFSPVLGCLGSECSKASSR<b>SAIEDLLFDKVKLSDVG</b>-FVEAYNNCTGGAEIRDLI
HKU1            DVDNINFKSLVGCLGPHCG-SSSR<b>SFFEDLLFDKVKLSDVG</b>-FVEAYNNCTGGSEIRDLL
229E            NLSSVIPSLPTSG-----SRVAGR<b>SAIEDILFSKLVTSGLG</b>TVDADYKKCTKGLSIADLA
NL63            NLSSVLPQRNIRS-----SRIAGR<b>SALEDLLFSKVVTSGLG</b>TVDVDYKSCTKGLSIADLA
                :...:  .             : ** :**:**.*:  :. * .   * .*  .    **
</code></pre>

## The future

T cells from patients who became infected with SARS-CoV in 2003 have been shown to cross-react with the SARS-CoV-2 virus 17 years later. Though the extent to which they can provide protection is not known. Clearly many people are not protected from infection by this means since the virus has spread widely. However many remain asymptomatic and pre-existing cellular immunity may be a significant factor in these cases. This is the next step in such research. There has also been much talk of reinfection, though it appears to be uncommon even though antibodies wane. More specific studies are needed in how pre-existing immunity to HCoVs may impact upon the outcome of infection. This will be important for estimating the actual thresholds for population wide immunity.

Note: Immunology can be an extremely complex and confusing subject with a lot of arcane terminology. The immune system iself is a sophisticated network of feedbacks and redundancies that often defies easy description. (This may be one reason why people are so confused about what they see in the media). An article like this will gloss over many details and is rather simplified. I am not an immunologist. For further details try the links below.

## Links

* [Immune responses and immunity to SARS-CoV-2](https://www.ecdc.europa.eu/en/covid-19/latest-evidence/immune-responses)
* [What is the role of T cells in COVID-19 infection?](https://www.cebm.net/covid-19/what-is-the-role-of-t-cells-in-covid-19-infection-why-immunity-is-about-more-than-antibodies/)
* [Not just antibodies: B cells and T cells mediate immunity to COVID-19](https://www.nature.com/articles/s41577-020-00436-4)
* [Why must T cells be cross-reactive?](https://www.nature.com/articles/nri3279)
* [T-Cells from recovered COVID-19 patients show promise to protect vulnerable patients from infection](https://www.prnewswire.com/news-releases/t-cells-from-recovered-covid-19-patients-show-promise-to-protect-vulnerable-patients-from-infection-301159600.html)
* [T-cell Covid immunity 'present in adults six months after first infection'](https://www.theguardian.com/world/2020/nov/02/t-cell-covid-immunity-present-in-adults-six-months-after-first-infection)

## References

1. Swain SL, McKinstry KK, Strutt TM. Expanding roles for CD4⁺ T cells in immunity to viruses. Nat Rev Immunol. 2012;12(2):136-148. Published 2012 Jan 20. doi:10.1038/nri3152
2. J. Mateus et al., “Selective and cross-reactive SARS-CoV-2 T cell epitopes in unexposed humans,” Science (80-. )., vol. 3871, no. August, p. eabd3871, Aug. 2020.
3. Braun J, et al. SARS-CoV-2-reactive T cells in healthy donors and patients with COVID-19. Nature. 2020 Jul 29. doi: 10.1038/s41586-020-2598-9. Epub ahead of print. PMID: 32726801.
4. Grifoni, A. et al. Targets of T cell responses to SARS-CoV-2 coronavirus in humans with COVID-19 disease and unexposed individuals. Cell 181, 1489–1501 (2020).
