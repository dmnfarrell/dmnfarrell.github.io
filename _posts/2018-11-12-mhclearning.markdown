---
layout: post
title:  "Create an MHC-Class I binding predictor in Python"
date:   2018-11-12 15:02:00
categories: bioinformatics
tags: [mhc]
---

## Background

<div style="width: 400px; float:right;">
<img src="/img/mhc-molecule.png" width="400px">
</div>

Peptide binding to MHC molecules is the key selection step in the Antigen-presentation pathway. This is essential for T cell immune responses. The 'epitope' is the peptide-MHC combination shown in the image at right. Key residues in the MHC contact the peptide and these differ between alleles. The prediction of peptide binding to MHC molecules has been much studied. The problem is simpler for class-I molecules since the binding peptide length is less variable (usually 8-11 but commonly 9). Typically binding predictors are based on training models with experimental binding affinity measurements with known peptide sequences. This data is available from the IEDB for many human alleles. New peptides can then be predicted based on their position specific similarity to the training data.

This requires encoding the peptide amino acid sequence numerically in a manner that captures the properties important for binding. Many possible encodings have been suggested and three are illustrated below.

### Peptide encoding

Several encoding techniques have been proposed for representing sequence of amino acids in multidimensional metric spaces. In particular in this work we are interested in a simple encoding that is suited to be coupled with a machine learning algorithm. We will use pandas dataframes to construct the encoding, though probably not the most optimal for speed, it is convenient. First we import the required packages.

```python
import os, sys, math
import numpy as np
import pandas as pd
%matplotlib inline
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("notebook", font_scale=1.4)
import epitopepredict as ep
```

## Encode peptides

To encode a peptide a few schemes are illustrated here. None of these methods take into account the interdependence of the amino acids in terms of their relative positions.

### One hot encoding
The first and simplest is a so-called 'one-hot encoding' of the amino acids producing a 20-column vector for each position that only contains a 1 where the letter corresponds to that amino acid.

The code below uses a pandas dataframe to construct the new features. The flatten command at the end re-arranges the 2-D matrix into a 1-D format so it can be used with the regression models. This applies to the other encoding methods also. The `show_matrix` method draws the 2D matrix as a table.


```python
codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def show_matrix(m):
    #display a matrix
    cm = sns.light_palette("seagreen", as_cmap=True)
    display(m.style.background_gradient(cmap=cm))

def one_hot_encode(seq):
    o = list(set(codes) - set(seq))
    s = pd.DataFrame(list(seq))    
    x = pd.DataFrame(np.zeros((len(seq),len(o)),dtype=int),columns=o)    
    a = s[0].str.get_dummies(sep=',')
    a = a.join(x)
    a = a.sort_index(axis=1)
    show_matrix(a)
    e = a.values.flatten()
    return e

pep='ALDFEQEMT'
e=one_hot_encode(pep)
```

<style  type="text/css" >
    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col0 {
            background-color:  #2e8b57;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col1 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col2 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col3 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col4 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col5 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col6 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col7 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col8 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col9 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col10 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col11 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col12 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col13 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col14 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col15 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col16 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col17 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col18 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col19 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col0 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col1 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col2 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col3 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col4 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col5 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col6 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col7 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col8 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col9 {
            background-color:  #2e8b57;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col10 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col11 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col12 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col13 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col14 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col15 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col16 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col17 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col18 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col19 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col0 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col1 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col2 {
            background-color:  #2e8b57;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col3 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col4 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col5 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col6 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col7 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col8 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col9 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col10 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col11 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col12 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col13 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col14 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col15 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col16 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col17 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col18 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col19 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col0 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col1 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col2 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col3 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col4 {
            background-color:  #2e8b57;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col5 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col6 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col7 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col8 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col9 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col10 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col11 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col12 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col13 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col14 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col15 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col16 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col17 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col18 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col19 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col0 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col1 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col2 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col3 {
            background-color:  #2e8b57;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col4 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col5 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col6 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col7 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col8 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col9 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col10 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col11 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col12 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col13 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col14 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col15 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col16 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col17 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col18 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col19 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col0 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col1 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col2 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col3 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col4 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col5 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col6 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col7 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col8 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col9 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col10 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col11 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col12 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col13 {
            background-color:  #2e8b57;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col14 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col15 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col16 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col17 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col18 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col19 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col0 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col1 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col2 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col3 {
            background-color:  #2e8b57;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col4 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col5 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col6 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col7 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col8 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col9 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col10 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col11 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col12 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col13 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col14 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col15 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col16 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col17 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col18 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col19 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col0 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col1 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col2 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col3 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col4 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col5 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col6 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col7 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col8 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col9 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col10 {
            background-color:  #2e8b57;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col11 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col12 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col13 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col14 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col15 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col16 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col17 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col18 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col19 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col0 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col1 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col2 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col3 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col4 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col5 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col6 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col7 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col8 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col9 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col10 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col11 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col12 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col13 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col14 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col15 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col16 {
            background-color:  #2e8b57;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col17 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col18 {
            background-color:  #ecf9f1;
        }    #T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col19 {
            background-color:  #ecf9f1;
        }</style>  
<table id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053" >
<thead>    <tr>
        <th class="blank level0" ></th>
        <th class="col_heading level0 col0" >A</th>
        <th class="col_heading level0 col1" >C</th>
        <th class="col_heading level0 col2" >D</th>
        <th class="col_heading level0 col3" >E</th>
        <th class="col_heading level0 col4" >F</th>
        <th class="col_heading level0 col5" >G</th>
        <th class="col_heading level0 col6" >H</th>
        <th class="col_heading level0 col7" >I</th>
        <th class="col_heading level0 col8" >K</th>
        <th class="col_heading level0 col9" >L</th>
        <th class="col_heading level0 col10" >M</th>
        <th class="col_heading level0 col11" >N</th>
        <th class="col_heading level0 col12" >P</th>
        <th class="col_heading level0 col13" >Q</th>
        <th class="col_heading level0 col14" >R</th>
        <th class="col_heading level0 col15" >S</th>
        <th class="col_heading level0 col16" >T</th>
        <th class="col_heading level0 col17" >V</th>
        <th class="col_heading level0 col18" >W</th>
        <th class="col_heading level0 col19" >Y</th>
    </tr></thead>
<tbody>    <tr>
        <th id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053level0_row0" class="row_heading level0 row0" >0</th>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col0" class="data row0 col0" >1</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col1" class="data row0 col1" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col2" class="data row0 col2" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col3" class="data row0 col3" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col4" class="data row0 col4" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col5" class="data row0 col5" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col6" class="data row0 col6" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col7" class="data row0 col7" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col8" class="data row0 col8" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col9" class="data row0 col9" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col10" class="data row0 col10" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col11" class="data row0 col11" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col12" class="data row0 col12" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col13" class="data row0 col13" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col14" class="data row0 col14" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col15" class="data row0 col15" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col16" class="data row0 col16" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col17" class="data row0 col17" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col18" class="data row0 col18" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row0_col19" class="data row0 col19" >0</td>
    </tr>    <tr>
        <th id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053level0_row1" class="row_heading level0 row1" >1</th>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col0" class="data row1 col0" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col1" class="data row1 col1" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col2" class="data row1 col2" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col3" class="data row1 col3" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col4" class="data row1 col4" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col5" class="data row1 col5" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col6" class="data row1 col6" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col7" class="data row1 col7" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col8" class="data row1 col8" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col9" class="data row1 col9" >1</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col10" class="data row1 col10" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col11" class="data row1 col11" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col12" class="data row1 col12" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col13" class="data row1 col13" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col14" class="data row1 col14" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col15" class="data row1 col15" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col16" class="data row1 col16" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col17" class="data row1 col17" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col18" class="data row1 col18" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row1_col19" class="data row1 col19" >0</td>
    </tr>    <tr>
        <th id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053level0_row2" class="row_heading level0 row2" >2</th>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col0" class="data row2 col0" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col1" class="data row2 col1" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col2" class="data row2 col2" >1</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col3" class="data row2 col3" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col4" class="data row2 col4" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col5" class="data row2 col5" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col6" class="data row2 col6" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col7" class="data row2 col7" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col8" class="data row2 col8" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col9" class="data row2 col9" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col10" class="data row2 col10" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col11" class="data row2 col11" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col12" class="data row2 col12" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col13" class="data row2 col13" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col14" class="data row2 col14" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col15" class="data row2 col15" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col16" class="data row2 col16" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col17" class="data row2 col17" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col18" class="data row2 col18" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row2_col19" class="data row2 col19" >0</td>
    </tr>    <tr>
        <th id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053level0_row3" class="row_heading level0 row3" >3</th>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col0" class="data row3 col0" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col1" class="data row3 col1" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col2" class="data row3 col2" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col3" class="data row3 col3" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col4" class="data row3 col4" >1</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col5" class="data row3 col5" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col6" class="data row3 col6" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col7" class="data row3 col7" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col8" class="data row3 col8" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col9" class="data row3 col9" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col10" class="data row3 col10" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col11" class="data row3 col11" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col12" class="data row3 col12" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col13" class="data row3 col13" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col14" class="data row3 col14" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col15" class="data row3 col15" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col16" class="data row3 col16" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col17" class="data row3 col17" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col18" class="data row3 col18" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row3_col19" class="data row3 col19" >0</td>
    </tr>    <tr>
        <th id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053level0_row4" class="row_heading level0 row4" >4</th>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col0" class="data row4 col0" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col1" class="data row4 col1" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col2" class="data row4 col2" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col3" class="data row4 col3" >1</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col4" class="data row4 col4" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col5" class="data row4 col5" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col6" class="data row4 col6" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col7" class="data row4 col7" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col8" class="data row4 col8" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col9" class="data row4 col9" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col10" class="data row4 col10" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col11" class="data row4 col11" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col12" class="data row4 col12" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col13" class="data row4 col13" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col14" class="data row4 col14" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col15" class="data row4 col15" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col16" class="data row4 col16" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col17" class="data row4 col17" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col18" class="data row4 col18" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row4_col19" class="data row4 col19" >0</td>
    </tr>    <tr>
        <th id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053level0_row5" class="row_heading level0 row5" >5</th>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col0" class="data row5 col0" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col1" class="data row5 col1" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col2" class="data row5 col2" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col3" class="data row5 col3" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col4" class="data row5 col4" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col5" class="data row5 col5" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col6" class="data row5 col6" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col7" class="data row5 col7" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col8" class="data row5 col8" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col9" class="data row5 col9" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col10" class="data row5 col10" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col11" class="data row5 col11" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col12" class="data row5 col12" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col13" class="data row5 col13" >1</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col14" class="data row5 col14" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col15" class="data row5 col15" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col16" class="data row5 col16" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col17" class="data row5 col17" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col18" class="data row5 col18" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row5_col19" class="data row5 col19" >0</td>
    </tr>    <tr>
        <th id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053level0_row6" class="row_heading level0 row6" >6</th>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col0" class="data row6 col0" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col1" class="data row6 col1" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col2" class="data row6 col2" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col3" class="data row6 col3" >1</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col4" class="data row6 col4" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col5" class="data row6 col5" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col6" class="data row6 col6" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col7" class="data row6 col7" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col8" class="data row6 col8" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col9" class="data row6 col9" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col10" class="data row6 col10" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col11" class="data row6 col11" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col12" class="data row6 col12" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col13" class="data row6 col13" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col14" class="data row6 col14" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col15" class="data row6 col15" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col16" class="data row6 col16" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col17" class="data row6 col17" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col18" class="data row6 col18" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row6_col19" class="data row6 col19" >0</td>
    </tr>    <tr>
        <th id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053level0_row7" class="row_heading level0 row7" >7</th>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col0" class="data row7 col0" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col1" class="data row7 col1" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col2" class="data row7 col2" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col3" class="data row7 col3" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col4" class="data row7 col4" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col5" class="data row7 col5" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col6" class="data row7 col6" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col7" class="data row7 col7" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col8" class="data row7 col8" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col9" class="data row7 col9" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col10" class="data row7 col10" >1</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col11" class="data row7 col11" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col12" class="data row7 col12" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col13" class="data row7 col13" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col14" class="data row7 col14" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col15" class="data row7 col15" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col16" class="data row7 col16" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col17" class="data row7 col17" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col18" class="data row7 col18" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row7_col19" class="data row7 col19" >0</td>
    </tr>    <tr>
        <th id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053level0_row8" class="row_heading level0 row8" >8</th>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col0" class="data row8 col0" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col1" class="data row8 col1" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col2" class="data row8 col2" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col3" class="data row8 col3" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col4" class="data row8 col4" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col5" class="data row8 col5" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col6" class="data row8 col6" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col7" class="data row8 col7" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col8" class="data row8 col8" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col9" class="data row8 col9" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col10" class="data row8 col10" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col11" class="data row8 col11" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col12" class="data row8 col12" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col13" class="data row8 col13" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col14" class="data row8 col14" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col15" class="data row8 col15" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col16" class="data row8 col16" >1</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col17" class="data row8 col17" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col18" class="data row8 col18" >0</td>
        <td id="T_618e1b9a_eddc_11e8_82f3_7085c28c1053row8_col19" class="data row8 col19" >0</td>
    </tr></tbody>
</table>

<br>

## NLF encoding

This method of encoding is detailed by Nanni and Lumini in their paper. It takes many physicochemical properties and transforms them using a Fisher Transform (similar to a PCA) creating a smaller set of features that can describe the amino acid just as well. There are 19 transformed features. This is available on the github link below if you want to try it. The result is shown below for our sample peptide ALDFEQEMT.


```python
#read the matrix a csv file on github
nlf = pd.read_csv('https://raw.githubusercontent.com/dmnfarrell/epitopepredict/master/epitopepredict/mhcdata/NLF.csv',index_col=0)

def nlf_encode(seq):    
    x = pd.DataFrame([nlf[i] for i in seq]).reset_index(drop=True)  
    show_matrix(x)
    e = x.values.flatten()
    return e

e=nlf_encode(pep)
```


<style  type="text/css" >
    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col0 {
            background-color:  #e1f3e9;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col1 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col2 {
            background-color:  #94c6aa;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col3 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col4 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col5 {
            background-color:  #a6d0b9;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col6 {
            background-color:  #9fccb3;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col7 {
            background-color:  #b5d9c5;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col8 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col9 {
            background-color:  #c2e1d0;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col10 {
            background-color:  #b7dac6;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col11 {
            background-color:  #81bb9a;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col12 {
            background-color:  #e3f3ea;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col13 {
            background-color:  #d6ecdf;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col14 {
            background-color:  #e3f3ea;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col15 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col16 {
            background-color:  #8dc2a4;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col17 {
            background-color:  #cde7d8;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col0 {
            background-color:  #91c4a8;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col1 {
            background-color:  #82bc9b;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col2 {
            background-color:  #d3eadd;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col3 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col4 {
            background-color:  #cfe8da;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col5 {
            background-color:  #e7f6ee;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col6 {
            background-color:  #e7f6ee;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col7 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col8 {
            background-color:  #a2ceb6;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col9 {
            background-color:  #4a9b6d;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col10 {
            background-color:  #81bb9b;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col11 {
            background-color:  #d3eadd;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col12 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col13 {
            background-color:  #b0d6c0;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col14 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col15 {
            background-color:  #3c9363;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col16 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col17 {
            background-color:  #d7ede1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col0 {
            background-color:  #bddecb;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col1 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col2 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col3 {
            background-color:  #70b18c;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col4 {
            background-color:  #d4ebde;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col5 {
            background-color:  #ddf0e5;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col6 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col7 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col8 {
            background-color:  #98c8ad;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col9 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col10 {
            background-color:  #98c8ae;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col11 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col12 {
            background-color:  #5ba57b;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col13 {
            background-color:  #a7d1ba;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col14 {
            background-color:  #b3d8c3;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col15 {
            background-color:  #77b592;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col16 {
            background-color:  #def1e7;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col17 {
            background-color:  #b7dac6;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col0 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col1 {
            background-color:  #e2f3ea;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col2 {
            background-color:  #ebf8f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col3 {
            background-color:  #a4cfb7;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col4 {
            background-color:  #cde7d8;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col5 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col6 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col7 {
            background-color:  #c7e4d4;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col8 {
            background-color:  #6db08a;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col9 {
            background-color:  #9eccb2;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col10 {
            background-color:  #bedecd;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col11 {
            background-color:  #e4f4eb;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col12 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col13 {
            background-color:  #e4f4eb;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col14 {
            background-color:  #a0cdb4;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col15 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col16 {
            background-color:  #9bcaaf;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col17 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col0 {
            background-color:  #78b694;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col1 {
            background-color:  #cae5d6;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col2 {
            background-color:  #76b592;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col3 {
            background-color:  #eaf8f0;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col4 {
            background-color:  #e1f3e9;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col5 {
            background-color:  #d6ecdf;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col6 {
            background-color:  #96c7ac;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col7 {
            background-color:  #90c4a7;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col8 {
            background-color:  #a2ceb6;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col9 {
            background-color:  #bbdcca;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col10 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col11 {
            background-color:  #d3eadd;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col12 {
            background-color:  #7fba99;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col13 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col14 {
            background-color:  #d9eee2;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col15 {
            background-color:  #94c6aa;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col16 {
            background-color:  #9bcaaf;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col17 {
            background-color:  #d7ede1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col0 {
            background-color:  #6aae88;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col1 {
            background-color:  #8cc1a3;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col2 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col3 {
            background-color:  #d0e9db;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col4 {
            background-color:  #dbefe4;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col5 {
            background-color:  #85bd9e;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col6 {
            background-color:  #c8e4d4;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col7 {
            background-color:  #bbdcca;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col8 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col9 {
            background-color:  #519f73;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col10 {
            background-color:  #b0d6c0;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col11 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col12 {
            background-color:  #37905e;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col13 {
            background-color:  #89bfa1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col14 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col15 {
            background-color:  #b2d7c2;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col16 {
            background-color:  #d1e9dc;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col17 {
            background-color:  #cde7d8;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col0 {
            background-color:  #78b694;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col1 {
            background-color:  #cae5d6;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col2 {
            background-color:  #76b592;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col3 {
            background-color:  #eaf8f0;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col4 {
            background-color:  #e1f3e9;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col5 {
            background-color:  #d6ecdf;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col6 {
            background-color:  #96c7ac;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col7 {
            background-color:  #90c4a7;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col8 {
            background-color:  #a2ceb6;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col9 {
            background-color:  #bbdcca;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col10 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col11 {
            background-color:  #d3eadd;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col12 {
            background-color:  #7fba99;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col13 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col14 {
            background-color:  #d9eee2;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col15 {
            background-color:  #94c6aa;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col16 {
            background-color:  #9bcaaf;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col17 {
            background-color:  #d7ede1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col0 {
            background-color:  #6aad87;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col1 {
            background-color:  #a5d0b8;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col2 {
            background-color:  #c5e2d2;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col3 {
            background-color:  #96c7ac;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col4 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col5 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col6 {
            background-color:  #b1d7c2;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col7 {
            background-color:  #cde7d9;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col8 {
            background-color:  #cde7d8;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col9 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col10 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col11 {
            background-color:  #cbe6d7;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col12 {
            background-color:  #daeee3;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col13 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col14 {
            background-color:  #d9eee2;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col15 {
            background-color:  #94c6aa;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col16 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col17 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col0 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col1 {
            background-color:  #b6dac6;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col2 {
            background-color:  #75b491;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col3 {
            background-color:  #c0dfce;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col4 {
            background-color:  #dcf0e5;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col5 {
            background-color:  #bcddcb;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col6 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col7 {
            background-color:  #e6f5ed;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col8 {
            background-color:  #6db08a;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col9 {
            background-color:  #90c4a7;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col10 {
            background-color:  #98c8ae;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col11 {
            background-color:  #57a378;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col12 {
            background-color:  #7fba99;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col13 {
            background-color:  #ddf0e5;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col14 {
            background-color:  #2e8b57;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col15 {
            background-color:  #b2d7c2;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col16 {
            background-color:  #ecf9f1;
        }    #T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col17 {
            background-color:  #c2e1d0;
        }</style>  
<table id="T_582c7baa_eddc_11e8_82f3_7085c28c1053" >
<thead>    <tr>
        <th class="blank level0" ></th>
        <th class="col_heading level0 col0" >1</th>
        <th class="col_heading level0 col1" >2</th>
        <th class="col_heading level0 col2" >3</th>
        <th class="col_heading level0 col3" >4</th>
        <th class="col_heading level0 col4" >5</th>
        <th class="col_heading level0 col5" >6</th>
        <th class="col_heading level0 col6" >7</th>
        <th class="col_heading level0 col7" >8</th>
        <th class="col_heading level0 col8" >9</th>
        <th class="col_heading level0 col9" >10</th>
        <th class="col_heading level0 col10" >11</th>
        <th class="col_heading level0 col11" >12</th>
        <th class="col_heading level0 col12" >13</th>
        <th class="col_heading level0 col13" >14</th>
        <th class="col_heading level0 col14" >15</th>
        <th class="col_heading level0 col15" >16</th>
        <th class="col_heading level0 col16" >17</th>
        <th class="col_heading level0 col17" >18</th>
    </tr></thead>
<tbody>    <tr>
        <th id="T_582c7baa_eddc_11e8_82f3_7085c28c1053level0_row0" class="row_heading level0 row0" >0</th>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col0" class="data row0 col0" >0.42</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col1" class="data row0 col1" >2.07</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col2" class="data row0 col2" >0.67</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col3" class="data row0 col3" >0.01</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col4" class="data row0 col4" >1.1</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col5" class="data row0 col5" >0.32</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col6" class="data row0 col6" >0.2</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col7" class="data row0 col7" >0.09</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col8" class="data row0 col8" >0.2</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col9" class="data row0 col9" >0.09</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col10" class="data row0 col10" >0.11</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col11" class="data row0 col11" >0.15</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col12" class="data row0 col12" >0.01</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col13" class="data row0 col13" >0.06</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col14" class="data row0 col14" >0.02</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col15" class="data row0 col15" >0.16</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col16" class="data row0 col16" >0.07</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row0_col17" class="data row0 col17" >0.03</td>
    </tr>    <tr>
        <th id="T_582c7baa_eddc_11e8_82f3_7085c28c1053level0_row1" class="row_heading level0 row1" >1</th>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col0" class="data row1 col0" >1.29</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col1" class="data row1 col1" >1.21</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col2" class="data row1 col2" >0.25</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col3" class="data row1 col3" >0.96</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col4" class="data row1 col4" >0.18</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col5" class="data row1 col5" >0.06</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col6" class="data row1 col6" >0.04</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col7" class="data row1 col7" >0</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col8" class="data row1 col8" >0.09</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col9" class="data row1 col9" >0.26</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col10" class="data row1 col10" >0.18</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col11" class="data row1 col11" >0.05</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col12" class="data row1 col12" >0</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col13" class="data row1 col13" >0.11</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col14" class="data row1 col14" >0.01</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col15" class="data row1 col15" >0.15</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col16" class="data row1 col16" >0.14</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row1_col17" class="data row1 col17" >0.02</td>
    </tr>    <tr>
        <th id="T_582c7baa_eddc_11e8_82f3_7085c28c1053level0_row2" class="row_heading level0 row2" >2</th>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col0" class="data row2 col0" >0.81</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col1" class="data row2 col1" >0.13</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col2" class="data row2 col2" >1.36</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col3" class="data row2 col3" >0.63</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col4" class="data row2 col4" >0.15</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col5" class="data row2 col5" >0.1</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col6" class="data row2 col6" >0.45</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col7" class="data row2 col7" >0.31</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col8" class="data row2 col8" >0.1</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col9" class="data row2 col9" >0.03</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col10" class="data row2 col10" >0.15</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col11" class="data row2 col11" >0.02</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col12" class="data row2 col12" >0.16</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col13" class="data row2 col13" >0.12</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col14" class="data row2 col14" >0.07</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col15" class="data row2 col15" >0.11</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col16" class="data row2 col16" >0.01</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row2_col17" class="data row2 col17" >0.05</td>
    </tr>    <tr>
        <th id="T_582c7baa_eddc_11e8_82f3_7085c28c1053level0_row3" class="row_heading level0 row3" >3</th>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col0" class="data row3 col0" >2.37</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col1" class="data row3 col1" >0.23</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col2" class="data row3 col2" >0.09</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col3" class="data row3 col3" >0.37</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col4" class="data row3 col4" >0.19</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col5" class="data row3 col5" >0.04</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col6" class="data row3 col6" >0.03</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col7" class="data row3 col7" >0.06</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col8" class="data row3 col8" >0.14</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col9" class="data row3 col9" >0.14</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col10" class="data row3 col10" >0.1</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col11" class="data row3 col11" >0.03</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col12" class="data row3 col12" >0.21</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col13" class="data row3 col13" >0.04</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col14" class="data row3 col14" >0.09</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col15" class="data row3 col15" >0.03</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col16" class="data row3 col16" >0.06</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row3_col17" class="data row3 col17" >0.18</td>
    </tr>    <tr>
        <th id="T_582c7baa_eddc_11e8_82f3_7085c28c1053level0_row4" class="row_heading level0 row4" >4</th>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col0" class="data row4 col0" >1.56</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col1" class="data row4 col1" >0.48</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col2" class="data row4 col2" >0.87</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col3" class="data row4 col3" >0.02</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col4" class="data row4 col4" >0.07</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col5" class="data row4 col5" >0.13</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col6" class="data row4 col6" >0.22</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col7" class="data row4 col7" >0.15</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col8" class="data row4 col8" >0.09</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col9" class="data row4 col9" >0.1</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col10" class="data row4 col10" >0.04</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col11" class="data row4 col11" >0.05</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col12" class="data row4 col12" >0.12</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col13" class="data row4 col13" >0.28</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col14" class="data row4 col14" >0.03</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col15" class="data row4 col15" >0.09</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col16" class="data row4 col16" >0.06</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row4_col17" class="data row4 col17" >0.02</td>
    </tr>    <tr>
        <th id="T_582c7baa_eddc_11e8_82f3_7085c28c1053level0_row5" class="row_heading level0 row5" >5</th>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col0" class="data row5 col0" >1.71</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col1" class="data row5 col1" >1.11</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col2" class="data row5 col2" >0.08</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col3" class="data row5 col3" >0.15</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col4" class="data row5 col4" >0.11</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col5" class="data row5 col5" >0.45</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col6" class="data row5 col6" >0.11</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col7" class="data row5 col7" >0.08</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col8" class="data row5 col8" >0.02</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col9" class="data row5 col9" >0.25</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col10" class="data row5 col10" >0.12</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col11" class="data row5 col11" >0.25</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col12" class="data row5 col12" >0.2</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col13" class="data row5 col13" >0.16</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col14" class="data row5 col14" >0.01</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col15" class="data row5 col15" >0.07</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col16" class="data row5 col16" >0.02</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row5_col17" class="data row5 col17" >0.03</td>
    </tr>    <tr>
        <th id="T_582c7baa_eddc_11e8_82f3_7085c28c1053level0_row6" class="row_heading level0 row6" >6</th>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col0" class="data row6 col0" >1.56</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col1" class="data row6 col1" >0.48</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col2" class="data row6 col2" >0.87</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col3" class="data row6 col3" >0.02</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col4" class="data row6 col4" >0.07</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col5" class="data row6 col5" >0.13</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col6" class="data row6 col6" >0.22</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col7" class="data row6 col7" >0.15</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col8" class="data row6 col8" >0.09</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col9" class="data row6 col9" >0.1</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col10" class="data row6 col10" >0.04</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col11" class="data row6 col11" >0.05</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col12" class="data row6 col12" >0.12</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col13" class="data row6 col13" >0.28</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col14" class="data row6 col14" >0.03</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col15" class="data row6 col15" >0.09</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col16" class="data row6 col16" >0.06</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row6_col17" class="data row6 col17" >0.02</td>
    </tr>    <tr>
        <th id="T_582c7baa_eddc_11e8_82f3_7085c28c1053level0_row7" class="row_heading level0 row7" >7</th>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col0" class="data row7 col0" >1.72</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col1" class="data row7 col1" >0.85</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col2" class="data row7 col2" >0.34</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col3" class="data row7 col3" >0.44</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col4" class="data row7 col4" >0.01</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col5" class="data row7 col5" >0.8</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col6" class="data row7 col6" >0.16</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col7" class="data row7 col7" >0.05</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col8" class="data row7 col8" >0.05</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col9" class="data row7 col9" >0.3</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col10" class="data row7 col10" >0.29</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col11" class="data row7 col11" >0.06</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col12" class="data row7 col12" >0.02</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col13" class="data row7 col13" >0.03</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col14" class="data row7 col14" >0.03</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col15" class="data row7 col15" >0.09</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col16" class="data row7 col16" >0.14</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row7_col17" class="data row7 col17" >0</td>
    </tr>    <tr>
        <th id="T_582c7baa_eddc_11e8_82f3_7085c28c1053level0_row8" class="row_heading level0 row8" >8</th>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col0" class="data row8 col0" >0.3</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col1" class="data row8 col1" >0.68</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col2" class="data row8 col2" >0.88</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col3" class="data row8 col3" >0.23</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col4" class="data row8 col4" >0.1</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col5" class="data row8 col5" >0.23</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col6" class="data row8 col6" >0.03</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col7" class="data row8 col7" >0.01</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col8" class="data row8 col8" >0.14</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col9" class="data row8 col9" >0.16</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col10" class="data row8 col10" >0.15</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col11" class="data row8 col11" >0.2</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col12" class="data row8 col12" >0.12</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col13" class="data row8 col13" >0.05</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col14" class="data row8 col14" >0.21</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col15" class="data row8 col15" >0.07</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col16" class="data row8 col16" >0</td>
        <td id="T_582c7baa_eddc_11e8_82f3_7085c28c1053row8_col17" class="data row8 col17" >0.04</td>
    </tr></tbody>
</table>

<br>

## Blosum encoding

BLOSUM62 is a substitution matrix that specifies the similarity of one amino acid to another by means of a score. This score reflects the frequency of substiutions found from studying protein sequence conservation in large databases of related proteins. The number 62 refers to the percentage identity at which sequences are clustered in the analysis. Encoding a peptide this way means we provide the column from the blosum matrix corresponding to the amino acid at each position of the sequence. This produces 21x9 matrix. see https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/blosum

The code to produce this array is shown below.


```python
blosum = ep.blosum62

def blosum_encode(seq):
    #encode a peptide into blosum features
    s=list(seq)
    x = pd.DataFrame([blosum[i] for i in seq]).reset_index(drop=True)
    show_matrix(x)
    e = x.values.flatten()    
    return e

def random_encode(p):
    return [np.random.randint(20) for i in pep]

e=blosum_encode(pep)
```

<style  type="text/css" >
    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col0 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col1 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col2 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col3 {
            background-color:  #c6e3d3;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col4 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col5 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col6 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col7 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col8 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col9 {
            background-color:  #a0cdb4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col10 {
            background-color:  #a4cfb7;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col11 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col12 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col13 {
            background-color:  #d7ede1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col14 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col15 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col16 {
            background-color:  #b5d9c5;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col17 {
            background-color:  #c6e3d3;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col18 {
            background-color:  #cde7d8;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col19 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col20 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col21 {
            background-color:  #b5d9c5;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col22 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col23 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col0 {
            background-color:  #cde7d8;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col1 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col2 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col3 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col4 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col5 {
            background-color:  #d4ebde;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col6 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col7 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col8 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col9 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col10 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col11 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col12 {
            background-color:  #75b491;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col13 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col14 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col15 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col16 {
            background-color:  #d1e9dc;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col17 {
            background-color:  #a0cdb4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col18 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col19 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col20 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col21 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col22 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col23 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col0 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col1 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col2 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col3 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col4 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col5 {
            background-color:  #a4cfb7;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col6 {
            background-color:  #75b491;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col7 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col8 {
            background-color:  #6db08a;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col9 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col10 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col11 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col12 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col13 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col14 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col15 {
            background-color:  #6db08a;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col16 {
            background-color:  #d1e9dc;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col17 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col18 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col19 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col20 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col21 {
            background-color:  #7fba99;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col22 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col23 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col0 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col1 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col2 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col3 {
            background-color:  #d9eee2;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col4 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col5 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col6 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col7 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col8 {
            background-color:  #6db08a;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col9 {
            background-color:  #7ab795;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col10 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col11 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col12 {
            background-color:  #a4cfb7;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col13 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col14 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col15 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col16 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col17 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col18 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col19 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col20 {
            background-color:  #d4ebde;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col21 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col22 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col23 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col0 {
            background-color:  #cde7d8;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col1 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col2 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col3 {
            background-color:  #7ab795;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col4 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col5 {
            background-color:  #75b491;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col6 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col7 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col8 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col9 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col10 {
            background-color:  #d4ebde;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col11 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col12 {
            background-color:  #d4ebde;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col13 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col14 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col15 {
            background-color:  #6db08a;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col16 {
            background-color:  #d1e9dc;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col17 {
            background-color:  #c6e3d3;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col18 {
            background-color:  #cde7d8;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col19 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col20 {
            background-color:  #75b491;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col21 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col22 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col23 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col0 {
            background-color:  #cde7d8;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col1 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col2 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col3 {
            background-color:  #a0cdb4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col4 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col5 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col6 {
            background-color:  #75b491;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col7 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col8 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col9 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col10 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col11 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col12 {
            background-color:  #a4cfb7;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col13 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col14 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col15 {
            background-color:  #6db08a;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col16 {
            background-color:  #d1e9dc;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col17 {
            background-color:  #a0cdb4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col18 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col19 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col20 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col21 {
            background-color:  #499a6d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col22 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col23 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col0 {
            background-color:  #cde7d8;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col1 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col2 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col3 {
            background-color:  #7ab795;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col4 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col5 {
            background-color:  #75b491;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col6 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col7 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col8 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col9 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col10 {
            background-color:  #d4ebde;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col11 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col12 {
            background-color:  #d4ebde;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col13 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col14 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col15 {
            background-color:  #6db08a;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col16 {
            background-color:  #d1e9dc;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col17 {
            background-color:  #c6e3d3;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col18 {
            background-color:  #cde7d8;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col19 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col20 {
            background-color:  #75b491;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col21 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col22 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col23 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col0 {
            background-color:  #cde7d8;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col1 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col2 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col3 {
            background-color:  #d9eee2;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col4 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col5 {
            background-color:  #a4cfb7;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col6 {
            background-color:  #d4ebde;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col7 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col8 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col9 {
            background-color:  #54a176;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col10 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col11 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col12 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col13 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col14 {
            background-color:  #6db08a;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col15 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col16 {
            background-color:  #d1e9dc;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col17 {
            background-color:  #7ab795;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col18 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col19 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col20 {
            background-color:  #d4ebde;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col21 {
            background-color:  #b5d9c5;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col22 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col23 {
            background-color:  #ecf9f1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col0 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col1 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col2 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col3 {
            background-color:  #b3d8c3;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col4 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col5 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col6 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col7 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col8 {
            background-color:  #add4be;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col9 {
            background-color:  #a0cdb4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col10 {
            background-color:  #a4cfb7;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col11 {
            background-color:  #8dc2a4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col12 {
            background-color:  #bcddcb;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col13 {
            background-color:  #d7ede1;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col14 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col15 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col16 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col17 {
            background-color:  #a0cdb4;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col18 {
            background-color:  #cde7d8;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col19 {
            background-color:  #5da67d;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col20 {
            background-color:  #a4cfb7;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col21 {
            background-color:  #b5d9c5;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col22 {
            background-color:  #2e8b57;
        }    #T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col23 {
            background-color:  #ecf9f1;
        }</style>  
<table id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053" >
<thead>    <tr>
        <th class="blank level0" ></th>
        <th class="col_heading level0 col0" >A</th>
        <th class="col_heading level0 col1" >R</th>
        <th class="col_heading level0 col2" >N</th>
        <th class="col_heading level0 col3" >D</th>
        <th class="col_heading level0 col4" >C</th>
        <th class="col_heading level0 col5" >Q</th>
        <th class="col_heading level0 col6" >E</th>
        <th class="col_heading level0 col7" >G</th>
        <th class="col_heading level0 col8" >H</th>
        <th class="col_heading level0 col9" >I</th>
        <th class="col_heading level0 col10" >L</th>
        <th class="col_heading level0 col11" >K</th>
        <th class="col_heading level0 col12" >M</th>
        <th class="col_heading level0 col13" >F</th>
        <th class="col_heading level0 col14" >P</th>
        <th class="col_heading level0 col15" >S</th>
        <th class="col_heading level0 col16" >T</th>
        <th class="col_heading level0 col17" >W</th>
        <th class="col_heading level0 col18" >Y</th>
        <th class="col_heading level0 col19" >V</th>
        <th class="col_heading level0 col20" >B</th>
        <th class="col_heading level0 col21" >Z</th>
        <th class="col_heading level0 col22" >X</th>
        <th class="col_heading level0 col23" >*</th>
    </tr></thead>
<tbody>    <tr>
        <th id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053level0_row0" class="row_heading level0 row0" >0</th>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col0" class="data row0 col0" >4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col1" class="data row0 col1" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col2" class="data row0 col2" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col3" class="data row0 col3" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col4" class="data row0 col4" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col5" class="data row0 col5" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col6" class="data row0 col6" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col7" class="data row0 col7" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col8" class="data row0 col8" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col9" class="data row0 col9" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col10" class="data row0 col10" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col11" class="data row0 col11" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col12" class="data row0 col12" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col13" class="data row0 col13" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col14" class="data row0 col14" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col15" class="data row0 col15" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col16" class="data row0 col16" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col17" class="data row0 col17" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col18" class="data row0 col18" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col19" class="data row0 col19" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col20" class="data row0 col20" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col21" class="data row0 col21" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col22" class="data row0 col22" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row0_col23" class="data row0 col23" >-4</td>
    </tr>    <tr>
        <th id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053level0_row1" class="row_heading level0 row1" >1</th>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col0" class="data row1 col0" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col1" class="data row1 col1" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col2" class="data row1 col2" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col3" class="data row1 col3" >-4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col4" class="data row1 col4" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col5" class="data row1 col5" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col6" class="data row1 col6" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col7" class="data row1 col7" >-4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col8" class="data row1 col8" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col9" class="data row1 col9" >2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col10" class="data row1 col10" >4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col11" class="data row1 col11" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col12" class="data row1 col12" >2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col13" class="data row1 col13" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col14" class="data row1 col14" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col15" class="data row1 col15" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col16" class="data row1 col16" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col17" class="data row1 col17" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col18" class="data row1 col18" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col19" class="data row1 col19" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col20" class="data row1 col20" >-4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col21" class="data row1 col21" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col22" class="data row1 col22" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row1_col23" class="data row1 col23" >-4</td>
    </tr>    <tr>
        <th id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053level0_row2" class="row_heading level0 row2" >2</th>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col0" class="data row2 col0" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col1" class="data row2 col1" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col2" class="data row2 col2" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col3" class="data row2 col3" >6</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col4" class="data row2 col4" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col5" class="data row2 col5" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col6" class="data row2 col6" >2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col7" class="data row2 col7" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col8" class="data row2 col8" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col9" class="data row2 col9" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col10" class="data row2 col10" >-4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col11" class="data row2 col11" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col12" class="data row2 col12" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col13" class="data row2 col13" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col14" class="data row2 col14" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col15" class="data row2 col15" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col16" class="data row2 col16" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col17" class="data row2 col17" >-4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col18" class="data row2 col18" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col19" class="data row2 col19" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col20" class="data row2 col20" >4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col21" class="data row2 col21" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col22" class="data row2 col22" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row2_col23" class="data row2 col23" >-4</td>
    </tr>    <tr>
        <th id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053level0_row3" class="row_heading level0 row3" >3</th>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col0" class="data row3 col0" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col1" class="data row3 col1" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col2" class="data row3 col2" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col3" class="data row3 col3" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col4" class="data row3 col4" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col5" class="data row3 col5" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col6" class="data row3 col6" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col7" class="data row3 col7" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col8" class="data row3 col8" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col9" class="data row3 col9" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col10" class="data row3 col10" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col11" class="data row3 col11" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col12" class="data row3 col12" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col13" class="data row3 col13" >6</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col14" class="data row3 col14" >-4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col15" class="data row3 col15" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col16" class="data row3 col16" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col17" class="data row3 col17" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col18" class="data row3 col18" >3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col19" class="data row3 col19" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col20" class="data row3 col20" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col21" class="data row3 col21" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col22" class="data row3 col22" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row3_col23" class="data row3 col23" >-4</td>
    </tr>    <tr>
        <th id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053level0_row4" class="row_heading level0 row4" >4</th>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col0" class="data row4 col0" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col1" class="data row4 col1" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col2" class="data row4 col2" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col3" class="data row4 col3" >2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col4" class="data row4 col4" >-4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col5" class="data row4 col5" >2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col6" class="data row4 col6" >5</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col7" class="data row4 col7" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col8" class="data row4 col8" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col9" class="data row4 col9" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col10" class="data row4 col10" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col11" class="data row4 col11" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col12" class="data row4 col12" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col13" class="data row4 col13" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col14" class="data row4 col14" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col15" class="data row4 col15" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col16" class="data row4 col16" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col17" class="data row4 col17" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col18" class="data row4 col18" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col19" class="data row4 col19" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col20" class="data row4 col20" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col21" class="data row4 col21" >4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col22" class="data row4 col22" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row4_col23" class="data row4 col23" >-4</td>
    </tr>    <tr>
        <th id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053level0_row5" class="row_heading level0 row5" >5</th>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col0" class="data row5 col0" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col1" class="data row5 col1" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col2" class="data row5 col2" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col3" class="data row5 col3" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col4" class="data row5 col4" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col5" class="data row5 col5" >5</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col6" class="data row5 col6" >2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col7" class="data row5 col7" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col8" class="data row5 col8" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col9" class="data row5 col9" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col10" class="data row5 col10" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col11" class="data row5 col11" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col12" class="data row5 col12" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col13" class="data row5 col13" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col14" class="data row5 col14" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col15" class="data row5 col15" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col16" class="data row5 col16" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col17" class="data row5 col17" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col18" class="data row5 col18" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col19" class="data row5 col19" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col20" class="data row5 col20" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col21" class="data row5 col21" >3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col22" class="data row5 col22" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row5_col23" class="data row5 col23" >-4</td>
    </tr>    <tr>
        <th id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053level0_row6" class="row_heading level0 row6" >6</th>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col0" class="data row6 col0" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col1" class="data row6 col1" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col2" class="data row6 col2" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col3" class="data row6 col3" >2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col4" class="data row6 col4" >-4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col5" class="data row6 col5" >2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col6" class="data row6 col6" >5</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col7" class="data row6 col7" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col8" class="data row6 col8" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col9" class="data row6 col9" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col10" class="data row6 col10" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col11" class="data row6 col11" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col12" class="data row6 col12" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col13" class="data row6 col13" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col14" class="data row6 col14" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col15" class="data row6 col15" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col16" class="data row6 col16" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col17" class="data row6 col17" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col18" class="data row6 col18" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col19" class="data row6 col19" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col20" class="data row6 col20" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col21" class="data row6 col21" >4</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col22" class="data row6 col22" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row6_col23" class="data row6 col23" >-4</td>
    </tr>    <tr>
        <th id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053level0_row7" class="row_heading level0 row7" >7</th>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col0" class="data row7 col0" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col1" class="data row7 col1" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col2" class="data row7 col2" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col3" class="data row7 col3" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col4" class="data row7 col4" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col5" class="data row7 col5" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col6" class="data row7 col6" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col7" class="data row7 col7" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col8" class="data row7 col8" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col9" class="data row7 col9" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col10" class="data row7 col10" >2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col11" class="data row7 col11" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col12" class="data row7 col12" >5</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col13" class="data row7 col13" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col14" class="data row7 col14" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col15" class="data row7 col15" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col16" class="data row7 col16" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col17" class="data row7 col17" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col18" class="data row7 col18" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col19" class="data row7 col19" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col20" class="data row7 col20" >-3</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col21" class="data row7 col21" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col22" class="data row7 col22" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row7_col23" class="data row7 col23" >-4</td>
    </tr>    <tr>
        <th id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053level0_row8" class="row_heading level0 row8" >8</th>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col0" class="data row8 col0" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col1" class="data row8 col1" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col2" class="data row8 col2" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col3" class="data row8 col3" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col4" class="data row8 col4" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col5" class="data row8 col5" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col6" class="data row8 col6" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col7" class="data row8 col7" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col8" class="data row8 col8" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col9" class="data row8 col9" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col10" class="data row8 col10" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col11" class="data row8 col11" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col12" class="data row8 col12" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col13" class="data row8 col13" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col14" class="data row8 col14" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col15" class="data row8 col15" >1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col16" class="data row8 col16" >5</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col17" class="data row8 col17" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col18" class="data row8 col18" >-2</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col19" class="data row8 col19" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col20" class="data row8 col20" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col21" class="data row8 col21" >-1</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col22" class="data row8 col22" >0</td>
        <td id="T_5d8b89ec_eddc_11e8_82f3_7085c28c1053row8_col23" class="data row8 col23" >-4</td>
    </tr></tbody>
</table>

<br>

## Create a regression model to fit encoded features

We are going to create a regression model that can fit the encoded peptide featuers to known affinity data. This model is trained on the known data and then used to predict new peptides. The data used for training is primarily from the IEDB and was curated by the authors of MHCflurry from various sources. A sample of the data is shown below. The affinity measure is usually given as IC50 which is  This scale is harder for regression so we linearize the scale using the equation log50k = 1-log(ic50) / log(50000).



```python
df = ep.get_training_set()
display(df[:5])
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th>allele</th>
      <th>peptide</th>
      <th>ic50</th>
      <th>measurement_inequality</th>
      <th>measurement_type</th>
      <th>measurement_source</th>
      <th>original_allele</th>
      <th>log50k</th>
      <th>length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>BoLA-1*21:01</td>
      <td>AENDTLVVSV</td>
      <td>7817.0</td>
      <td>=</td>
      <td>quantitative</td>
      <td>Barlow - purified MHC/competitive/fluorescence</td>
      <td>BoLA-1*02101</td>
      <td>0.171512</td>
      <td>10</td>
    </tr>
    <tr>
      <th>1</th>
      <td>BoLA-1*21:01</td>
      <td>NQFNGGCLLV</td>
      <td>1086.0</td>
      <td>=</td>
      <td>quantitative</td>
      <td>Barlow - purified MHC/direct/fluorescence</td>
      <td>BoLA-1*02101</td>
      <td>0.353937</td>
      <td>10</td>
    </tr>
    <tr>
      <th>2</th>
      <td>BoLA-2*08:01</td>
      <td>AAHCIHAEW</td>
      <td>21.0</td>
      <td>=</td>
      <td>quantitative</td>
      <td>Barlow - purified MHC/direct/fluorescence</td>
      <td>BoLA-2*00801</td>
      <td>0.718615</td>
      <td>9</td>
    </tr>
    <tr>
      <th>3</th>
      <td>BoLA-2*08:01</td>
      <td>AAKHMSNTY</td>
      <td>1299.0</td>
      <td>=</td>
      <td>quantitative</td>
      <td>Barlow - purified MHC/direct/fluorescence</td>
      <td>BoLA-2*00801</td>
      <td>0.337385</td>
      <td>9</td>
    </tr>
    <tr>
      <th>4</th>
      <td>BoLA-2*08:01</td>
      <td>DSYAYMRNGW</td>
      <td>2.0</td>
      <td>=</td>
      <td>quantitative</td>
      <td>Barlow - purified MHC/direct/fluorescence</td>
      <td>BoLA-2*00801</td>
      <td>0.935937</td>
      <td>10</td>
    </tr>
  </tbody>
</table>
</div>

<br>

We wish to create a regression model to fit the peptide features against the true affinity values. The code below sets up a train/test scheme for a neural network regressor using scikit-learn. scikit-learn provides many possible models that can be created with relatively little knowledge of machine learning. I chose a neural network here because it seemed most appropriate after some experimentation. The parameters used were arrived at through a grid search which is not shown here as it is somewhat beyond the scope of this article. We then make scatter plots comparing the test to predicted values. The diagonal line is the ideal one to one value.


```python
from sklearn import metrics
from sklearn.model_selection import train_test_split,cross_val_score,ShuffleSplit
from sklearn.neural_network import MLPRegressor
import epitopepredict as ep

def auc_score(true,sc,cutoff=None):

    if cutoff!=None:
        true = (true<=cutoff).astype(int)
        sc = (sc<=cutoff).astype(int)        
    fpr, tpr, thresholds = metrics.roc_curve(true, sc, pos_label=1)
    r = metrics.auc(fpr, tpr)
    #print (r)
    return  r

def test_predictor(allele, encoder, ax):

    reg = MLPRegressor(hidden_layer_sizes=(20), alpha=0.01, max_iter=500,
                        activation='relu', solver='lbfgs', random_state=2)
    df = ep.get_training_set(allele, length=9)
    #print (len(df))
    X = df.peptide.apply(lambda x: pd.Series(encoder(x)),1)
    y = df.log50k
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1)
    reg.fit(X_train, y_train)
    sc = reg.predict(X_test)
    x=pd.DataFrame(np.column_stack([y_test,sc]),columns=['test','predicted'])
    x.plot('test','predicted',kind='scatter',s=20,ax=ax)
    ax.plot((0,1), (0,1), ls="--", lw=2, c=".2")
    ax.set_xlim((0,1));  ax.set_ylim((0,1))
    ax.set_title(encoder.__name__)    
    auc = ep.auc_score(y_test,sc,cutoff=.426)
    ax.text(.1,.9,'auc=%s' %round(auc,2))
    sns.despine()

sns.set_context('notebook')
encs=[blosum_encode,nlf_encode,one_hot_encode,random_encode]
allele='HLA-A*03:01'
fig,axs=plt.subplots(2,2,figsize=(10,10))
axs=axs.flat
i=0
for enc in encs:
    test_predictor(allele,enc,ax=axs[i])
    i+=1
```

<div style="width: 400px;">
<img src="/img/peptide_encoding_plots.png" width="600px">
</div>

## Build a predictor

It can be seen from the plots that blosum_encoding gives superior results (highest auc value) for this particular test/train set and allele. Repeating the process generally shows the same outcome. Now that we know how to make encode the peptides and make a model we can do it for any allele. We can write a routine that build the predictor from any training set.


```python
def build_predictor(allele, encoder):

    data = ep.get_training_set(allele)
    if len(data)<200:
        return
    from sklearn.neural_network import MLPRegressor
    reg = MLPRegressor(hidden_layer_sizes=(20), alpha=0.01, max_iter=500,
                        activation='relu', solver='lbfgs', random_state=2)    
    X = data.peptide.apply(lambda x: pd.Series(encoder(x)),1)
    y = data.log50k
    print (allele, len(X))
    reg.fit(X,y)       
    return reg
```

## Train and save models

Now we can just use the code above for all alleles in which we have training data (>200 samples) and produce a model for each one. Each model is saved to disk for later use so we don't have to re-train every time we want to predict peptides for that allele. sklearn uses the joblib library to persist models, which is similar to the pickle module.


```python
def get_allele_names():
    d = ep.get_training_set(length=9)
    a = d.allele.value_counts()
    a =a[a>200]
    return list(a.index)

al = get_allele_names()
path = 'models'
for a in al:
    fname = os.path.join(path, a+'.joblib')
    reg = build_predictor(a, blosum_encode)
    if reg is not None:
        joblib.dump(reg, fname, protocol=2)
```

## Comparison to other prediction algorithms

The two best MHC-class I prediction tools are currently netMHC/netMHCpan-4.0 and MHCFlurry. It is useful to compare our results to those. In order to do this I have implemented the algorithm we made (which uses the code above) as a predictor in the `epitopepredict` package. The object is created using `P = ep.get_predictor('basicmhc1')`. This standardizes the calls to prediction and the results returned so it can be directly compared to the other tools. The code below creates 3 predictor objects and evaluates their performance on some new data for each allele available. The score is then recoreded each time. Importantly, the peptides in the evaluation set are not present in the training set.

Note that predictors all return an ic50 value but internally use the log50k value for fitting. The predictors are evaluated using the roc auc metric with a threshold of 500nM. The auc is a common metric for classiifcation and can be used for regression if a threshold is chosen. Though others measures may be used such as pearson correlation co-efficient.


```python
def evaluate_predictor(P, allele):

    data = ep.get_evaluation_set(allele, length=9)
    #print (len(data))
    P.predict_peptides(list(data.peptide), alleles=allele, cpus=4)
    x = P.get_scores(allele)
    x = data.merge(x,on='peptide')
    auc = auc_score(x.ic50,x.score,cutoff=500)
    return auc, data

preds = [ep.get_predictor('basicmhc1'),
         ep.get_predictor('netmhcpan',scoring='affinity'),
         ep.get_predictor('mhcflurry')]
comp=[]
evalset = ep.get_evaluation_set(length=9)
test_alleles = evalset.allele.unique()

for P in preds:    
    m=[]
    for a in test_alleles:        
        try:
            auc,df = evaluate_predictor(P, a)
            m.append((a,auc,len(df)))            
        except Exception as e:
            print (a,e)
            pass
    m=pd.DataFrame(m,columns=['allele','score','size'])
    m['name'] = P.name
    comp.append(m)
```

## Results

Finally we use the resultant dataframe to plot the auc scores for each method per allele.


```python
c=pd.concat(comp)
x=pd.pivot_table(c,index=['allele','size'],columns='name',values='score')#.reset_index()
def highlight_max(s):
    is_max = s == s.max()
    return ['background-color: yellow' if v else '' for v in is_max]
#display(x.style.apply(highlight_max,1))
#print(c)

ax=sns.boxplot(data=c,y='score',x='name')#,hue='allele')
g=sns.catplot(data=c,y='score',x='allele',hue='name',
              kind='bar',aspect=3,height=5,legend=False)
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.setp(g.ax.get_xticklabels(), rotation=90);
plt.tight_layout()
plt.savefig('benchmarks.png')
x.to_csv('benchmarks.csv')
```

<div style="width: 400px;">
<img src="/img/mhc_benchmarks1.png" width="400px">
</div>

<div style="width: 800px;">
<img src="/img/mhc_benchmarks2.png" width="800px">
</div>


You can see clearly that our method is quite inferior to the other two for many of the alleles! Our model is limited in several ways:

* it only works with 9-mers
* some alleles have preference for different peptide lengths and this is not accounted for
* limited training data for some alleles, the other tools use methods to overcome this
* the neural network is not likely not sophisticated enough

Why the blosum62 encoding is better than the others tested is not entirely clear. However it should be better than one hot encoding since there is information loss when simply using 0 and 1 in a sparse matrix. As mentioned the predictor is available from within the epitopepredict package. It is used as a basic simple model only. To create the basic MHC predictor in epitopepredict we would use the code below. In the future this approach will be used as a basis for improving the predictor. Improvements would include handling different length peptides and changing the regression model.

## Using epitopepredict

```python
P = ep.get_predictor('basicmhc1')
from epitopepredict import peptutils
seqs = peptutils.create_random_sequences(10)
df = pd.DataFrame(seqs,columns=['peptide'])
res = P.predict_peptides(df.peptide, alleles=ep.get_preset_alleles('mhc1_supertypes'), cpus=1)
```

## References

* L. Nanni and A. Lumini, “A new encoding technique for peptide classification,” Expert Syst. Appl., vol. 38, no. 4, pp. 3185–3191, 2011.
* V. Jurtz, S. Paul, M. Andreatta, P. Marcatili, B. Peters, and M. Nielsen, “NetMHCpan-4.0: Improved Peptide–MHC Class I Interaction Predictions Integrating Eluted Ligand and Peptide Binding Affinity Data,” J. Immunol., vol. 199, no. 9, 2017.
* T. J. O’Donnell, A. Rubinsteyn, M. Bonsack, A. B. Riemer, U. Laserson, and J. Hammerbacher, “MHCflurry: Open-Source Class I MHC Binding Affinity Prediction,” Cell Syst., vol. 7, no. 1, p. 129–132.e4, 2018.
* J. Hu and Z. Liu, “DeepMHC : Deep Convolutional Neural Networks for High-performance peptide-MHC Binding Affinity Prediction,” bioRxiv, pp. 1–20, 2017.
* https://www.iedb.org/
