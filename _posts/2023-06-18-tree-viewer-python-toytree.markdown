---
layout: post
title:  "A phylogenetic tree viewer with PyQt and ToyTree"
date:   2023-06-18 13:30:00
categories: bioinformatics
tags: [genomics,python]
thumbnail: /img/phylo_tree.png
---

## Background

[Toyplot](https://toyplot.readthedocs.io/) is a Python toolkit for plotting using HTML, SVG, and Javascript to create embeddable graphs. [Toytree](https://toytree.readthedocs.io/) is based on toyplot and provides the ability to plot tree structures. It is a fairly minimalist package but can create useful tree plots that are better than some alternatives in Python. There are relatively few packages in Python that support rendering of phylogenetic trees that are much use. I have tried to use ete3 but found it frustrating and overly complex. Biopython is the opposite, it provides very basic functionality that has not been extended. Basically none of these can match the equivalents in R like ape, phytools or ggtree.

This post describes a simple GUI tool that uses Toytree and PyQt to plot trees and color the nodes according to associated meta data that is loaded from a csv file. This allows basic annotation with categorical or continuous values matched to each tip in the tree. The trees can be loaded from newick format files. The program can be loaded with a test tree and data for testing.

## Code

Currently this program is included as part of the [SNiPgenie](https://github.com/dmnfarrell/snipgenie) package where it can be invoked using the `snipgenie-treeviewer` command. 

The `TreeViewer` class is a QWidget so you can also add this tool to your own application by importing the module and creating the widget like this:

```python
from snipgenie import treeview
tv = treeview.TreeViewer(tree='tree.newick', meta='meta.csv')
#load tree separately if needed
tv.load_tree('tree.newick')
#redraw tree 
tv.update()
#add the widget to your app as you require e.g:
parent.addWidget(tv)
```

## Example

<div style="width: auto;">
 <a href="/img/tree_viewer.gif"> <img class="scaled" src="/img/tree_viewer.gif"></a>  
  <p class="caption"></p>
</div>

## Links

* [code on github](https://github.com/dmnfarrell/snipgenie/blob/master/snipgenie/treeview.py)
