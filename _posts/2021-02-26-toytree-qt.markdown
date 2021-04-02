---
layout: post
title:  "A phylogenetic tree viewer with Qt and Toytree"
date:   2021-02-26 10:18:00
categories: software
tags: [pyside2,toytree]
thumbnail: /img/toytree_qt.png
---

## Background

<div style="width: 320px; float:right;">
<img src="/img/toytree_qt.png" width="280px">
</div>

Toytree is a Python library for tree manipulation and plotting. So it can be used to view phylogenetic trees. The goal of Toytree is to provide a light-weight Python equivalent to commonly used tree manipulation and plotting libraries in R. It's based on [toyplot](https://toyplot.readthedocs.io/) which uses Javascript and html, making plots self-contained and embeddable in browsers.

## Implementation

We can exploit these features to create an interactive tree viewer with the Qt toolkit in Python. This uses the `QWebEngineView()`, a widget that is used to view and edit web documents. Since the plots are javascript we can simply embed the html directly in the widget. The view can then be updated whenever we replot. Toytree objects support many styling attributes. A style dictionary can be created and then passed to the `draw` method when updating. This dict is populated from a widget using a special function included in the code that creates a dialog from a set of widget options. We also add a list view at the bottom that allows individual tips to be selected and colored or set as the root node. The checkboxes shown aren't functional here as yet but they could be used to make tips invisible.

The code is available [here](https://github.com/dmnfarrell/teaching/blob/master/gui/toytreeviewer.py). You can download this python file and run it from the command line. This will work with either pyside2 or pyqt5 which can be installed with pip. To install dependencies use: `pip install toytree pyside2 numpy`

<div style="width: auto;">
 <a href="/img/toytree_viewer.gif"> <img class="small-scaled" src="/img/toytree_viewer.gif"></a>  
</div>

## Links

* [Toytree](https://toytree.readthedocs.io/)
