---
layout: post
title:  "Plot phylogenies with annotation in R using ggtree and gheatmap"
date:   2022-08-02 12:19:00
categories: R
tags: [phylogeny,plotting]
thumbnail: /img/ggtree_example.png
---

## Background

There are lots of online examples of how to draw phylogenetic trees using various R tools. One is ggtree, based on the ggplot packages, which provides a side range of annotation and layouts. This example shows how to write some functions that can plot trees with an arbitrary number of heatmap annotations, given the correct meta data. The columns to be used are provided as a list along with colormaps for each. The first column is used for the tip colors. Note that you will see a warning: _Scale for 'fill' is already present. Adding another scale for 'fill', which will replace the existing scale._ when this is run. This is because a new `scale_fill_manual` is called each time we add a heatmap. It can be ignored.

## Code

```R
library("ape")
library(RColorBrewer)
library(dplyr)
library('ggplot2')
library('ggtree')
library(tidytree)
library(ggnewscale)

gettreedata <- function(tree, meta){
    #get treedata object
    d<-meta[row.names(meta) %in% tree$tip.label,]
    d$label <- row.names(d)
    y <- full_join(as_tibble(tree), d, by='label')
    y <- as.treedata(y)
    return(y)
}

get_color_mapping <- function(data, col, cmap){
    #get a color mapping
    labels <- (data[[col]])  
    names <- levels(as.factor(labels))
    n <- length(names)
    if (n<10){      
        colors <- suppressWarnings(c(brewer.pal(n, cmap)))[1:n]
    }
    else {
        colors <- colorRampPalette(brewer.pal(8, cmap))(n)
    }
    names(colors) = names
    return (colors)
}

ggplottree <- function(tree, meta, cols=NULL, cmaps=NULL, layout="rectangular",
                       offset=10, tiplabel=FALSE, tipsize=3,
                       title='') {

    y <- gettreedata(tree, meta)
    p <- ggtree(y, layout=layout)
    if (is.null(cols)){
        return (p)
    }

    col <- cols[1]
    cmap <- cmaps[1]
    df<-meta[tree$tip.label,][col]
    colors <- get_color_mapping(df, col, cmap)

    #tip formatting    
    p1 <- p + new_scale_fill() +    
          geom_tippoint(mapping=aes(fill=.data[[col]]),size=tipsize,shape=21) +
          scale_fill_manual(values=colors, na.value="white")

    p2 <- p1  
    #heatmaps
    if (length(cols)>1){
        for (i in 2:length(cols)){
            col <- cols[i]
            cmap <- cmaps[i]
            df <- meta[tree$tip.label,][col]            
            colors <- get_color_mapping(df, col, cmap)       
            p2 <- p2 + new_scale_fill()
            p2 <- gheatmap(p2, df, offset=i*offset, width=.08,
                      colnames_angle=0, colnames_offset_y = .05)  +
                  scale_fill_manual(values=colors, name=col)
        }
    }

    p2 <- p2 + theme_tree2(legend.text = element_text(size=20), legend.key.size = unit(1, 'cm'),
                        legend.position="left", plot.title = element_text(size=40))    
    return(p2)
}
```

## Usage

A simple tree and meta data is created in the code below for illustration. We then plot by calling the above function.

```R
tree <- ape::read.tree(text='((A, B), ((C, D), (E, F)));')
options(repr.plot.width=8, repr.plot.height=6)
df <- data.frame (name = c('A','B','C','D','E','F'),
                  label1 = c('X','X','Y','Y','Y','Y'),
                  species = c('dog','cat','dog','dog','cat','cat'),
                  country = c('Ireland','Ireland','France','UK','France','UK')
                  )
row.names(df) <- df$name
#call the function
ggplottree(tree, df, cols=c('label1','species','country'),
          cmaps=c('Set1','Set2','Set3'), tipsize=8, offset=.5)

```

Which makes a table like this to use as the label data:

```
    name	label1	species	country
A	A	X	dog	Ireland
B	B	X	cat	Ireland
C	C	Y	dog	France
D	D	Y	dog	UK
E	E	Y	cat	France
F	F	Y	cat	UK
```

Finally, the tree looks like this:

<div style="width: auto;">
 <img class="small-scaled" src="/img/ggtree_example.png">
   <p class="caption">Example tree.</p>
</div>

<div style="width: auto;">
 <img class="small-scaled" src="/img/ggtree_example_circ.png">
   <p class="caption">Circular layout.</p>
</div>

## Links

* https://yulab-smu.top/treedata-book/chapter10.html
