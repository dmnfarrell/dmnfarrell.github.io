---
layout: post
title:  "Plot a phylogenetic tree with bokeh and biopython"
date:   2024-10-25 20:00:00
categories: bioinformatics
tags: [python,plotting]
thumbnail: /img/bokeh_tree.png
---

## Background

Bokeh can be used to make a wide variety of interactive plots. Here we use it to plot an interactive phylogenetic tree. We skip the more complicated part of tree drawing by using the Biopython `Phylo` module to read a newick tree and then get the tip coordinates from the tree object. From here it's a simple matter of connecting the tips with the branches and adding labels and so on. This code accepts a dataframe with metadata for the tip colors.

## Code

In the code below the function let's biopython do the tip arrangement work. The function `calculate_positions` gets the x and y positions for the tips. The branch segments are built using this information and vertical and horizontal lines drawn. Finally the tips are drawn using `scatter`. If a dataframe is provided this is merged with the base metadata derived from the tree tips. The `ColumnDataSource` is made from this. The toolips are built from the metadata columns but this could be customised.

```python
import numpy as np
import pandas as pd
from bokeh.plotting import figure
from bokeh.models import (ColumnDataSource,
                            HoverTool, BoxZoomTool, TapTool,
                            Legend, LegendItem, LabelSet,
                            Arrow, NormalHead, OpenHead, VeeHead)
from bokeh.io import output_notebook, output_file
output_notebook()
from bokeh.io import show, save
from Bio import Phylo
from io import StringIO
import pylab as plt
import matplotlib as mpl

def plot_phylogeny(tree, df=None, tip_size=10, lw=1, font_size='10pt',
                   tip_labels=True, label_offset=10, **kwargs):
    """
    Plots a phylogenetic tree with tips colored according to metadata.

    Parameters:
    - tree (Bio.Phylo.BaseTree.Tree): Phylogenetic tree from Biopython
    - df (pd.DataFrame): DataFrame containing metadata for each tip with a column for coloring
    - tip_size: size of tips
    - kwargs are passed to figure init
    """

    def calculate_positions(tree):
        x_positions = {}
        y_positions = {}

        def assign_y_positions(clade, y_base=0):
            if clade.is_terminal():
                y_positions[clade] = y_base
                return y_base + 1
            else:
                child_positions = []
                for child in clade.clades:
                    y_base = assign_y_positions(child, y_base)
                    child_positions.append(y_positions[child])
                y_positions[clade] = sum(child_positions) / len(child_positions)
                return y_base

        def assign_x_positions(clade, x_position=0):
            x_positions[clade] = x_position
            for child in clade.clades:
                branch_length = child.branch_length if child.branch_length else 0
                assign_x_positions(child, x_position + branch_length)

        assign_y_positions(tree.root)
        assign_x_positions(tree.root)
        return x_positions, y_positions

    tip_names = [clade.name for clade in tree.get_terminals()]
    # Get positions
    x_positions, y_positions = calculate_positions(tree)

    # Data for Bokeh
    x_h0, y_h0, x_h1, y_h1 = [], [], [], []
    x_v0, y_v0, x_v1, y_v1 = [], [], [], []

    for clade in tree.find_clades(order='level'):
        if clade.clades:
            children_y = [y_positions[child] for child in clade.clades]
            x_v0.append(x_positions[clade])
            x_v1.append(x_positions[clade])
            y_v0.append(min(children_y))
            y_v1.append(max(children_y))

            for child in clade.clades:
                x_h0.append(x_positions[clade])
                x_h1.append(x_positions[child])
                y_h0.append(y_positions[child])
                y_h1.append(y_positions[child])

    # Define ranges for plot
    x_min, x_max = min(x_positions.values()), max(x_positions.values())
    y_min, y_max = min(y_positions.values()), max(y_positions.values())
    x_padding = (x_max - x_min) * 0.1
    y_padding = (y_max - y_min) * 0.05
    x_range = (x_min - x_padding, x_max + x_padding)
    y_range = (y_min - y_padding, y_max + y_padding)
    if len(tip_names)>80:
        tip_labels = False
    if tip_labels == True:
        x_range = (x_range[0], x_range[1]*1.2)

    p = figure(tools="pan,ywheel_zoom,box_zoom,reset,save",
               x_axis_label="Branch length", y_axis_label="Clade",
               x_range=x_range,
               y_range=y_range, **kwargs)
    # Add segments for branches
    p.segment(x_h0, y_h0, x_h1, y_h1, line_width=lw, line_color="black")
    p.segment(x_v0, y_v0, x_v1, y_v1, line_width=lw, line_color="black")

    # Add terminal nodes with color mapping
    tip_coords = [(x_positions[clade], y_positions[clade]) for clade in tree.get_terminals()]
    x_tip, y_tip = zip(*tip_coords)

    #make metadata
    metadata = pd.DataFrame(zip(tip_names,x_tip,y_tip),columns=['name','x','y'])
    metadata['marker'] = 'circle'
    if df is not None:
        metadata = metadata.merge(df, left_on='name', right_index=True)
    if not 'color' in metadata:
        metadata['color'] = 'black'

    source = ColumnDataSource(metadata)
    #draw tips
    r = p.scatter('x', 'y', source=source, size=tip_size, color='color',
                  marker='marker')
    # Add tip labels
    if tip_labels == True:
        labels = LabelSet(x='x', y='y', text='name', level='glyph', text_color='color',
                        text_baseline='middle', x_offset=label_offset, text_font_size=font_size,
                        source=source)
        p.add_layout(labels)

    #hover tool
    tooltips=[]
    for c in metadata.columns:
        tooltips.append([c,f'@{c}'])
    h = HoverTool(renderers=[r], tooltips=(tooltips))
    p.add_tools(h)
    p.yaxis.visible = False
    p.xgrid.visible = False
    p.ygrid.visible = False
    p.toolbar.logo = None
    return p

 def get_color_mapping(df, col, cmap='Set1', seed=12):
    """Get color map for categorical dataframe column"""

    import matplotlib.colors as colors
    c = df.sort_values(col)[col].unique()
    c_map = mpl.colormaps[cmap]
    clrs = [colors.rgb2hex(c_map(i)) for i in range(len(c))]
    colormap = dict(zip(c, clrs))
    newcolors =  [colormap[i] if i in colormap else 'Black' for i in df[col]]
    return newcolors, colormap
```

## Usage

Here is a simple example. We read a newick string into a Biopython tree object and that is provided to the function. A color column is assigned from the dataframe using a color mapping function. If there is no color column the value is set to black.

```python
treedata = "((A:0.5,B:0.6),(C:0.7,(D:0.4,E:0.5):0.6),((F:0.8,G:0.9):0.4,H:0.7):0.5);"
# Sample metadata DataFrame with color groups
df = pd.DataFrame({
    'name': ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
    'group': ['Group1', 'Group1', 'Group2', 'Group2', 'Group3', 'Group3', 'Group4', 'Group4']
})
df = df.set_index('name',drop=False)
df['color'],c = get_color_mapping(df, 'group', 'Set1')
handle = StringIO(treedata)
tree = Phylo.read(handle, "newick")

p = plot_phylogeny(tree, df, font_size='14pt', lw=2, width=600, height=500)
show(p)
```

Here is the result saved as standalone html:

{% include bokeh_tree_plot.html %}
<br>
Another example with more tips but no coloring:

{% include bokeh_tree_plot2.html %}

## Links

[Bokeh](https://bokeh.org/)