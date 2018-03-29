---
layout: post
title:  "New grid layout dialog for sub-plots"
date:   2018-03-28 11:20:20
categories: dataexplore
tags: [plotting]
---

Plotting in grid layouts is now easier with a new grid layout tool using an interactive table that you can lay out your plots in.

## How to use

The gif animation below shows how to use the tool to generate subplots by clicking and dragging in the grid to select the area for your next plot. Note that subplots will be overwritten if you select the same cell as one currently occupied but if you drag over this cell the region will be plotted over. The tool assumes the user will know how to avoid overlaps. So it's best to have a good idea of how to layout the plots beforehand, or just use trial and error. You can remove subplots from the drop down menu, listed according to their positions.

<div style="width: 400px;">
<img src="/img/grid_layout_example.gif" width="600px">
</div>

## Automatic multiple views

Grid layout includes a 'multiview' feature. This allows you to auto-generate different kinds of plot in the grid for the same data every time you plot. This could be useful for quickly previewing regions of data repeatedly without having to set the plot type each time. This will overwrite whatever plot you currently have displayed. The feature is also illustrated in the gif above.
