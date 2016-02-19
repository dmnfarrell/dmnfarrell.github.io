---
layout: post
title:  "Combining different plots in a grid layout"
date:   2016-01-16 13:07:00
categories: dataexplore
tags: [plotting]
---

## Background

The default plotting mode in DataExplore is to continually overwrite the plot window based on the currently plotted data. Some of the plotting options will generate multiple subplots but these are re-created in one step. For some applications you will want to combine different ways of looking at the data in one plot or plots from different tables.

## Usage

This can be done by selecting the 'use grid layout' option in the plot viewer and then using the Grid Layout tab dialog (shown below) to adjust the rows and column settings to place the subplots. Every time you re-plot the data will be plotted at the current row/column location, clearing the previous subplot at that location. You will need to clear the current figure before starting unless you want to overlay on this.

<div style="width: 400px;">
<a href="/img/plot_layout_options.png"><img src="/img/plot_layout_options.png" width="500px"></a>
</div>

Changing the number of rows/columns makes a finer grid then adjusting the row/column spans allow fitting varying sizes. Below the example plot uses a grid of 3 columns and 2 rows. The title on each plot has been added to show it's row/col co-ordinate. You can remove subplots using the dialog by selecting them in the plot list. Plots can only be changed by re-selecting the data used to create them.

<div style="width: 400px;">
<a href="/img/grid_layout.png"><img src="/img/grid_layout.png" width="500px"></a>
</div>

It is also possible to use the grid to make an inset plot. You simply plot as normal, then don't clear the plot but just select 'use grid layout' and add your subplot(s) where you want. The example below uses an inset plot at row 1 and column 3.

<div style="width: 400px;">
<a href="/img/inset_plot.png"><img src="/img/inset_plot.png" width="500px"></a>
</div>

If you want to combine plots from different tables the current way is to create a sub-table and paste the other table into it. Then plot onto the grid from this table.

## Links

* [Matplotlib Subplot Using GridSpec](http://matplotlib.org/users/gridspec.html)
