---
layout: post
title:  "Categorical plotting with the factor plots plugin"
date:   2015-12-05 13:05:00
categories: dataexplore
---

<iframe width="320" height="240" style="float: right; padding:12px;"
src="https://www.youtube.com/embed/DyHO7XfBW4o"
frameborder="0" allowfullscreen>
</iframe>

##Background

Sometimes you might have data with various categorical columns that you want to summarise. This will often be to view statistical information like the spread of values in a column. This can be done using histograms or boxplots in dataexplore. However the *factor plots plugin* provides more advanced  statistical graphics when you want to see variability over multiple related categories. Seaborn is a Python visualization library based on matplotlib. Since the functions operate on dataframes they are ideally suited for use in dataexplore.

##Factor plotting

Factor plots allow multiple comparisons to be made in a single graph. That is, you can split data by more than one variable along an axis or between plots. In seaborn these dimensions are called row/col (the plot dimensions) x,y (axes) and hue (grouping/color within plots). These concepts are illustrated in the [documentation on factor plotting](http://stanford.edu/~mwaskom/software/seaborn/generated/seaborn.factorplot.html) and below.

Although this can be a very useful tool, be aware of the caveat given in the seaborn documentation: "although it is possible to make rather complex plots using this function, in many cases you may be better served by created several smaller and more focused plots than by trying to stuff many comparisons into one figure". Plots like this get quickly complicated and may often confuse the viewer.

##Usage

<div style="width: 400px;">
<img src="/img/factorplot_dialog.png" width="500px">
</div>

Open the plugin from the plugins menu using plugins->factor plots. The dialog appears below the main window. The controls are mostly obvious apart from the 'factor' widgets. The combo boxes allow you to select which columns are used as a 'factor' in the plot. Usually you should select at least the col value. This will split the data by that variable (column) into subplots with n columns (use the column wrap option to set this). Hue will split the data in one plot and color by this variable with a legend.

##Examples

You can quickly try out the plotting function using the sample dataset. This is simply some random columns with a categorical 'label' variable. To plot distributions of values over this variable in multiple plots simple select 'label' in the col menu then re-plot. You will see a bar graph with error bars representing the spread in each group. Try the other kinds of plot like box, violin etc. You will get plots like those shown below.

<div style="width: 400px;">
<img src="/img/sample_factorplot_formats.png" width="800px">
</div>

Try using 'label' in the hue and x options instead (delete the other factor selections or you will get weird results) and you will see the plots arranged as below. On the left the the label column as hue, so split in a single plot into groups. On the right is the result when 'label' is x. This makes each bar represent all data for label a, b etc.. Note that when the x value is not selected the data is always plotted on the x-axis by the columns remaining. So _you should limit your table selection to only include the columns you want in the plot_.

<div style="width: 400px;">
<img src="/img/sample_factorplot_factors.png" width="600px">
</div>

##Tips data

You can see from the above that such plots get complicated quickly and picking the correct 'factor' to plot may require some thought. As another, more realistic example, the tips data available in DataExplore can also be used to test the plugin. This has more categorical variables that can be factors and hence even more complex. The plot below is split on three factors - sex, smoker and day. You may try a combination of selections in the dialog to see which ones correspond to this plot.

<div style="width: 400px;">
<img src="/img/factorplot_tips1.png" width="400px">
</div>

## Links

* [Seaborn](http://stanford.edu/~mwaskom/software/seaborn/)

* [Factor analysis](https://en.wikipedia.org/wiki/Factor_analysis)