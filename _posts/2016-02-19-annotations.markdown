---
layout: post
title:  "Adding annotations to plots"
date:   2016-02-19 15:35:00
categories: dataexplore
tags: [annotation]
---

Version 0.7.2 of DataExplore adds the ability to create text annotations on plots. This is still an experimental feature but the basic functionality is to allow users to add text labels to plots. This may be useful for publication or presentation purposes. Normally the idea would be to add the annotations when you have prepared the final plot to your satisfaction. However the labels will be saved and can be moved around the plot by dragging the the mouse, so you can add labels and adjust them later. Individual labels can also be deleted.

##Adding a label

In the annotations tab the options are generally self explanatory. Font size, colors,  box shape etc. can be selected. Note the annotations can't be formatted or edited 'in place' (except to move them) so pick the options you want first and press 'create'. If you don't like it just delete and then add another. The 'clear' button deletes all current labels.

##Coordinate systems

Labels can be added according to three different coordinate systems:

1. their position relative to the data points ('data')
2. their position in the current subplot ('axes fraction')
3. their position in the whole figure ('figure fraction')

Adding using data coordinates (default) allows you to label specific points and the label will move with the data if the plot is redrawn with more or less data. The label will not appear if it is outside the range of the current plot but will be remembered. If you want to keep the label constant in the plot use the other two options. 'axes fraction' means the label x,y position is stored as a fraction of the current axis or subplot. This is useful to label subplots in the grid layout. The last option keeps the position the same regardless as it is in the reference frame of the whole figure.

##Example

In the plot below, labels A and B are added as data labels.

<div style="width: 400px;">
<a href="/img/annotations1.png">
<img src="/img/annotations1.png" width="400px"></a>
</div>
When the plot is redone with more points the x and y scales change and the labels shift also and stay next to the data points. An axis and figure label have also been added here. These will not move when replotting with new data.

<div style="width: 400px;">
<a href="/img/annotations2.png">
<img src="/img/annotations2.png" width="400px"></a>
</div>
Finally, if we clear the plot and use the grid layout, the same plot as above is placed as a subplot. The data labels are in the correct place. The axis label tracks the subplot and the figure label remains in the upper right.

<div style="width: 400px;">
<a href="/img/annotations3.png">
<img src="/img/annotations3.png" width="400px"></a>
</div>

##Caveats

This feature is not fully complete and there may be some shifting of labels upon replotting. Therefore you should finalize your plot and adjust your labels with the mouse before saving the final result.

##Links

* [Annotating Axes in matplotlib](http://matplotlib.org/users/annotations_guide.html/)
