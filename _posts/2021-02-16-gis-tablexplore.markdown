---
layout: post
title:  "A simple GIS plugin for Tablexplore "
date:   2021-02-16 10:18:00
categories: software
tags: [maps,tablexplore]
thumbnail: /img/tablexplore_gis.png
---

## Background

<div style="width: 320px; float:right;">
<img src="/img/tablexplore_gis.png" width="280px">
</div>

[Tablexplore](https://github.com/dmnfarrell/tablexplore) is an application for data analysis and plotting built in Python using the PySide2/Qt toolkit. The interface allows quick visualization of data with convenient plotting. The program is intended mainly for educational/scientific use but should be useful for a variety of general applications. The program also has a plugin system that allows extra functionality to be added arbitrarily. These could be run in separate windows or placed in a pane inside the main application, under the table. From the plugin you could use the table or plotting tools or just make a plugin that does a utility function like file manipulation.

## A small GIS plugin

As an example this post describes a GIS plugin. [GIS](https://www.esri.com/en-us/what-is-gis/overview) tools are for analysis in geographic information systems. ArcGIS (commercial) and QGIS (free) are two such programs. Here we implement some of the more elementary aspects of GIS in a plugin that can use the plotting and table frames in Tablexplore. The interface consists only of a small dialog underneath the main table. Here maps are loaded from shapefiles into the list view. Each map layer can be plotted, moved up and down and properties like colors, labels changed. Layers can be overlaid or plotted in different subplots. Currently you can only import shapefiles which are sets of coordinates for points or polygons. These can be edited in the table and new columns added for coloring purposes for example.

<div style="width: auto;">
 <a href="/img/tablexplore_gis_usage.gif"> <img class="small-scaled" src="/img/tablexplore_gis_usage.gif"></a>
</div>

## Tools

Some basic tools are available for geometric manipulation. A shown below we can perform some typical transform operations like rotate, buffer and convex hull. Each result is put into a new layer. The Simulate Shapes option lets you add test polygons for practice or teaching use.

<div style="width: auto;">
 <a href="/img/tablexplore_gis_tools.gif"> <img class="small-scaled" src="/img/tablexplore_gis_tools.gif"></a>
</div>

Geopandas is used for this plugin so if you are using tablexplore from the pip installation you should install this library yourself using `pip install geopandas`. If using the Windows executable, AppImage or snap then the plugin should be available without you needing to do anything extra.
Plugin code examples are available [here](https://github.com/dmnfarrell/tablexplore/tree/master/tablexplore/plugins).

## Links

* [github page](https://github.com/dmnfarrell/tablexplore)
* [homepage](https://dmnfarrell.github.io/tablexplore/)
* [QGIS - A Free and Open Source Geographic Information System](https://www.qgis.org/en/site/)
