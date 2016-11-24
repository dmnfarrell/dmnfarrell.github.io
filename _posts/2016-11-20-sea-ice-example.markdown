---
layout: post
title:  "Example: importing and plotting arctic sea ice extent data"
date:   2016-11-20 11:47:00
categories: dataexplore
tags: [plotting,arctic]
---

<iframe width="320" height="240" style="float: right; padding:12px;"
src="https://www.youtube.com/embed/SHyADz4Rm0A"
frameborder="0" allowfullscreen>
</iframe>

## Background

Version 0.73 adds a few minor features and fixes that are illustrated here using data from NOAA/NSIDC on arctic sea ice extent. The accompanying screencast goes through the importing and plotting steps. The data used here was calculated using raw data provided by NSIDC (the National Snow and Ice Data Center). This combines measures of north and south sea ice extent to estimate the global total. It has been in the news due to the dramatic anomaly in the winter 2016 trend of ice formation. Explanations of the data can be found in the links below.
The purpose of this example is to get the data into dataexplore and reproduce the plot below:

<div style="width: 700px;">
<img src="https://sites.google.com/site/arctischepinguin/home/sea-ice-extent-area/grf/nsidc_global_extent_byyear_b.png" width="500px">
</div>
(Image taken from the artic sea ice forum created by user Wipneus). Note that dataexplore can't currently plot error bands as in this graph.

## Links

* [2016 sea ice area and extent forum](https://forum.arctic-sea-ice.net/index.php/topic,1457.msg93338.html)
* [weather underground blog post](https://www.wunderground.com/blog/JeffMasters/crazy-cryosphere-record-low-sea-ice-an-overheated-arctic-and-a-snow)
* [ftp site with original data](ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/south/daily/data/)
