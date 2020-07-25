---
layout: post
title:  "Death causes in England and Wales comparison - Winton Centre"
date:   2020-07-24 11:36:00
categories: plotting
tags: [covid19,seaborn]
thumbnail: /img/winton_page_scr.png
---

## Background

The flood of information about Covid-19 deaths and cases in the media makes it sometimes very hard for people to get information in context. Any number can look big or small without a relative measure. As such it becomes very hard to get a realistic estimate of risk on the part of the public. The **Winton Centre for Risk and Evidence Communication** have [compiled an interesting table](https://wintoncentre.maths.cam.ac.uk/coronavirus/covid-19-resources-make-sense-numbers-test/how-have-covid-19-fatalities-compared-other-causes-death/) of causes of mortality in England and Wales. They have used this to compare typical deaths in the period up to July from a selection of causes to covid-19 estimated deaths. It is worth reading the text in the original page explaining how they made their data table. It is also pointed out that there is no completely 'neutral' presentation and that this data is merely showing comparative causes by age group. Though it is obvious that age is the major risk factor, it says nothing about possible long-term health effects of those who became ill. The plots below are made using the second table, deaths scaled to a typical 16 weeks.

## Comparisons of 16 weeks of Covid deaths vs other causes in England and Wales

Here the table is broken into three age groups across multiple 5 year bands to compress into one readable plot. The x-axis is a log scale for easier comparison.

<div style="width: auto; float:center;">
 <a href="/img/eng_wales_deaths_causes_all.png"> <img class="scaled" src="/img/eng_wales_deaths_causes_all.png"></a>
</div>

## Comparison by age group

<div style="width: auto; float:center;">
 <a href="/img/eng_wales_deaths_causes_byage.png"> <img class="scaled" src="/img/eng_wales_deaths_causes_byage.png"></a>
</div>

## The same data grouped by cause

<div style="width: auto; float:center;">
 <a href="/img/eng_wales_deaths_causes_bycause.png"> <img class="scaled" src="/img/eng_wales_deaths_causes_bycause.png"></a>
</div>

The jupyter notebook and csv file for making these plots can be found [here](https://github.com/dmnfarrell/teaching/tree/master/covid_stats/).

## Links

* [notebook on github](https://github.com/dmnfarrell/teaching/blob/master/covid_stats/uk.ipynb)
* [How have Covid-19 fatalities compared with other causes of death?](https://wintoncentre.maths.cam.ac.uk/coronavirus/covid-19-resources-make-sense-numbers-test/how-have-covid-19-fatalities-compared-other-causes-death/)
* [Our World in Data Coronavirus site](https://ourworldindata.org/coronavirus)
