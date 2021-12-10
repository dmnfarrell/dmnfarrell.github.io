---
layout: post
title:  "Excess deaths in Ireland for Summer 2021"
date:   2021-10-20 09:08:00
categories: general
tags: [plotting,rip.ie,python]
thumbnail: /img/ireland_deaths_ripie_summary_v3.png
---

## Background

In a [previous post](/plotting/ireland-deaths-reanalysis) I showed how to retrieve daily notices from RIP.ie. These can be used to count the number of daily deaths and compare between years. This is a useful proxy for the official registered deaths and provides an up to date estimate of current trends since registered data is delayed by up to 3 months. Deaths in Ireland have shown an unusual increase over the Summer months of 2021. This is higher for the months of June to September than any previous year. These are raw figures and not adjusted for population increase or age profile, but that would make relatively little difference for the past few years. **The cause of this increase is unknown** yet it should surely be a concern for public health systems. Unlike the UK, Ireland has poor reporting on live mortality statistics. RIP.ie was not designed for this task and therefore we cannot further analyse this data by age group for example.

Note that the method used to extract RIP.ie data here attempts to clean the data of duplicates and extract correct dates. This process may not be perfect so the figures are close approximations. Also note that is possible this increase in 2021 months is an artifact such as the result of changes in how deaths are recorded to RIP.ie, though this seems unlikely.

## Year on year comparison for summer months

<div style="width: auto;">
 <a href="/img/ireland_deaths_ripie_summary_v3.png"> <img class="small-scaled" src="/img/ireland_deaths_ripie_summary_v3.png"></a>  
   <p class="caption">Yearly totals compared for months of June to September.</p>
</div>

This increase is not reflected in the current Euromomo data which uses the official registered death figures which have a significant lag in updating.

<div style="width: auto;">
 <a href="/img/euromomo_ireland_2021.png"> <img class="scaled" src="/img/euromomo_ireland_2021.png"></a>  
   <p class="caption">Euromomo deaths for Ireland 2018-present.</p>
</div>

## Weekly values

The following shows the deaths by week. This data is approximately up to date to week 41 (Mid-October) of 2021, though there are probably some missing entries for the latest week.

<div style="width: auto;">
 <a href="/img/ireland_deaths_ripie_byweek.png"> <img class="scaled" src="/img/ireland_deaths_ripie_byweek.png"></a>  
   <p class="caption">Weekly death totals compared from 2018-2021.</p>
</div>

## Eurostat using RIP.ie?

It has been pointed out by Seamus Coffey [here](https://twitter.com/seamuscoffey/status/1450537528063406099) that Ireland has recently been added to the [Eurostat](https://ec.europa.eu/eurostat/cache/statistics_explained/visualisations/weeklydeath/?lang=en) deaths data series and they appear to be using the same data source as indicated here:

<div style="width: auto;">
 <a href="/img/eurostat_ireland_deaths_twitter.jpg"> <img class="scaled" src="/img/eurostat_ireland_deaths_twitter.jpg"></a>  
   <p class="caption">Eurostat Ireland weekly deaths compared to RIP.ie (posted by Seamus Coffey).</p>
</div>

I can therefore also compare my estimates with the Eurostat data as below. As you can see there is an overestimate in my data by a consistent value of about 50. I don't know why this is. However relative values remain the same between years so the conclusions stand that there are increases in 2021 Summer months over other years.

<div style="width: auto;">
 <a href="/img/eurostat_ireland_deaths_compared.png"> <img class="scaled" src="/img/eurostat_ireland_deaths_compared.png"></a>  
   <p class="caption">Eurostat Ireland weekly deaths compared to our RIP.ie estimate.</p>
</div>

## Links

* [RIP.ie](https://rip.ie/)
* [Eurostat deaths](https://ec.europa.eu/eurostat/cache/statistics_explained/visualisations/weeklydeath/?lang=en)
* [GRO](https://www.gov.ie/en/service/49c66f-registering-a-death-in-ireland/)
* [Our World in Data dataset](https://covid.ourworldindata.org/data/owid-covid-data.csv)
