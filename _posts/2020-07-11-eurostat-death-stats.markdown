---
layout: post
title:  "Eurostat deaths from all causes data"
date:   2020-07-11 12:18:00
categories: plotting
tags: [eurostat,seaborn]
thumbnail: /img/eurostat_scr.png
---

## Background

<div style="width: 250px; float:right;">
  <img src="/img/eurostat_scr.png" width="220px">
</div>

The number of deaths by week in most countries of the European Union are stored by [Eurostat](https://ec.europa.eu/eurostat) and accessible via  their data explorer tool. This is a bit of a clunky interface but fairly straightforward to use. The data on total deaths is useful for assessing the possible effects of the COVID-19 pandemic on the European population deaths. In particular excess deaths can be estimated by comparing to the mean values over a given period. The data used here was downloaded for all years >2005 with age and sex breakdown for all available countries. The United Kingdom has limited data.

## Flu season deaths create a regular cycle

The [flu season](https://en.wikipedia.org/wiki/Flu_season) accounts for a substantial amount of the cyclical changes in deaths in the Northern and Southern hemisphere. The cycle is regular and predictable. Influenza can be more severe depending on the circulating major influenzavirus subtypes. The exact mechanism behind the seasonal nature of influenza outbreaks is unknown.

<div style="width: auto; float:center;">
 <a href="/img/eurostat_flu_cycle.png"> <img class="scaled" src="/img/eurostat_flu_cycle.png"></a>
</div>

## Totals per year by age group during 'covid peak'

If we take the covid peak (see below) as being around weeks 10-20 of 2020, we can roughly compare the deaths for the same period each year. Remember these are total deaths, not for specific causes. This clearly shows the majority of the increase over other years or excess is due to the >70 age group.

<div style="width: auto; float:center;">
 <a href="/img/eurostat_fluseason_deaths.png"> <img class="scaled" src="/img/eurostat_fluseason_deaths.png"></a>
</div>

## Total deaths per year up to end of June

This plot shows the totals for the first half of each year, up to June for selected countries.

<div style="width: auto; float:center;">
 <a href="/img/eurostat_4countries_totaldeaths.png"> <img class="scaled" src="/img/eurostat_4countries_totaldeaths.png"></a>
</div>

## Total deaths per year up to end of June, individual plots

The same plot broken down for multiple countries. Italy has incomplete data for 2020 so seems lower.

<div style="width: auto; float:center;">
 <a href="/img/eurostat_totaldeaths_bycountry.png"> <img class="scaled" src="/img/eurostat_totaldeaths_bycountry.png"></a>
</div>

## 2020 weekly trend compared to mean

To compare the actual peaks during the worst covid period we simply plot the trend line for each country with a mean line for all years for reference. This kind of plot is produced by [Euromomo](https://www.euromomo.eu/graphs-and-maps). This shows the clear peak of deaths in older age groups from March.

<div style="width: auto; float:center;">
 <a href="/img/eurostat_2020peak_trend.png"> <img class="scaled" src="/img/eurostat_2020peak_trend.png"></a>
</div>

## 2020 weekly trend comparison by age

How do the above trends break down by age group? This is seen below for Sweden with 2020 in red and 2018 in orange.

<div style="width: auto; float:center;">
 <a href="/img/eurostat_2020peak_trend_byage_sweden.png"> <img class="scaled" src="/img/eurostat_2020peak_trend_byage_sweden.png"></a>
</div>

However for Spain the peak in deaths is seen for all older age groups even at ages 50-59, though the total numbers are much smaller compared to the >70 group.

<div style="width: auto; float:center;">
 <a href="/img/eurostat_2020peak_trend_byage_spain.png"> <img class="scaled" src="/img/eurostat_2020peak_trend_byage_spain.png"></a>
</div>

Note that the above plots do not use any covid related dataset but purely data on total deaths regardless of cause.

These plots were made using pandas and seaborn. The code can be found in a Jupyter notebook [here](https://github.com/dmnfarrell/teaching/tree/master/covid_stats)

## Links

* [source of data](https://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=demo_r_mweek3&lang=en)
* [About Eurostat Weekly death statistics](https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Weekly_death_statistics)
