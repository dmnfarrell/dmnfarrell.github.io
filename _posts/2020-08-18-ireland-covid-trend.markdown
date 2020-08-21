---
layout: post
title:  "Ireland COVID-19 trend in positive rate"
date:   2020-08-18 10:26:00
categories: plotting
tags: [covid19]
thumbnail: /img/ireland_covid_tests_prate.png
---

## Background

Irelands Covid-19 2020 outbreak has led to 1,775 official total deaths (HPSC data) with 27,499 positive tests. The first case identified was on 29 February and first identified fatality was on 11 March. Since then, like most European countries, Ireland has been in a delaying phase  with varying degrees of restrictions. Some of these have been lifted gradually but the policy remains one of suppression. Therefore any spikes in detected cases are met with alarm. Increases in test positives in August were largely caused by localised clusters, some located in meat factories. In response, on August 7 the government locked down three counties containing these clusters. On 18 August, the government announced [six new measures](https://www.thejournal.ie/government-advice-gatherings-5178822-Aug2020/) that again restrict gatherings. These responses are centred on a policy of complete suppression which will have to continue indefinitely according to the current logic, as long as the virus circulates and is tested for.

The data used here is provided at [data.gov.ie](https://data.gov.ie/dataset?q=covid&sort=score+desc%2C+metadata_created+desc). The plots below can be updated at any time by downloading the notebook and running it locally.

## Tests performed vs positive rates

With any ongoing transmission of a virus, continual community testing will eventually reveal some clusters. Many will go likely go undetected and may not spread if sufficient immunity exists or if distancing measures prevent it (or both). If we plot the tests performed and actual positives we can see the rate (positives/total tests) has steadily declined since May. (The data is plotted here as a rolling mean with 7 day window to smooth the curves). This is despite the gradual re-opening. Tests have increased rapidly again in August in response to the clusters with contact tracing allowing other cases to be discovered. Thus by definition more cases are found in these areas. Such fluctuations are inevitable and a result of testing and tracing. However the total positive rate continues downward.

<div style="width: auto; float:center;">
 <a href="/img/ireland_covid_tests_prate.png"> <img class="scaled" src="/img/ireland_covid_tests_prate.png"></a>
</div>

The by-county nature of the new cases can be seen in these plots (Covid19CountyStatisticsHPSCIreland dataset). The three counties of interest show spikes but whether they are in clusters (they largely were) can't be seen in these plots.

<div style="width: auto; float:center;">
 <a href="/img/ireland_covid_tests_bycounty.png"> <img class="scaled" src="/img/ireland_covid_tests_bycounty.png"></a>
</div>

## How to respond to spikes in test positives?

The reporting of cases due to positive tests without reference to localisation, symptomatic status, age group and testing numbers is not very helpful. HPSC does produce [epidemiology reports](https://www.hpsc.ie/a-z/respiratory/coronavirus/novelcoronavirus/casesinireland/epidemiologyofcovid-19inireland) which have some more detail but few will read them. Unfortunately, details matter. For example, the majority of cases in August have been in people under the age of 45. Many cases have been due to obvious clusters probably where repeated transmissions occur within groups as in the meat plants. The public health response appears bizarre in this context. If these known locations can be monitored and restrictions placed there in the short term, what is the purpose of regional 'lockdowns'?

## Links

* [notebook on github](https://github.com/dmnfarrell/teaching/blob/master/covid_stats/ireland.ipynb)
* [data.gov.ie covid datasets](https://data.gov.ie/dataset?q=covid&sort=score+desc%2C+metadata_created+desc)
* [hpsc epidemiology reports](https://www.hpsc.ie/a-z/respiratory/coronavirus/novelcoronavirus/casesinireland/epidemiologyofcovid-19inireland/)
* [Ireland covid-19 data hub](https://covid19ireland-geohive.hub.arcgis.com)
