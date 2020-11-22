---
layout: post
title:  "Epidemics, PCR and the dangers of mass testing"
date:   2020-11-21 12:06:00
categories: general
tags: [covid19,sarscov2]
thumbnail: /img/dna_amp.png
---

## Background

<div style="width: 400px; float:right;">
<img src="/img/dna_amp.png" width="300px">
</div>

Very early on when the SARS-CoV-2 virus was detected in China, PCR tests were developed for detecting the presence of viral RNA in sputum samples. [PCR](https://en.wikipedia.org/wiki/Polymerase_chain_reaction) stands for Polymerase chain reaction. The method used in this case is RT-PCR becuase the viral RNA has to be converted into DNA first. It is a technique very widely used in molecular biology to detect sometimes tiny amounts of DNA isolated from a biological sample. It does this by amplifying a piece of DNA using a template that you know exists uniquely in the virus. The amplification is repeated in successive cycles until there is enough DNA. The more cycles used, the more likely it will detect a small fragment. If you use enough cycles eventually you will get a signal, but it won't be reliable beyond a certain point. So a cutoff point is selected where you stop amplifying. You can also quantify the amount of DNA present by counting the cycles needed to get a specific signal, the lower the cycles needed the more DNA is present. This equates approximately to viral load.

## False positives

The sensitivity of the PCR method can be very useful but it is also a liability when it comes to mass testing. False positives occur for a few different reasons. Even if the test is perfectly done, they will happen a small percentage of the time. This is the false positive rate (FPR). The exact value isn't known and depends on the specifics of the test. However it could be from 0.8% to 4% according to [the Lancet](https://www.thelancet.com/journals/lanres/article/PIIS2213-2600(20)30453-7/fulltext). The FPR in practice will also depend on the threshold of cycles used as mentioned above. You can change your sensitivity by adjusting this threshold and that makes it have less false positives but might miss some true ones. So there is a trade-off always in testing. The threshold is somewhat arbitrary.

There are additional source of error though. If someone has had the virus some time previously it may amplify traces of the virus that are left in the upper respiratory tract even though no live virus is present. Another important source is cross contamination, where another sample that is positive contaminates the negative one during sample collection or processing. Then both samples are positive. This problem would be made worse in difficult conditions were issues with handling, hygiene and overworked staff might all combine to make it worse. So the effective FPR is a combination of things and will probably go up in mass testing situations.

## The problem

Testing has been carried out in the wealthier countries on an industrial scale since the Summer. The UK is now carrying out about 300,000 tests per day. It's important to realise that the total false positives will be a percentage (the FPR) of total tests done. The more tests you do, the more false positives. Below is a simplistic illustration of the effect of false positives with large scale testing. The orange is the 'true' signal and represents some portion of actual detectable cases. The red line is roughly what would be seen from tests. As testing, in green, is increased it can capture the cases better. With more testing early on you get a better reflection of the true epidemic. However as testing is increased but prevalence goes down you will get a significant proportion of false detections thus inflating the figures. These four cases are the same except there is more testing being done in total in each example. The first plot we could consider as analagous to a more targeted testing regime were few tests are being done but on symptomatic patients. False positives are not a big problem here. In the final plot 20 times more testing is being done and at late stages we are doing so many that it will be now in situations were people are not symptomatic. We have now altered the so-called 'pretest probability' of disease and could be virtually sampling randomly. When low pretest probability exists, positive results should be [interpreted with caution](https://pubmed.ncbi.nlm.nih.gov/32398230/). These plots assume an FPR of **0.5%** in all cases. Obviously reality is much more complex than this so this is just a crude example. But if the FPR is higher than we think it will make the problem even worse.

<div style="width: auto; float:center;">
 <a href="/img/fp_problem_example.jpg"> <img class="scaled" src="/img/fp_problem_example.jpg"></a>
 <p class="caption">Illustration of the potential problem with false positives in the presence of large scale testing in a low prevalence regime. This is purely hypothetical. </p>
</div>

## The virus is probably endemic and seasonal

No one is claiming that all the positives seen in testing results are false. The virus is clearly seasonal and has made a return in the Autumn. However, since operational FPR varies in practice it's not possible to know the true signal in the noise when we do so much testing. In particular it would be important to know the statistics on 'weak' positives. There is [no convincing data](https://www.ecdc.europa.eu/en/covid-19/latest-evidence/transmission) to suggest that detection of low levels of viral RNA by PCR equates with infectivity unless confirmed by viral culture. Despite Covid deaths rising excess deaths for this Winter thus far look the same as 2018. Does this mean that such deaths are being erroneously labelled?

There is a very real possibility that mass testing with contact tracing has led to an over estimate of true cases of infection. This then feeds back onto itself and drives more testing because every positive is treated as a case. Cases (positives) appear to drive much of government policy, modelling and public support for the measures. Since the restrictions have such enormous effects on society we must look very seriously at this possibility and act accordingly. New tests that do not suffer from such problems like the antigen tests could be considered to replace PCR for large scale use. In any case, testing of asymptomatic individuals is not advisable unless the test is extremely accurate and the policy of mass testing itself should be questioned.

## Links

* [Dr. Clare Craig on COVID testing problems](https://www.youtube.com/watch?v=380DLg-nAqI)
* [Ch4 Dispatches report on a UK COVID-Testing Lab](https://www.channel4.com/press/news/dispatches-uncovers-serious-failings-one-uks-largest-covid-testing-labs)
* [Impact of false-positives and false-negatives in the UK's COVID-19 RT-PCR testing programme](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/895843/S0519_Impact_of_false_positives_and_negatives.pdf)
* [Dartmouth-Hitchcock: Epidemic That Wasn’t](https://www.nytimes.com/2007/01/22/health/22whoop.html)
* [Ireland's COVID-19 Data Hub](https://covid19ireland-geohive.hub.arcgis.com/)
* [HPSC report: Proposal for the management of weak positive PCR results](https://www.hpsc.ie/a-z/respiratory/coronavirus/novelcoronavirus/guidance/outbreakmanagementguidance/PCR%20weak%20results%20guidance.pdf)
