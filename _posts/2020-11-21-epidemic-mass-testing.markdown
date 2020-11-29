---
layout: post
title:  "Epidemics, PCR and the dangers of mass testing"
date:   2020-11-21 12:06:00
categories: general
tags: [covid19,sarscov2]
thumbnail: /img/dna_amp.png
---

## Background

<div style="width: 320px; float:right;">
<img src="/img/dna_amp.png" width="300px">
</div>

Very early on when the SARS-CoV-2 virus was detected in China, PCR tests were developed for detecting the presence of viral RNA in sputum samples. This is the so-called swab test. [PCR](https://en.wikipedia.org/wiki/Polymerase_chain_reaction) stands for Polymerase chain reaction. The method used in this case is RT-PCR becuase the viral RNA has to be converted into DNA first. It is a technique very widely used in molecular biology to detect sometimes tiny amounts of DNA/RNA isolated from a biological sample. It does this by amplifying a piece of DNA using a template that you know exists uniquely in the target (i.e. the virus). The amplification is repeated in successive cycles until there is enough DNA to detect. The more cycles used, the more likely it will detect a small fragment. If you use enough cycles eventually you will get a signal, but it won't be reliable beyond a certain point. So a cutoff point is selected where you stop amplifying. You can also quantify the amount of DNA present by counting the cycles needed to get a specific signal, the lower the cycles needed the more DNA is present. This equates approximately to viral load.

## False positives

The sensitivity of the PCR method can be very useful but it is also a liability when it comes to mass testing. False positives occur for a few different reasons. Even if the test is perfectly done, they will happen a small percentage of the time. This is the false positive rate (FPR). The exact value isn't known and depends on the specifics of the test. You may see figures ranging from about 0.1-0.5% for this. There are additional source of error though that are context dependent. If someone has had the virus some time previously it may amplify traces of the virus that are left in the upper respiratory tract even though no live virus is present. Another important source is cross contamination, where another sample that is positive contaminates the negative one during sample collection or processing. Then both samples are positive. This problem would be made worse in difficult conditions were issues with handling, hygiene and overworked staff might all combine to make it worse. So the effective or operational FPR is a combination of things and will probably go up in mass testing situations. This actual FPR could be from 0.8% to 4% according to [the Lancet](https://www.thelancet.com/journals/lanres/article/PIIS2213-2600(20)30453-7/fulltext) though that is just an educated guess.

There is another sense in which 'false positives' are being talked about that have less to do with the technical details of the test but the interpretation - what constitutes a 'case'. It depends on the threshold of cycles used as mentioned above. You can change your sensitivity by adjusting this threshold and that produces less false positives but might miss some true ones. So there is a trade-off always in testing. The threshold is somewhat arbitrary. It is a policy decision. This last meaning is probably the one most people mean really when they talk about the subject and it implies all the factors mentioned above as well as the interpretation part. That is at least how I understand it.

## The problem

Testing has been carried out in the wealthier countries on an industrial scale since the Summer. The UK is now carrying out about 300,000 tests per day. It's important to realise that the total false positives will be a percentage (the FPR) of total tests done. The more tests you do, the more false positives. Below is a simplistic illustration of the effect of false positives with large scale testing. The orange is the 'true' signal and represents some portion of actual detectable cases. The red line is roughly what would be seen from tests. As testing, in green, is increased it can capture the cases better. With more testing early on you get a better reflection of the true epidemic. However as testing is increased but prevalence goes down you will get a significant proportion of false detections thus inflating the figures. These four cases are the same except there is more testing being done in total in each example. The first plot we could consider as analagous to a more targeted testing regime were few tests are being done but on symptomatic patients. False positives are not a big problem here. In the final plot 20 times more testing is being done and at late stages we are doing so many that it will be now in situations were people are not symptomatic. We have now altered the so-called 'pretest probability' of disease and could be virtually sampling randomly. When low pretest probability exists, positive results should be [interpreted with caution](https://pubmed.ncbi.nlm.nih.gov/32398230/). These plots assume an FPR of **0.5%** in all cases. But if the FPR is higher than we think it will make the problem even worse. Obviously reality is much more complex than this so this is just a crude example.

<div style="width: auto; float:center;">
 <a href="/img/fp_problem_example.jpg"> <img class="scaled" src="/img/fp_problem_example.jpg"></a>
 <p class="caption">Illustration of the potential problem with false positives in the presence of large scale testing in a low prevalence regime. This is purely hypothetical. </p>
</div>

## Actual FPR could vary greatly depending on the testing regime

The virologist Ian Mackay has discussed this subject in a [blogpost](https://virologydownunder.com/the-false-positive-pcr-problem-is-not-a-problem/). However the title of that article might be a bit misleading. He correctly states that the FPR is not a problem in specific contexts. In Australia the positive rates have been extremely low and this is said to be evidence of low FPR. Certainly it must be low in Australia. However, as he points out, prevalence is virtually zero there. This means two things: 1) sources of cross contamination in labs are very minimal as there are so few real positives and 2) there are next to no individuals with recent infections from which viral fragments could be falsely detected as positives. Either or both could be significant factors in other places like the UK where there is far more testing being done in labs that seem to have had problems with sampling and handling. These things alone could be enough to alter the actual FPR in the UK, though it will vary.

## The virus is probably endemic and very likely seasonal

No one is claiming that all the positives seen in testing results are false. The virus is clearly seasonal and has made a return in the Autumn all across Europe. However, since operational FPR varies in practice it's not possible to know the true signal in the noise when we do so much testing as in the UK. The PCR test as a means to detect virus and even diagnose infections is perfectly good. It is the interpretation of results that is important. In particular it could be useful to know the statistics on 'weak' positives. There is [no convincing data](https://www.ecdc.europa.eu/en/covid-19/latest-evidence/transmission) to suggest that detection of low levels of viral RNA by PCR equates with infectivity unless confirmed by viral culture. So counting weak positives (those at high cycle thresholds) as clinical cases is problematic. Despite Covid deaths rising excess deaths for this Winter thus far look the same as 2018 across most of Europe. Does this mean that such deaths are being erroneously labelled?

At the start of the pandemic using PCR was justified as it was important to detect anyone who might have the virus in a new outbreak. False negatives were more of a concern. Now the same strategy is counter-productive when applied in vast numbers with poor quality control. There is a very real possibility that mass testing with contact tracing has led to an over estimate of true cases of infection. This then feeds back onto itself and drives more testing because every positive is treated as a case. Cases (positives) appear to drive much of government policy, modelling and public support for the measures. Since the restrictions have such enormous effects on society we must look very seriously at this possibility and act accordingly. New tests that do not suffer from such problems like the antigen tests could be considered to replace PCR for large scale use. Though this could be incredibly expensive. In any case, testing of asymptomatic individuals is not advisable unless the test is extremely accurate and the policy of mass testing itself should be questioned.

## Links

* [Dr. Clare Craig on COVID testing problems](https://www.youtube.com/watch?v=380DLg-nAqI)
* [Ch4 Dispatches report on a UK COVID-Testing Lab](https://www.channel4.com/press/news/dispatches-uncovers-serious-failings-one-uks-largest-covid-testing-labs)
* [Impact of false-positives and false-negatives in the UK's COVID-19 RT-PCR testing programme](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/895843/S0519_Impact_of_false_positives_and_negatives.pdf)
* [Dartmouth-Hitchcock: Epidemic That Wasn’t](https://www.nytimes.com/2007/01/22/health/22whoop.html)
* [Ireland's COVID-19 Data Hub](https://covid19ireland-geohive.hub.arcgis.com/)
* [HPSC report: Proposal for the management of weak positive PCR results](https://www.hpsc.ie/a-z/respiratory/coronavirus/novelcoronavirus/guidance/outbreakmanagementguidance/PCR%20weak%20results%20guidance.pdf)
