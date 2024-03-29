---
layout: post
title:  "Educational software for data analysis"
date:   2015-06-14 12:00:00
categories: education
thumbnail: /img/education_computer.jpg

---

In developing a tool for simple data analysis that could be used in teaching it is interesting to look at the current options. The near universal use of spreadsheets, largely in the form of Excel, has meant they are being used for tasks not originally designed for. There is plenty of evidence of this even if one is not an educator. Just walk around a research lab in a university or a workplace office. If you were rude enough to stop and look over peoples shoulders you might witness Excel being used in turn as a statistical tool, database or batch programming tool. See this article on  [spreadsheet addiction](http://www.burns-stat.com/documents/tutorials/spreadsheet-addiction/) by Patrick Burns if the topic intrigues you. The author does not hate spreadsheets. In fact they are a wonderful idea but stretched to beyond breaking point in many cases. The macro in particular is a tool being used in place of actual programming in sometimes absurd ways.

Due to the changing nature and size of datasets businesses are transitioning away from spreadsheets to using proper business analytics tools, where required. Scientists (those not specialised in informatics) too are very slowly beginning to adapt to new tools like R. However it's easier said than done because people want something that can replace Excel and do the other stuff too. We don't really need to replace Excel, just start using a few other applications alongside it. Even most professional data scientists likely still use spreadsheets for pre-processing data or presenting reports.

So there is a gap to be filled - but with what? To get an idea of the rationale educators might use I would strongly recommend reading [this blog discussion](https://learnandteachstatistics.wordpress.com/2013/02/11/excel-spss-minitab-or-r/) on stats teaching software. The comments section in particular are informative. Though it concentrates on statistics many of the points apply to data analysis in general. The dichotomy between programming and visual tools is a concern for teachers. Why waste time confusing students with new programming concepts if you are teaching statistics (the reverse can also be said). This is one reason R might not gain traction with some students since it requires at least some programming skill, though R-studio is clearly an excellent free tool. Actually there is not as big a gap between Excel and fully fledged programming skills as one might think. This blog on carrying out [common excel tasks in pandas](http://pbpython.com/excel-pandas-comp.html) shows the same analyses applied in parallel by using Excel and Python/Pandas. It indicates a useful way to each programming to students already very familiar with spreadsheets. Admittedly understanding that lesson requires fairly advanced spreadsheet knowledge and understanding of data analysis already.

The above linked discussion also contains the oft repeated reason for using a popular (commercial) tool - that it is 'common in business' or 'industry standard'. This is a very poor reason to choose a piece of software for teaching in schools or universities. In fact it is likely half the reason spreadsheets are still so misused since it perpetuates the same mistakes. However many educators are forced (or prefer) to turn to commercial solutions. In data analytics there are not that many free programs available for non-programmers (see links below). There are now many web based solutions, most being commercial. The web app is certainly the new frontier for such tools but it still isn't quite that easy to build highly interactive applications. So there is a space for desktop tools too.

## Free software for data analysis

* [R-studio](http://www.rstudio.com/) Requires learning R language basics
* [iNZight](https://www.stat.auckland.ac.nz/~wild/iNZight/index.php) Built for schools
* [Tableau Public](https://public.tableau.com/s/) Free version of the business tool
