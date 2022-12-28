---
layout: post
title:  "Can ChatGPT solve bioinformatic problems with Python?"
date:   2022-12-21 13:00:00
categories: bioinformatics
tags: [python]
thumbnail: /img/man-at-computer.png
---

## Background

<div style="width: 320px; float:right;">
 <img src="/img/man-at-computer.png" width="300px">
</div>

There has been a lot of talk about OpenAI's new chatbot, ChatGPT. This is basically a very advanced chatbot. How it works is beyond my ability to explain. It is trained on a huge amount of information from the internet and can answer general questions or write poems and essays. It can also code in virtually any language quite well. You will see plenty of youtube videos marvelling at it's ability to produce (sometimes) usable code upon description of a specific coding task. There is plenty of talk of this technology replacing human programmers. This might not be true yet but it surely does have many applications in coding education, bioinformatics, statistics and many other fields. Not least of which is for learning purposes. It could even become a replacement for using search engines.

A few points to note:

* ChatGPT will often not give the same answer to an identical question so it's sometimes worth asking it twice. It will also correct itself if you point out an error in the code it's given you.
* It doesn't seem to distinguish between different versions of APIs or software versions. So some solutions will be correct but out of date for the version you might have.

## Example 1

Let's start with something simple. We ask **"write python function to get the reverse complement of a dna sequence"**. The result returned below works fine:

<div style="width: auto;">
 <a href="/img/chatgpt-python.png"> <img class="small-scaled" src="/img/chatgpt-python.png"></a>
</div>

## Example 2

Here's another simple one: **"write a python function that counts the number of reads in a gzipped fastq file"**. This one works quite well though it might not be optimal. It's certainly something a human would write.

<div style="width: auto;">
 <a href="/img/chatgpt-python4.png"> <img class="small-scaled" src="/img/chatgpt-python4.png"></a>
</div>

## Example 3

A bit more challenging is this one: **"write a python function that predicts all the open read frames in a genome sequence"**. This request is self explanatory. Here is the answer:

<div style="width: auto;">
 <a href="/img/chatgpt-python3.png"> <img class="small-scaled" src="/img/chatgpt-python3.png"></a>
</div>

This actually returns a result without errors. I didn't specify the input so it assumes just a long string of DNA which is reasonable. However it returns 50,000 ORFs from a bacterial sequence in which there are about 4000! I also asked it to predict the protein coding sequences and it gave much the same answer. In this cases we have probably under-specified the problem and need to give it more details.

## Example 4

To finish with something more complex, here is what happens when we ask ChatGPT to **"write a python function that converts a genbank file into a pandas dataframe and stores the gene features in each column"**. Genbank files are just text files with records in non tabular format. There is a 'feature' for each gene with 'qualifiers' (fields) being information about the gene like it's name. This is a quite a short prompt but it makes a good try at providing the answer, shown below. I put the code as images here so that someone doesn't come to this site and copy/paste it by mistake.

<div style="width: auto;">
 <a href="/img/chatgpt-python1.png"> <img class="small-scaled" src="/img/chatgpt-python1.png"></a>
</div>

In fact this solution doesn't work. This might be because the prompt is ambiguous. It returns a single row with each column storing the fields for all genes. So if we ask this: **"write a python function that converts a genbank file into a pandas dataframe with a row for each feature and each qualifier in a column"**. It returns the following:

<div style="width: auto;">
 <a href="/img/chatgpt-python2.png"> <img class="small-scaled" src="/img/chatgpt-python2.png"></a>
</div>

This is a better try but also fails because it seems confused over how to construct the dataframe in the loop. Also `append` is deprecated in the Pandas library. This could be a general problem of training on sites like stackoverflow with outdated answers.

A programmer could use these functions to build a working one but it's arguable that is easier than doing it from scratch.

## What people are saying

We are probably prone to over hype the impact of this technology because it's so startling at first glance. The truth is that this cannot currently replace an experienced programmer. As shown in [this video](https://www.youtube.com/watch?v=McXjNZFkQFU) by Conner Ardman you cannot just copy and paste code. You might spend more time trying to figure out if the answer it gives you is correct than to do it yourself. This is especially true of the challenging problems that you'd really want it to solve for you. The comments on this video are also revealing and a sample of what people are saying:

```
"I think I will use this a lot, but for simple things where I don't remember the exact syntax or implementation and just need a reminder. Also, sometimes when I'm programming it can be something simple that I overlooked, and then I could just copy&paste into chatGPT and have it find it for me (which is what I used it for yesterday).."
```

```
"As much as I've played around with it, I've realized that you have to be super specific in the way you submit your request. Include everything you'd want it to do first-hand."
```

```
"After playing for hours I can safely say that the main problem with the code provided by ChatGPT is that it gives incorrect answers not all the times but "sometimes". This "sometimes" is what causes trust issues..."
```

```
"As much I love the tech behind chatGPT, it still requires somebody with enough brains to handle its output. An secondly, If you pass in new tech frameworks or something just new, it can't handle that because according to its reinforcement training, it can't really learn from new upcoming things."
```

## Bionformatics type things chatGPT can help you with

Software engineers, scientists who have no coding experience are not going to be using chatGPT to solve novel or complex problems. It will help beginners but **to solve real world problems you do have to know enough to ask the right questions**.

Here's the kind of things it will good for:

* Ask it clearly defined functions with fairly simple tasks and it might produce a finished product, but you need to check it.
* Providing a scaffold for functions very quickly which you can then bug fix, add details or augment to suit. This would take care of the work of naming variables and writing logic for example.
* Building mundane things like regular expressions which humans don't find it easy to construct.
* Helping to clean or import data tables.

ChatGPT will be a great complement to learning coding. The challenge is to recognize the tasks where it can save you time. In future it will be able to do even more and maybe coders will be redundant but that won't be for some time. Also remember that it's free for now but probably will be monetized in the future. There will be free and open source versions though.

## Links

* [Can AI Do My Software Engineering Job For Me?](https://www.youtube.com/watch?v=McXjNZFkQFU)
