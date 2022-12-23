---
layout: post
title:  "How to host your podcast with github"
date:   2022-10-15 11:00:00
categories: general
tags: [github,jekyll]
thumbnail: /img/podcast_feed.jpg
---

## Background

<div style="width: 320px; float:right;">
 <img src="/img/podcast_feed.jpg" width="300px">
</div>

If you want to make your own podcasts the usual method is probably to use a [commercial service](https://www.podcastinsights.com/best-podcast-hosting/) that handles all the file hosting and distribution. However the process is actually surprisingly simple for those with a bit of technical knowledge and patience.

Web standards are really the key to how podcasting works. RSS (Really Simple Syndication) and its offshoot
[Atom](https://en.wikipedia.org/wiki/Atom_(web_standard)) are how sites create a 'feed' that indicates when updates are made. RSS is just XML-formatted plain text. This format can be read by a program that presents site content and updates in a convenient way. In the same way, podcast managers retrieve new episodes. Details like episode title, length, file location and description are stored in the feed xml. So a podcast manager is very like an RSS reader but deals with audio. So all we really need is to make the feed and link to where you have store the associated audio files. We can do all this inside a simple github hosted site as long as your audio files are not larger than 100MB.

## How it works

Here we use Jekyll for creating the website and generating the feed.xml file. It's commonly used for blogs. Jekyll is a static site generator that will generate the feed and other episode updates dynamically when you add new posts. You can use any web site generator though. In this example we use the jekyll-skeleton template provided by [Tim Klapdor](https://github.com/timklapdor). This was forked and some changes made for this example. The final website is [here](https://dmnfarrell.github.io/podcast-example). The feed url can be seen linked in that page under RSS. This is what other software will use to retrieve your podcast.

## Steps

<div style="width: 200px;float:right;">
 <a href="/img/podcast-scr1.png"> <img class="scaled" src="/img/podcast-scr1.png"></a>  
</div>

* Create a new repository with your chosen Jekyll template or clone/fork the example site linked here.
* Edit the _config.yml file to enter your own site details.
* Make a new post for your episode (see below) and add an audio file to the audio folder for it.
* Add an image for the episode, place in img folder or where suits best.
* Run `jekyll serve` to test your site.
* Commit your changes
* Activate github pages hosting to use the root of the repo. The site will then be generated. Go to the site address to test.
* You can check the feed is working by simply adding the feed url manually in any podcast manager. Try download/stream an episode to check it works.
* Register your podcast on whichever service you wish e.g. itunes, google podcast, spotify etc.

## Adding a new post

To add a new post just create a new file in the `_posts` folder. The post should contain at least the details below. You can use the example as a guide. Length is the size of the audio in bytes. The post doesn't have to contain any other text but you would probably have episode details in here. Your audio can be committed and pushed to github like any other file and linked in the file keyword. Notice the audio link is to the raw file, not the page on github. Add an image for the episode, place in img folder. This is the thumbnail keyword.

```
layout: post
title:  "What is Atom?"
date:   2022-09-19 11:23:34 +1000
file: https://github.com/dmnfarrell/podcast-example/raw/master/audio/track1.mp3
categories: atom update
description: "All about atom and rss - test."
duration: "00:50"
length: "1029183"
keywords: "atom"
thumbnail: /img/logo.png
```

## Links

* [Use GitHub Pages & Jekyll to host a podcast](https://wiobyrne.com/use-github-pages-jekyll-to-host-a-podcast/)
* [Apple podcast feed requirements](https://podcasters.apple.com/support/823-podcast-requirements)
* [Github repo for this example](https://github.com/dmnfarrell/podcast-example)
