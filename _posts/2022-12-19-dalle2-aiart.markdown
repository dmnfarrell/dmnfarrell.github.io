---
layout: post
title:  "DALLE-2 and AI generated art."
date:   2022-12-19 18:00:00
categories: general
tags: [ai]
thumbnail: /img/dalle-abstract.png
---

## Background

DALL-E 2 (or just dalle) is an AI model from OpenAI used to generate art and imagery from text entered by humans. It is 'a generative language model' that takes sentences (called prompts) and creates corresponding original images. The actual method I cannot explain as I do not understand it well enough. Basically, it appears to be made of two models. One that converts the semantic meaning of some text into a vector space that is an image representation (CLIP). This is very roughly analagous to how a human would encode image representations of objects and can then generalise based on this. CLIP is trained on hundreds of millions of images and their associated captions. The other model called unCLIP 'decodes' the first model and converts this general representation into a specific image. This process is not deterministic so each time you give it a prompt, dalle will give you something different. A much more detailed explanation is [here](https://www.assemblyai.com/blog/how-dall-e-2-actually-works/).

## Example

You can try this yourself by registering on the OpenAI website which at the moment gives you 50 free credits. 'Prompts' are currently restricted to 400 words. Coming up with the appropriate prompt to get the image you want back is tricky. Supposedely being more specific is better to get what you want but it will depend on how well dalle can recognise and put together specific items. Names of real people are not allowed for example. The training set has also been filtered to remove certain material.

Let's say we give it the prompt: **"oil painting, irish people"**. I got the following results. You will get something else if you provide the exact same phrase.

<div style="width: auto;">
 <a href="/img/dalle-irish-people1.png"> <img class="small-scaled" src="/img/dalle-irish-people1.png"></a>
</div>

So when I repeated the prompt I got the following images, all different but you can tell they are trying to represent the same concepts. Notice how they vary in number of people and how close faces are etc. This is because I gave only a general phrase and only constrained it by asking for an oil painting hence the styles are similar. We can also request the image be in a known artistic style or of a specific artist. In fact there is a whole pdf book guide on how to make prompts, called the [DALLE2 prompt book](https://dallery.gallery/the-dalle-2-prompt-book/).

<div style="width: auto;">
 <a href="/img/dalle-irish-people2.png"> <img class="small-scaled" src="/img/dalle-irish-people2.png"></a>
</div>

## Variations

Another feature is to make variations on a generation you have done. So in theory you can iterate each time and end up with something completey different.

<div style="width: auto;">
 <a href="/img/dalle-variations-example1.png"> <img class="small-scaled" src="/img/dalle-variations-example1.png"></a>
</div>

## Other examples

Specific prompts certainly help to get more tailored responses. Most of the ones below were what I had in mind when I wrote the prompt. The more abstract the request, the less predictable the result may be. Click to see the full image.

<div class ="image-gallery">
    <div class="box"><a href="/img/dalle-liffey-sunrise.png" title="liffey">
       <img src="/img/dalle-liffey-sunrise.png" alt=file  class="img-gallery" />
     </a>
     <p class="caption">"oil painting of liffey dublin style of turner in high detail with sunrise"</p>
     </div>
     <div class="box"><a href="/img/dalle-xenomorph-vangogh.png" title="alien van gogh">
        <img src="/img/dalle-xenomorph-vangogh.png" alt=file  class="img-gallery" />
      </a>
      <p class="caption">"alien xenomorph, van gogh style"</p>
      </div>
     <div class="box"><a href="/img/dalle-mango-chair.png" title="mango chair">
        <img src="/img/dalle-mango-chair.png" alt=file  class="img-gallery" />
      </a>
      <p class="caption">"An armchair in the shape of a mango"</p>
      </div>
     <div class="box"><a href="/img/dalle-martian-julesverne.png" title="jules verne">
        <img src="/img/dalle-martian-julesverne.png" alt=file  class="img-gallery" />
      </a>
      <p class="caption">"fantasy oil painting of martian invasion jules verne, style neoclassicism 19th century"</p>
      </div>
     <div class="box"><a href="/img/dalle-hobbit-house.png" title="hobbit house">
        <img src="/img/dalle-hobbit-house.png" alt=file  class="img-gallery" />
      </a>
      <p class="caption">"hobbit house with garden, oil painting"</p>
      </div>
     <div class="box"><a href="/img/dalle-chessboy.png" title="chess boy">
        <img src="/img/dalle-chessboy.png" alt=file  class="img-gallery" />
      </a>
      <p class="caption">"African boy in karate costume playing chess. 6 years old. Soft lit background. painting"</p>
      </div>      
 </div>

## Limitations

On the one hand dalle is incredibly impressive in the results it can produce so rapidly. It creates a completed picture with disparate elements combined together in a visually meaningful way that we can generally interpret. This means the model can adequately encode the separate elements in 'image space'. Remember that these are new images that have never been seen before even if they are very similar in some cases to those in the enormous training set. However it won't always put images together in a way that we intend. In particular imparting spatial relations or getting it to include specific numbers of items/people is hard. These things humans have an intuitive understanding of, given a description. The model doesn't have an understanding of the physical world. This is probably because it is not trained to encode such information in its model. Perhaps it is too obvious to state that Dalle is good at composition in paintings only because it has copied from human artists who are.

For example if we give dalle this prompt: **"oil painting of dog kicking soccer ball at huge goal, cloudy sky"**, it is quite clear what we mean. A child could draw this simple picture even if crudely. But dalle could not produce these elements in an image accurately and typically gives the image below, with a huge ball instead. Perhaps the prompt was not specific enough of course.

<div style="width: auto;">
 <a href="/img/dalle-dog-ball.png"> <img class="small-scaled" src="/img/dalle-dog-ball.png"></a>
</div>

Here is an example from David Madras about how hard spatial composition is for dalle:

<div style="width: auto;">
<blockquote class="twitter-tweet"><p lang="en" dir="ltr">The ways in which <a href="https://twitter.com/hashtag/dalle?src=hash&amp;ref_src=twsrc%5Etfw">#dalle</a> is so incredible (and it is) really put a fine point on the ways in which compositionality is so hard <a href="https://t.co/I6DC4g53MK">pic.twitter.com/I6DC4g53MK</a></p>&mdash; David Madras (@david_madras) <a href="https://twitter.com/david_madras/status/1512573390896480267?ref_src=twsrc%5Etfw">April 8, 2022</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>
</div>

## Ethical problems

The huge advancement represented by this tool is very real but there is a tendency to over hype this technology especially on first seeing it work. It is amazing but can this really 'replace' human artists? Obviously it can be used to make art of a kind. That art may be indirectly derived from the work of real artists on which the models are trained. But arguably it can't go beyond it's training set. I don't think artists will be replaced but they are right to feel under threat of losing their income. There is also an element of exploitation here. Many living artists rail against the use of their work this way as shown by protests on the ArtStation website. DeviantArt [doesn't include](https://www.deviantart.com/team/journal/Create-AI-Generated-Art-Fairly-with-DreamUp-933537821) users art by default in training its **DreamUp** AI image-generation tool. It's hard to see how the AI plagiarism problem can be solved ultimately though.

There is also the problem of fake and offensive imagery. Dalle can make photorealistic images of places and even people. Though sometimes they look downright weird. In the future faked generated images will be of better quality and we won't be able to tell them apart from real ones. Images of someone famous performing revolting or illegal acts might be easy to make. This will have to dealt with as the problem arises. Whether we like it or not AI art is here and will only get better.

## Links

* [DALL·E 2, Explained: The Promise and Limitations of a Revolutionary AI](https://towardsdatascience.com/dall-e-2-explained-the-promise-and-limitations-of-a-revolutionary-ai-3faf691be220)
* [Hierarchical Text-Conditional Image Generation with CLIP Latents](https://arxiv.org/abs/2204.06125)
* [Artists stage mass protest against AI-generated artwork on ArtStation](https://arstechnica.com/information-technology/2022/12/artstation-artists-stage-mass-protest-against-ai-generated-artwork/)
* [OpenAI labs](https://labs.openai.com/)
* [MidJourney](https://www.midjourney.com/home)
* [Will AI replace artists?](https://www.youtube.com/watch?v=Ao5aV4I4Ivw)
