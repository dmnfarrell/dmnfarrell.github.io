---
layout: post
title:  "Diarization with OpenAI whisper and pyannote.audio"
date:   2023-11-12 19:00:00
categories: general
tags: [ai,python]
thumbnail: /img/whisper.png
---

## Background

OpenAIs [whisper](https://github.com/openai/whisper) library is an effective and free means of doing speech to text analysis. It's easy to use once installed and will output a set of files with timestamps for each sentence spoken. This is ideal for subtitling videos. It does identify individual speakers however, or group the conversation into passages according to who is speaking. This process is called diarization and can be acchieved using the [pyannote-audio](https://github.com/pyannote/pyannote-audio) library. This is based on PyTorch and hosted on the huggingface site. Here is some code for using it, mostly adapted from code from [Dwarkesh Patel](https://www.youtube.com/channel/UCXl4i9dYBrFOabk0xGmbkRA). To do this you need a recent GPU probably with at least 6-8GB of VRAM to load the medium model. I used an nvidia RTX 3060 for this.

## Install

It's recommended you make a new virtual environment for installing the packages. This will isolate them from your other Python packages if needed.

```
pip install openai-whisper pyannote.audio
```

## Code

This code will take an audio file and convert it to mono using ffmeg, then use whisper to transcribe it. The voice segments are delineated using the `PretrainedSpeakerEmbedding` model. The clustering algorithm then fits the embeddings to assign each segments to a speaker accordingly. Finally the output can be written as a transcript. This is the process as I understand it at least.

```python
import glob, os, subprocess
import torch
import whisper
import pyannote.audio
from sklearn.cluster import AgglomerativeClustering
from pyannote.audio import Audio
from pyannote.core import Segment
import wave
import contextlib
import numpy as np
import pandas as pd

from pyannote.audio.pipelines.speaker_verification import PretrainedSpeakerEmbedding
embedding_model = PretrainedSpeakerEmbedding( 
    "speechbrain/spkrec-ecapa-voxceleb",
    device=torch.device("cuda"))

def extract_speakers(model, path, num_speakers=2):
    """Do diarization with speaker names"""
    
    mono = 'mono.wav'
    cmd = 'ffmpeg -i {} -y -ac 1 mono.wav'.format(path)
    subprocess.check_output(cmd, shell=True)
    result = model.transcribe(mono)
    segments = result["segments"]
    
    with contextlib.closing(wave.open(mono,'r')) as f:
      frames = f.getnframes()
      rate = f.getframerate()
      duration = frames / float(rate)
        
    audio = Audio()
    def segment_embedding(segment):
        start = segment["start"]
        # Whisper overshoots the end timestamp in the last segment
        end = min(duration, segment["end"])
        clip = Segment(start, end)
        waveform, sample_rate = audio.crop(mono, clip)
        return embedding_model(waveform[None])

    embeddings = np.zeros(shape=(len(segments), 192))
    for i, segment in enumerate(segments):
      embeddings[i] = segment_embedding(segment)
    embeddings = np.nan_to_num(embeddings)
    
    clustering = AgglomerativeClustering(num_speakers).fit(embeddings)
    labels = clustering.labels_
    for i in range(len(segments)):
      segments[i]["speaker"] = 'SPEAKER ' + str(labels[i] + 1)
    return segments    

def write_segments(segments, outfile):
    """write out segments to file"""
    
    def time(secs):
      return datetime.timedelta(seconds=round(secs))
    
    f = open(outfile, "w")    
    for (i, segment) in enumerate(segments):
      if i == 0 or segments[i - 1]["speaker"] != segment["speaker"]:
        f.write("\n" + segment["speaker"] + ' ' + str(time(segment["start"])) + '\n')
      f.write(segment["text"][1:] + ' ')
    f.close()
```

We can then just use a few lines of code to run the process. You can select another model (e.g. 'small', 'large'). This example uses the English model (.en). I found that if you don't specify this results are of lower quality (assuming that's the language being spoken). You can also use a 'large' model for better accuracy but it will take more GPU memory.

```python
model = whisper.load_model('medium.en')
seg = extract_speakers(model, 'voices.wav')
write_segments(seg, 'transcript.txt')
```

## Example

Note the results won't be perfect of course. Sometimes whisper will overshoot so you can garbage at the end. Here is an example of a short interview excerpt and the transcript it produces below without editing.

<figure>
  <figcaption>Listen to excerpt:</figcaption>
  <audio controls src="/other/vidal.mp3">
    <a href="/other/vidal.mp3"> Download audio </a>
  </audio>
</figure>

### Transcript:

SPEAKER 1 0:00:00
A strange and volatile state, a whole Pacific nation with 22 million people, and we have a governor who is a Zen Buddhist Jesuit. Now, this is practically a resistible combination. He's called Jerry Brown, and his second go at the governorship, his second term is over. He wants to go to the United States Senate, and I would like to go to the United States Senate, so we may be, we may be, we have to always say this very carefully, we are both thinking seriously about running in 1982, which is a little over a year away, but you plan early, and you run hard. In England, you stand for election, we run. That is one of the many differences between our great lands. \
SPEAKER 2 0:00:46
Every American boy is told that he could become president. Do you still believe you could become president? \
SPEAKER 1 0:00:53
Every American girl now is also told she could become president. The presidency, I don't think is anything you really want to have now, because the political situation is such that the actual governing of the United States is not done by the president or by the Congress. It's done by the Chase Manhattan Bank, a very important bank, the oil companies, the defense industries, people in boardrooms around the United States. We have a law, really, that the better known you are, the more publicity you get in a country like the United States, the less power you have. Countries run by people with names like Mellon, that nobody knows anything about. So what they do is they hire lawyers, sneaky lawyers like Richard Nixon, and they put him in to be the president. He gets all the attention, all the publicity, but they do the actual governing. \
SPEAKER 2 0:01:47
Why do you want to be in politics? \
SPEAKER 1 0:01:49
Well, because first of all, I'd be the candidate against the Chase Manhattan Bank and the candidate against the oil companies, which means an awful lot of money would be put up to see that I was defeated, and there might be enough put up to affect that. On the other hand, occasionally you have a chance for a citizen politician to come in. One person isn't going to do much of anything in the United States Senate, but you have a chance to sort of say who's stealing what from whom and try to expose the system from within. I think it's a useful thing to do. I'm so genetically arranged that I'm obliged to be a candidate. It isn't really that it's a thing that you would sanely choose. Anyway, I'm a natural politician, alas. \
SPEAKER 2 0:02:33
Then do you think you've wasted your life? \
SPEAKER 1 0:02:35
I haven't wasted my life, but I would like to have been two people. I certainly never wanted to be a writer, but I was born one. It's just like being double-jointed or having perfect pitch. It's something that you're born with. You don't necessarily want to do it. As I began to read my first book, called The Duck and the Kangaroo, I was about six years old. It was a tale of unnatural vice between a duck and a kangaroo. I found myself writing a book at the same time I was reading one. \
SPEAKER 2 0:03:07
So you were born not only with a silver spoon in your mouth, but a fountain pen in your right baby fist. \
SPEAKER 1 0:03:12
It wasn't a silver spoon, I think silver tongue. \
SPEAKER 2 0:03:14
That must have been very uncomfortable, wouldn't it? \
SPEAKER 1 0:03:17
I think silver tongue is the phrase you're looking for. My family had no money. We were political and connected with just about everybody by marriage, but my end was penniless. So all we had were our wits to live by. \
SPEAKER 2 0:03:31
You weren't a privileged child. I always thought you were. \
SPEAKER 1 0:03:33
No, I was privileged in the sense that I was brought up in a family of senators and military people and politicians. My family managed to marry just about everybody, but there were certain sections of the family that were very rich and others that were not. Ms. Onassis, as she is known now, is also in my same situation. We both had a stepfather who was extremely rich, but neither of us had a penny, so she had to pursue wealth in her way, and I had to write for a living and do other things to make a living. \
SPEAKER 2 0:04:10
Is wealth the only class or caste that exists in America? \
SPEAKER 1 0:04:15
I would say no. We have a complex class system, like all democratic societies that have no classes. To Mars it's so complicated that you almost can't describe it. I think the nice thing about the political arrangement in America is that you can go all through your life without ever knowing anything about what class you belong to, or if you get interested in it, it's like Proust. It's all Byzantium. \
SPEAKER 2 0:04:39
But it's based upon what, then, power? \
SPEAKER 1 0:04:41
It's based mostly upon money, as it is practically anywhere else, even though the great English families, as they lose their money, then they have their glorious names and their schools, which maintains your class system. We have our schools, we have our manager class. I just went back to my old school, Exeter, which is what you would call a public school here. It must be not a grand one like Eton. It's something like Winchester. It's sort of for rather bright boys on the make who will then become presidents of banks and things. I realized I hadn't been back since I graduated, and they were having their 200th birthday, the school. As I looked around, I said, you know, it hadn't occurred to me what a deliberate thing the school was, how deliberately we were being shaped to be the managers of the United States, of the banks, and of the brokerage houses, and of the law firms. Not so much in politics. That was a little vulgar. We were to be the people behind the scenes running the country. And I thought, really, they knew what they were doing, the American ruling class, since we were not the rich boys, but we were the ones who were considered rather sharp. 

## Links

* [notebook](https://github.com/dmnfarrell/teaching/blob/master/machine_learning/diarization.ipynb)
* [colab example](https://colab.research.google.com/drive/1V-Bt5Hm2kjaDb4P1RyMSswsDKyrzc2-3?usp=sharing#scrollTo=buGt4moR5Mac)
* [whisper](https://github.com/openai/whisper)
* [pyannote](https://github.com/pyannote/pyannote-audio)
* [Example on youtube](https://www.youtube.com/watch?v=MVW746z8y_I&ab_channel=1littlecoder)