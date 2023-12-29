---
layout: post
title:  "A Panel app for image-to-image generation"
date:   2023-12-29 11:00:00
categories: general
tags: [ai,python]
thumbnail: /img/img2img_app.png
---

## Background

[Previously](/general/stable-diff-img2img) we saw how to implement the Stable Diffusion image-to-image model using the Python diffusers library. There are plenty of websites now that offer AI image generation since it has become so popular. This post simply shows how we can make your own basic web dashboard with [Panel](https://panel.holoviz.org/) that does something similar. The app allows someone to upload an image and generate new ones with a prompt and some of the settings previously demonstrated. Images are placed in a tabbed container and can be removed when needed.

## Code

First we create our dashboard app by defining a function that returns the widgets. The `execute` function calls the img2imgprompt method. In Panel we can lay out our widgets using `pn.Column` and `pn.Row` panes. Note that other required methods are placed in another tools.py file and imported. This is kept in the same directory.

```python
import os, glob
import random, math
import numpy as np
import pandas as pd
from PIL import Image
from tools import *

import panel as pn
import panel.widgets as pnw
pn.extension('tabulator', css_files=[pn.io.resources.CSS_URLS['font-awesome']])

def dashboard():
    
    def run(**kwargs):
        """Run prompt""" 
      
        filename = img2imgprompt(path='temp', n=1, **kwargs)
        name = kwargs['seed']
        add_image(filename)
        #print (filename)

    w=210    
    styles = ['','oil','impressionist','pencil','ink','watercolor','crayon drawing','digital art','pop art','cubism',
              'sculpture','craft clay','linocut','engraving','anime','studio photography','analog film',
              'abstract','pixel art','paper collage','isometric','lowpoly','origami']
    #title = html_pane = pn.pane.HTML("""<h2>image-to-image app</h2>""")
    file_input = pnw.FileInput(width=w,accept='.png,.jpg')    
    go_btn = pnw.Button(name='run',width=w,button_type='success')
    stop_btn = pnw.Button(name='stop',width=w,button_type='danger')
    prompt_input = pnw.TextAreaInput(name='prompt',value='',width=w)
    style_input = pnw.Select(name='style',options=styles,width=w)
    strength_input = pnw.FloatSlider(name='strength',value=.8,step=.01,start=.01,end=.99,width=w)
    guidance_input = pnw.IntSlider(name='guidance',value=5,step=1,start=0,end=10,width=w)
    seed_input = pnw.IntInput(name='seed',value=0,step=1,start=0,end=10000,width=w)
    progress = pn.indicators.Progress(name='Progress', value=0, width=w, bar_color='primary')    
    widgets = pn.Column(pn.WidgetBox(file_input,prompt_input,style_input,strength_input,
                                     guidance_input,seed_input,go_btn,progress), height=700, width=230)
                
    img_pane = pn.pane.Image(caption='start image', sizing_mode="stretch_width")
    tabs = pn.Tabs(closable=True, tabs_location='left',sizing_mode="stretch_width")

    def update_image(event):
        imgfile = file_input.value
        img_pane.object = imgfile
        
    def add_image(imgfile):
        #add new image 
        name = os.path.basename(imgfile)
        new = pn.pane.Image(imgfile, caption=name, sizing_mode="stretch_both")        
        tabs.append(new)
        tabs.active = len(tabs)-1
              
    def execute(event):
        #run the model with widget         
        img = file_input.filename
        seed = seed_input.value
        if seed == 0:
            seed = None
        if img == None:
            return
        progress.value=-1
        run(prompt=prompt_input.value, style=style_input.value, 
            init_images=[img],strength=strength_input.value, seed=seed)
        progress.value=0
        
    file_input.param.watch(update_image, 'value')
    go_btn.param.watch(execute, 'clicks')
    
    app = pn.Column(
                 pn.Row(widgets,tabs,pn.Column(img_pane,width=240),
                 sizing_mode='stretch_both',
                 styles={'background': 'WhiteSmoke'}))

    return app
```

Finally we create the app object and place it in a bootstrap template and call the `servable` method to launch it. This code is all placed in a python file. The full script is [here](https://github.com/dmnfarrell/teaching/blob/master/machine_learning/img2imgapp.py). If you wanted to run this in the background you could launch it using the `nohup` command. BY default Panel apps launched like this will run on localhost:5000.

```python
bootstrap = pn.template.BootstrapTemplate(title='image-to-image app',
                logo='flowers.jpg',header_color='blue')
pn.config.sizing_mode = 'stretch_width'
app = dashboard()
bootstrap.main.append(app)
bootstrap.servable()
```

Here is an example of the running application: 

<div style="width: auto;">
 <a href="/img/img2img_app.gif"> <img class="scaled" src="/img/img2img_app.gif"></a>  
  <p class="caption"></p>
</div>

For options on deplying such apps via the internet, see [this link](https://docs.bokeh.org/en/latest/docs/user_guide/server/deploy.html#ug-server-deploy).

## Links

* [Panel](https://panel.holoviz.org/gallery/index.html)
* [Bootstrap Panel examples](https://bootsnipp.com/tags/panel)
* [img2img](https://huggingface.co/docs/diffusers/using-diffusers/img2img)