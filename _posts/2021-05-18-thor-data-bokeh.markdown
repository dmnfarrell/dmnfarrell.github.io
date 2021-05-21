---
layout: post
title:  "Viewing the THOR dataset with Bokeh and Panel"
date:   2021-05-18 12:16:00
categories: plotting
tags: [plotting,bokeh,panel,THOR]
thumbnail: /img/F-105Ds_bomb_vietnam.jpg
---

## Background

Following from the [previous post](/plotting/thor-seasia), we can view the same THOR data interactively with Bokeh and Panel. These are Python libraries for interactive web visualization of data. The data is read in as in the previous post and we then write a function to plot the map data points using Bokeh with a Tile Renderer. We also use a function here to convert lat and long coordinates into web mercartor. The dashboard then consists of a date selector and some filter widgets that are connected to the `update_map` function via the `param.watch` mechanism. The code is all below in one block.

## Result

The resulting widget is animated here.

<div>
 <a href="/img/thor_bokeh_dashboard.gif"> <img class="small-scaled" src="/img/thor_bokeh_dashboard.gif"></a>  
   <p class="caption">.</p>
</div>

The complete code can be found [here](https://github.com/dmnfarrell/teaching/tree/master/vietnam_thor). You can also view a demo of the dashboard hosted [here](http://bola.ucd.ie/thor). If this no longer available you can run the notebook from [mybinder](https://mybinder.org/v2/gh/dmnfarrell/teaching/tree/master/vietnam_thor/HEAD) too.

## Code

 ```python
 colormap={'NORTH VIETNAM':'brown','SOUTH VIETNAM':'orange','LAOS':'red',
                 'CAMBODIA':'green','THAILAND':'blue','UNKNOWN':'gray'}

 providers = ['CARTODBPOSITRON','STAMEN_TERRAIN','OSM','ESRI_IMAGERY']
 cats = ['TGTCOUNTRY','WEAPONTYPE','MFUNC_DESC']

def wgs84_to_web_mercator(df, lon="LON", lat="LAT"):
    """convert mat long to web mercartor"""

    k = 6378137
    df.loc[:,"x"] = df[lon] * (k * np.pi/180.0)
    df.loc[:,"y"] = np.log(np.tan((90 + df[lat]) * np.pi/360.0)) * k
    return df

 def draw_map(df=None, long=None, lat=None, height=500, colorby='TGTCOUNTRY',
             point_size=5,
              tile_provider='CARTODBPOSITRON'):
    tile_provider = get_provider(tile_provider)
    tools = "pan,wheel_zoom,box_zoom,hover,tap,lasso_select,reset,save"
    sizing_mode='stretch_both'

    # range bounds supplied in web mercator coordinates
    k = 6378137
    pad = 700000
    if lat == None:
        lat = 16
    if long == None:
        long = 108
    x = long * (k * np.pi/180.0)
    y = np.log(np.tan((90 + lat) * np.pi/360.0)) * k

    p = figure(x_range=(x-pad, x+pad), y_range=(y-pad, y+pad),
               x_axis_type="mercator", y_axis_type="mercator", tools=tools,
               plot_width=height, plot_height=height, sizing_mode=sizing_mode)
    p.add_tile(tile_provider)
    if df is None:
        return
    df.loc[:,'color'] = [colormap[i] if i in colormap else 'gray' for i in df[colorby]]

    source = ColumnDataSource(df)    
    p.circle(x='x', y='y', size=point_size, alpha=0.7, color='color', source=source)
    p.toolbar.logo = None    
    p.title.text = "date"
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("TGTCOUNTRY", "@TGTCOUNTRY"),
        ("MSNDATE", "@MSNDATE"),
        ("TAKEOFFLOCATION", "@TAKEOFFLOCATION"),
        ("WEAPONTYPE", "@WEAPONTYPE"),
        ("MFUNC_DESC", "@MFUNC_DESC")     
    ])
    hover.formatters={'@MSNDATE': 'datetime'}
    return p

 def dashboard():
     cols = list(x.columns)
     colorby='TGTCOUNTRY'
     map_pane=pn.pane.Bokeh(width=700)
     df_pane = pn.pane.DataFrame(width=600,height=600)
     date_picker = pnw.DatePicker(name='Pick Date',width=200)
     from datetime import date  
     date_picker.value=date(1965, 1, 1)    
     date_slider = pnw.DateSlider(name='Date', start=dt.datetime(1965, 1, 1),
                                  end=dt.datetime(1973, 10, 31), value=dt.datetime(1968, 1, 1))      
     tile_select = pnw.Select(name='tile layer',options=providers,width=200)
     filterby_select = pnw.Select(name='filter by',value='',options=['']+cols[1:4],width=200)
     value_select = pnw.Select(name='value',value='',options=[],width=200)
     find_btn = pnw.Button(name='find in region',button_type='primary',width=200)

     def update_tile(event=None):
         p = map_pane.object
         p.renderers = [x for x in p.renderers if not str(x).startswith('TileRenderer')]
         rend = renderers.TileRenderer(tile_source= get_provider(tile_select.value))
         p.renderers.insert(0, rend)

     def update_filter(event):
         col=filterby_select.value
         if col=='':
             value_select.options = []
         else:
             value_select.options = sorted(list(x[col].dropna().unique()))

     def find_in_region(event):
         #get points in selected map area
         p = map_pane.object
         source = p.renderers[1].data_source
         d = x[(x.x>p.x_range.start) & (x.x<p.x_range.end) & (x.y>p.y_range.start) & (x.y<p.y_range.end)]
         #add any filter
         d = do_filter(d)
         if len(d)==0:
             return
         elif len(d)>25000:           
             p.title.text = 'too many points!'
         else:
             d.loc[:,'color'] = [colormap[i] if i in colormap else 'gray' for i in d[colorby]]        
             source.data = dict(d)
             p.title.text = 'selected %s points' %len(d)
         map_pane.param.trigger('object')     
         return

     def do_filter(d):
         col = filterby_select.value
         val = value_select.value
         if col != '':
             d = d[d[col]==val]
         return d

     def update_date(event):
         date_slider.value = date_picker.value     

     def update_map(event=None, date=None):
         p = map_pane.object
         source = p.renderers[1].data_source
         if date == None:
             date = str(date_slider.value)
         d = x[x.MSNDATE==date]
         d = do_filter(d)
         if len(d)==0:
             return  
         d.loc[:,'color'] = [colormap[i] if i in colormap else 'gray' for i in d[colorby]]
         source.data = dict(d)
         p.title.text = date

     sdate='1968-01-01'
     d = x[x.MSNDATE==sdate]
     map_pane.object=draw_map(d)

     date_slider.param.watch(update_map,'value')
     date_picker.param.watch(update_date,'value')
     tile_select.param.watch(update_tile,'value')
     filterby_select.param.watch(update_filter,'value')
     value_select.param.watch(update_map,'value')
     find_btn.on_click(find_in_region)

     dashboard = pn.Column(date_slider,pn.Row(
                pn.Column(date_picker,tile_select,filterby_select,
                value_select,find_btn),map_pane))
     return dashboard

 app=dashboard()
 ```

## Links

* [Panel](https://panel.holoviz.org/)
* [Bokeh](https://bokeh.pydata.org/)
* [THOR dataset](https://data.world/datamil/vietnam-war-thor-data)
