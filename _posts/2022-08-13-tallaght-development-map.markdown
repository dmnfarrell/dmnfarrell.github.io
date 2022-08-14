---
layout: post
title:  "Mapping the historical development of Tallaght"
date:   2022-08-13 07:00:00
categories: general
tags: [maps,geopandas,python]
thumbnail: /img/tallaght_square_91.jpg
---

## Background

<div style="width: 320px; float:right;">
<img src="/img/tallaght_square_91.jpg" width="280px">
</div>

Tallaght is large surburb of Dublin about 13 km southwest of Dublin city, near the foothills of the Wicklow Mountains. Originally founded as a monastic settlement in 769 AD, it later became an important defensive outpost along the 'Pale' boundary. It remained a rural village until the 1960s when the Irish government commissioned town planner Myles Wright to devise an expansion plan for Dublin City. The Myles Wright Plan which was broadly adopted resulted in the creation of the three new towns of Tallaght, Clondalkin/Lucan and Blanchardstown to the west of the City. About 1970 onwards saw very rapid housing development over a large area around the old village. There have been successive waves of expansion since. Early and even subsequent development was plagued by frankly awful planning decisions (or shoddy implementation) leading to severe lack of amenities, insufficient transport links and limited retail services. A rail link was proposed but never constructed. It was only in the late 1980s that a large town centre was built, including the Square shopping centre, large library and hospital. A tram connection (Luas) to the city was also built in the 2000s.

Rapidity of development in the 70s is shown in the aerial views below. The old village and Priory is in the top left corner with estate to the south being Old Bawn and Bancroft to the west. Bancroft was first to be built in the 1960s.

<div style="width: auto;">
 <a href="/img/tallaght_aerial.jpg"> <img class="small-scaled" src="/img/tallaght_aerial.jpg"></a>
   <p class="caption">Tallaght aerial views (taken by Irish Air Corps - via Tallaght historical society).</p>
</div>

## Mapping

To map the development I used QGIS to download all OSM elements in the geographical area. These were grouped by a label such as name of the housing estate. Each area could then by assigned a year of construction using geopandas. Random dates within each group were also assigned for animation purposes (below). See the [Jupyter notebook](https://github.com/dmnfarrell/teaching/blob/master/geo/tallaght.ipynb) for associated code if interested.

## What was built when

Construction dates are not in the OSM data. It's not always easy to find out when individual areas were built though admittedly I didn't try very hard. So some of the times are rough guesses from knowledge of the area. Tallaght estates and commercial developments are shown below and coloured by date. Note that Firhouse is not part of Tallaght. Citywest, Saggart and Rathcoole are included to show the ongoing westward development. Roads are omitted for clarity. Note the complete lack of any substantial green belts separating any areas within Tallaght. The green space between Rathcoole/Saggart and Citywest is farm or 'waste' land that represent the last green banks between these urban areas. They are currently earmarked for more houses or a private 18 hole golf course!

<div style="width: auto;">
 <a href="/img/tallaght_areas.png"> <img class="scaled" src="/img/tallaght_areas.png"></a>
   <p class="caption">Tallaght main estates and approximate construction years.</p>
</div>

The image below shows roads. The M50 is to the far right.

<div style="width: auto;">
 <a href="/img/tallaght_detailed.png"> <img class="scaled" src="/img/tallaght_detailed.png"></a>
   <p class="caption">Tallaght main housing estates and nearby areas.</p>
</div>

## Animation

We can then attempt to crudely animate the development since 1960 to the present day by plotting structures built up to a given date. Note again this isn't exact as some years are guesses. It also isn't fine grained in the sense that some buildings are lumped in to common areas and given a year. It would be much easier if OSM editors added a year built attribute! To edit the OSM data adding correct dates would be a lot of work.

<div style="width: auto;">
  <video controls width="100%">
    <source src="/img/tallaght_growth.mp4" type="video/mp4">
    </video>
</div>


## Links

* [Tallaght](https://en.wikipedia.org/wiki/Tallaght)
* [Tallaght historical society](https://www.facebook.com/tallaghthistoricalsociety/)
* [Tallaght: The Planning and Development of an Irish New Town](/other/JIUSv3n12004_4.pdf)
* [Jupyter notebook](https://github.com/dmnfarrell/teaching/blob/master/geo/tallaght.ipynb)
