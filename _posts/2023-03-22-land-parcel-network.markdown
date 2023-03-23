---
layout: post
title:  "Simulate land parcels and fragmentation with geopandas"
date:   2023-03-22 14:00:00
categories: python
tags: [maps,geopandas,python]
thumbnail: /img/farm_hedgerows.png
---

## Background

Land fragmentation is the division of land holdings into discrete parcels that may be dispersed over a wide area. It is common in an agricultural context where a single farmer can own multiple parcels of contiguous (separated by roads, hedgerows or fencing) or non-contiguous land. For cattle farming this means animals could spend their time in different portions of land at different times of the year. This presents a challenge for something like contact networks where determining connections between farms relies on knowing which parcels are close to others.
Here geopandas used to make groups of polygons that approximate such land parcels. We then use networkx to generate crude contact networks taking into account farms with neighbouring fragments.

## Code

The method used to generate random fragments here is to make a set of random points, derive the voronoi cells from these and then cluster them into 'farms'. We can merge contiguous polygons and optionally sparsify them to get less dense parcels. The more cells/points used the more sides the polygons will have. This is code verbose not very satisfactory. This is a classic case of quickly conceiving a way to solve a problem and then realising it's hard to add new functionality, namely adding more fragments, on top of the solution. So there are undoubtedly cleaner ways to accomplish this.

Extra methods not shown here are in the [notebook on github](https://github.com/dmnfarrell/teaching/blob/master/geo/land_parcels.ipynb) with the complete example.

```python
def generate_land_parcels(cells=100,herds=10,empty=0,fragments=0,seed=None):
    """
    Simulate land parcels with fragmentation.
    Args:
        cells: number of points to make initial cell polyons
        herds: number of farms
        empty: fraction of fragments that are empty
        fragments: number of empty cells to add back as fragments
    """

    from shapely.geometry import Point,MultiPoint,MultiPolygon,Polygon
    from libpysal import weights, examples
    from libpysal.cg import voronoi_frames
    from sklearn.cluster import KMeans

    n = cells
    k = herds
    if seed != None:
        np.random.seed(seed)
    x,y = np.random.randint(1,1000,n),np.random.randint(1,1000,n)
    coords = np.column_stack((x, y))
    cells, generators = voronoi_frames(coords, clip="extent")
    centroids = cells.geometry.centroid
    #cluster parcels into 'herds'
    kmeans = KMeans(n_clusters=k,n_init='auto').fit(coords)
    cells['cluster'] = kmeans.labels_

    #remove some proportion of cells randomly
    e = cells.sample(frac=empty, random_state=seed)
    cells.loc[e.index,'cluster'] = 'empty'

    #create new GeoDataFrame
    poly=[]
    data = {'cluster':[]}
    for c,g in cells.groupby('cluster'):
        if c == 'empty':
            continue
        poly.append(MultiPolygon(list(g.geometry)))
        data['cluster'].append(c)
    farms = gpd.GeoDataFrame(data=data,geometry=poly)

    #merge contiguous fragments of same herd
    farms = farms.dissolve(by='cluster').reset_index()
    #remove polygons with 'holes'
    def no_holes(x):
        if type(x) is MultiPolygon:
            return MultiPolygon(Polygon(p.exterior) for p in x.geoms)
        else:
            return Polygon(x.exterior)
    farms.geometry = farms.geometry.apply(no_holes)

    #remove empty cells inside parcels
    e = e[e.within(farms.unary_union)==False]

    #assign some of the empty cells as fragments
    for i,r in cells.loc[e.index].sample(fragments).iterrows():    
        poly = farms.iloc[0].geometry
        if type(poly) is MultiPolygon:
            geom = poly.geoms
            new = MultiPolygon(list(geom) + [r.geometry])
        else:
            geom = poly
            new = MultiPolygon([geom,r.geometry])
        farms.loc[0,'geometry'] = new

    #merge contiguous fragments again in case we added fragments
    farms = farms.dissolve(by='cluster').reset_index()

    farms['herd'] = farms.apply(lambda x: get_short_uid(4),1)
    farms['fragments'] = farms.geometry.apply(count_fragments)
    farms['color'] = random_hex_colors(len(farms),seed=seed)
    return farms
```

## Examples

The parameters used in combination might provide an appropriate result depending on what is required.

### Change number of farms/cells

To illustrate parameter usage we show the effect of increasing the number of initial cells. This determines the number of edges parcels will have. You can see this in the plots below where the number of farms is also varied.

<div style="width: auto;">
 <a href="/img/land_parcels_vary_cells.png"> <img class="scaled" src="/img/land_parcels_vary_cells.png"></a>  
  <p class="caption"></p>
</div>

### Change empty parameter

The `empty` parameter sparsifies the polygons as below. This also has the effect of fragmenting farms by removing a fraction of the original cells.

<div style="width: auto;">
 <a href="/img/land_parcels_vary_empty.png"> <img class="scaled" src="/img/land_parcels_vary_empty.png"></a>  
  <p class="caption"></p>
</div>

### Add more fragmentation

We can further fragment the farms by adding the `fragments` parameter. You need to use this in combination with the empty parameter since the fragments are added by allocating the empty cells to other random farms. If there aren't enough empty cells this will break as it currently is written. The number of fragments made won't always be exact either.

<div style="width: auto;">
 <a href="/img/land_parcels_vary_fragments.png"> <img class="scaled" src="/img/land_parcels_vary_fragments.png"></a>  
</div>

## Contact network from parcels

Finally we can generate a contact network from the parcels which is why I was doing this originally. This produces a network graph with edges between all nodes (centroids of parcels) with contiguous fragments. Therefore we would expect that very fragmented farms are highly connected in the network. That can be seen below. Not the code currently just adds all the fragments to the first farm but we can easily change to add randomly or otherwise.

<div style="width: auto;">
 <a href="/img/land_parcels_contact_network.png"> <img class="scaled" src="/img/land_parcels_contact_network.png"></a>
  <p class="caption">Contact network from connected parcels. The red farm on top is fragmented and most connected in the network</p>
</div>

Same with more farms below. This isn't very realistic but illustrative. Node locations are centroid of original parcels which can become meaningless with a lot of fragmentation.

<div style="width: auto;">
 <a href="/img/land_parcels_contact_network2.png"> <img class="scaled" src="/img/land_parcels_contact_network2.png"></a>
  <p class="caption">The red farm is by far most connected in the network.</p>
</div>

## Links

* [code on github](https://github.com/dmnfarrell/teaching/blob/master/geo/land_parcels.ipynb)
