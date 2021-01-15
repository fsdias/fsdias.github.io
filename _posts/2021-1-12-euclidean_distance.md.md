---
title: "Calculate euclidean distances in R"
date: 2021-1-12
tags: [gis,r,Riverscale]
header:
  image: "images/net_distance/banner_net_dist.png"
excerpt: "GIS, River networks, Euclidean distance,Riverscale"
mathjax: "true"
---

In a previous post, I explained how to calculate network distances using point locations of vegetation censuses in a river network. Now, I'm going to do the same with euclidean distances. Euclidean distance is the straight line distance between two sites.

As in the previous post, we will need two shapefiles to start, one with the river network and another with the point location of the vegetation censuses. Important notes about the shapefiles:

* the river shapefile can't have any topological problems. I recommend that you use QGIS and run the "Fix geometries" algorithm to fix any issues. It's available in the Processing toolbox.
* the points shapefile needs to have at least the following columns: id) a unique identifier, x) a column with the "x" coordinates, y) a column with "y" coordinates. 
* Use a metric coordinate reference system. Do not use WGS84 (EPSG 4326).

We will be working with "rivers.shp" and "points.shp"([download](https://github.com/fsdias/blog_data/tree/main/network_distance)). Map below.

<img src="{{ site.url }}{{ site.baseurl }}/images/net_distance/river_map.png" alt="linearly separable data">


We start by loading the following packages:
```r
library(raster)
library(igraph)
library(shp2graph)
library(dplyr)
library(reshape2)
library(sp)
```

and import both shapefiles using the shapefile() function from the "raster" package
```r
river<-shapefile("river.shp")
points<-shapefile("points.shp")
```

To calculate euclidean distances we use the "spDists" from the "sp" package.

```r
dst<-spDists(points,points)
```

Then we edit the resulting matrix, remove the diagonal and convert it to a data.frame.

```r
dst<-as.matrix(dst)
rownames(dst) <- points$id
colnames(dst) <- points$id
dst[upper.tri(dst)]=NA
euc_distances<-reshape2::melt(dst, na.rm=T)
colnames(euc_distances)<-c("item1","item2","euc_dist")
euc_distances<-filter(euc_distances,item1!=item2)
```

<img src="{{ site.url }}{{ site.baseurl }}/images/euc_dist/euc_dist_table.png" alt="linearly separable data">


**Filipe Dias**

