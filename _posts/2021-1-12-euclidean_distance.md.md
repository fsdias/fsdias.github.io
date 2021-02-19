---
title: "Calculate Euclidean distances in R"
date: 2021-1-12
tags: [gis,r,Riverscale]
header:
  image: "images/net_distance/banner_net_dist.png"
excerpt: "GIS, River networks, Euclidean distance,Riverscale"
mathjax: "true"
---

In a [previous post](https://fsdias.github.io/network_distance/), I explained how to calculate network distances using point locations of vegetation censuses in a river network. Now, I'm going to do the same with euclidean distances. Euclidean distance is the straight line distance between two sites.

We need a shapefile with the point location of the vegetation censuses. Important notes about this shapefile:

* it has to have at least the following columns: id) a unique identifier, x) a column with the "x" coordinates, y) a column with "y" coordinates. 
* Use a metric coordinate reference system. Do not use WGS84 (EPSG 4326).

Click here to download the shapefile([download](https://github.com/fsdias/blog_data/tree/main/network_distance)). This map shows the location of the points.

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

and import the shapefile using the shapefile() function from the "raster" package
```r
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

