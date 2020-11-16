---
title: "Calculate network distances in R"
date: 2018-01-28
tags: [gis,r]
header:
  image: "/images/perceptron/percept.jpg"
excerpt: "GIS, River networks, Network distance"
mathjax: "true"
---

A while ago, I had to calculate network distances between vegetation samples collected in different river network locations. There are some tools in QGIS and GRASS GIS that can perform these calculations. Still, I decided to use R because it allows me to format the output structure.

We need two shapefiles to start, one with the river network and another with the vegetation samples' point location. Important notes about the shapefiles:

* the river shapefile can't have any topological problems. I recommend that you use QGIS and run the "Fix geometries" algorithm to fix any issues. It's available in the Processing toolbox.
* the points shapefile needs to have at least the following columns, one with a unique identifier (e.g., ID), a column with the "x" coordinates, and another one with the "y" coordinates. Use a metric coordinate reference system, do not use WGS84 (EPSG 4326).

In this post we will be working with "rivers.shp" and "points.shp"([download](https://github.com/fsdias/blog_data/tree/main/network_distance)). Map below.

<img src="{{ site.url }}{{ site.baseurl }}/images/net_distance/river_map.png" alt="linearly separable data">


We start by loading the following packages:
```r
library(raster)
library(igraph)
library(shp2graph)
library(dplyr)
library(reshape2)
```

Next, we import both shapefiles using the shapefile() function from the "raster" package
```r
river<-shapefile("river.shp")
points<-shapefile("points.shp")
```

and we check the connectivity of the river network.
```r
river = nt.connect(river)
```

Now we add the vegetation sample points to the river network. We use "approach=2" to map each point to the nearest point on the network and add them as new nodes/vertices.

```r
  rivernel = points2network(ntdata=river,
                            pointsxy=coordinates(points),
                            ELComputed=TRUE,
                            approach=2, Detailed=TRUE, ea.prop=c(1,0))
```
The next step is to convert the "rivernel" object into an "igraph" object".
```r
    rivg = nel2igraph(rivernel[[1]], rivernel[[2]], weight=rivernel[[8]])
```

Each vegetation sample point was added to the river network and received a unique identifier. We can access and save those unique identifiers by running:

```r
      ipoints = rivernel[[3]]
```

Now, we create a table that with the "x,y"" coordinates of each of our vegetation samples alongside their original "id"" (the one in the shapefile) and the "id" that points2network() assigned them. I called this new "id" "cod_p2n".

```r
      a<-data.frame(cod_p2n=seq_along(1:length(rivernel[[4]])),x=floor(rivernel[[4]]),y=floor(rivernel[[5]]))
      a<-filter(a, cod_p2n >= min(ipoints)) #this works because new points are added to the end
      df<-data.frame(x=floor(points@data$x),y=floor(points@data$y),id=points@data$id)    
      table<-merge(a,df,by=c("x","y")) 
      table<-arrange(table,cod_p2n) 
```

In the last step we create a distances matrix and convert it into a data.frame.

```r
  mat_dist<-distances(rivg,ipoints,ipoints)
  mat_dist<-as.matrix(mat_dist)
  rownames(mat_dist) <- table$cod_dqa
  colnames(mat_dist) <- table$cod_dqa
  mat_dist[upper.tri(mat_dist)] = NA
  net_distances<-reshape2::melt(mat_dist, na.rm=T)
  colnames(net_distances)<-c("item1","item2","net_dist")
  rownames(net_distances) <- c()
  net_distances<-filter(net_distances,item1!=item2)
```

Here's the result.

<img src="{{ site.url }}{{ site.baseurl }}/images/net_distance/net_dist_table.png" alt="linearly separable data">

