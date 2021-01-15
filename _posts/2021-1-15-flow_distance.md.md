---
title: "Calculate flow distances in R"
date: 2021-1-15
tags: [gis,r,Riverscale]
header:
  image: "images/net_distance/banner_net_dist.png"
excerpt: "GIS, River networks, flow distance, Riverscale"
mathjax: "true"
---

So far, I have explained how to calculate [network distances](https://fsdias.github.io/network_distance/) and [Euclidean distances](https://fsdias.github.io/euclidean_distance.md/) in a river network. Now, it's time to look at flow distance, which unlike the above distances considers the network's directionality. Flow distance measures the distance between sites that are connected by water flow.

We will be working with two shapefiles, one containing the river network and another the point location of the vegetation censuses. The points shapefile has to have a point located at the river's mouth. The purpose of this point is to identify the location where the water flows to. In the example shapefile, it is identified as "End" in the "id_code" column


Important notes about the shapefiles:

* the river shapefile can't have any topological problems. I recommend that you use QGIS and run the "Fix geometries" algorithm to fix any issues. It's available in the Processing toolbox.
* the points shapefile needs to have at least the following columns: id) a unique identifier, x) a column with the "x" coordinates, y) a column with "y" coordinates. 
* Use a metric coordinate reference system. Do not use WGS84 (EPSG 4326).

We will be working with "rivers.shp" and "points_flow.shp"([download](https://github.com/fsdias/blog_data/tree/main/flow_dist)). Map below.

<img src="{{ site.url }}{{ site.baseurl }}/images/flow_distance/map_flow_dist.png" alt="linearly separable data">


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

We start by checking the connectivity of the river network

```r
river = nt.connect(river)
```

and proceed by adding the vegetation sample points to the river network. We use "approach=2" to map each point to the nearest point on the network and add them as new nodes/vertices.

```r
  rivernel = points2network(ntdata=river,
                            pointsxy=coordinates(points),
                            ELComputed=TRUE,
                            approach=2, Detailed=TRUE, ea.prop=c(1,0))
```

Next, we convert the "rivernel" object into an "igraph" object".

```r
    rivg = nel2igraph(rivernel[[1]], rivernel[[2]], weight=rivernel[[8]])
```

Each vegetation sample point was added to the river network and received a unique identifier. We can access and save those unique identifiers by running:

```r
      ipoints = rivernel[[3]]
```

Now, to the fun part. We need to identiy all "paths" linking vegetation samples and the "End"" point (river mouth). Below, I exemplify this procedure. The red line shows that path that links samples 4,2, 1, 6 and the "End" point.


<img src="{{ site.url }}{{ site.baseurl }}/images/flow_distance/path.png" alt="flow paths">


We add each "path" to a list (conveniently) called "list".

```r
n <- c(1:length(ipoints))
  list<- list()
  for (i in n){
    v = sps$vpath[[i]]
    v1<-intersect(v,ipoints)
    list[[i]]<-v1
  }
```

Next, we select paths with two or more points with the following commands:

```r
list<-list %>% keep(function(x) length(x) > 1) 
```

Then, we need to create a table that links the unique identifier "id_code" wih the internal unique identifier that igraph generates for each point. I named the latter "cod_p2n".

```r
a<-data.frame(cod_p2n=seq_along(1:length(rivernel[[4]])),x=floor(rivernel[[4]]),y=floor(rivernel[[5]]))
id_num<-as.numeric(rivernel[[3]]) #vector with p2n codes for each added point
a<-filter(a, cod_p2n >= min(id_num)) #this works because new points are added to the end
df<-data.frame(x=floor(points@data$x),y=floor(points@data$y),id_code=points@data$id_code) #create data frame with xy coordinates and id_code
table<-merge(a,df,by=c("x","y")) #merged based on xy coordinates
table<-arrange(table,cod_p2n) #order points so that "End" is the last point
```

In the last step, we calculate the network distance between the vegetation samples of each "path"

```r
flowlist = list()
for (i in 1:length(list)){ #Number of
    path<-list[[i]]
    data_select<-filter(table, cod_p2n %in% path)
    data_select<-filter(data_select, id_code != "End")
    mat_dist<-distances(rivg,data_select$cod_p2n,data_select$cod_p2n) 
    rownames(mat_dist) <-data_select$id_code
    colnames(mat_dist) <- data_select$id_code
    mat_dist[upper.tri(mat_dist)] = NA
    flow_distances<-reshape2::melt(mat_dist, na.rm=T)
    colnames(flow_distances)<-c("item1","item2","flow_dist")
    flow_distances<-filter(flow_distances,item1!=item2)
    flowlist[[i]] <- distinct(flow_distances) #don't keep repeated comparisons
  }
```

and create with a data.frame.

```r
  output = do.call(rbind, flowlist)
```

Here's what it looks like.

<img src="{{ site.url }}{{ site.baseurl }}/images/flow_distance/flow_dist_table.png" alt="flow paths">

*Filipe Dias
