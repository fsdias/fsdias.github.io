---
title: "Distances between polygons with PosGIS and QGIS"
date: 2021-1-12
tags: [gis,PostGIS,QGIS]
header:
  image: "images/net_distance/banner_net_dist.png"
excerpt: "GIS, PostGIS, QGIS,ST_Distance,"
mathjax: "true"
---

In this post, I'm going to show you how to calculate the minimum distance between polygons (edge to edge) using QGIS and PostGIS. We will be working with a shapefile called "polygons.shp" that can be downloaded [here](https://github.com/fsdias/blog_data/tree/main/distances).

The first step is to load "polygons.shp" into your PostGIS database with the QGIS plugin DB Manager. Next, we run the following query that calculates the minimum distance between every polygon. My database is called "filipe" and the schema I'm using is called "riverscale"

```r
CREATE TABLE riverscale.distances AS
SELECT row_number() over () as id, ST_ShortestLine(g1.geom,g2.geom)::GEOMETRY(LINESTRING,3763)
AS geom, ST_Distance(g1.geom,g2.geom)
AS distance, g1.id AS id_from, g2.id
AS id_to FROM riverscale.polygons
AS g1, riverscale.polygons AS g2
```

This query creates a LINESTRING geometry which contains the shortest paths between pairs of features. Unfortunately this query create a geometry without primary key. We can solve this by running:

```r
ALTER TABLE riverscale.distances
ADD constraint distances_pkey PRIMARY KEY (id);
```

This is the output:

<img src="{{ site.url }}{{ site.baseurl }}/images/distances_postgis/distances_output.png" alt="">


**Filipe Dias**

