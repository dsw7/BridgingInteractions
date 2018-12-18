## Description

Algorithm for finding the closest protein surface coordinate relative to either a   
metal (**MT**) or bridging SD (**SD**) coordinate. The algorithm assumes that  
the input already contains both a bridging interaction and a metal.  

---  

## Step 1  
The input. I approximate a "protein surface" using a set of equidistant coordinates  
spaced a fixed radius from an origin (i.e. _on a circle!_).  
<img src = "https://github.com/dsw7/BridgingInteractions/blob/master/scalene-triangle/coordinate-to-surface/scalene_step1.png" width="400">

---  

## Step 2  
I project a series of vectors from the SD coordinate to all surface coordinates.  
<img src = "https://github.com/dsw7/BridgingInteractions/blob/master/scalene-triangle/coordinate-to-surface/scalene_step2.png" width="400">  
  
---  

## Step 3  
I project a series of vectors from the MT coordinate to all surface coordinates.  
<img src = "https://github.com/dsw7/BridgingInteractions/blob/master/scalene-triangle/coordinate-to-surface/scalene_step3.png" width="400">  

---  

## Step 4  
The set of coordinates yielding the minimum MT / surface distance are collected. The  
set of coordinates yielding the minimum SD / surface distance are also collected. Of the two pairs, the pair  
yielding the shortest distance is chosen. The surface coordinate yielding this shortest distance is  
labelled coordinate **SF**.  
<img src = "https://github.com/dsw7/BridgingInteractions/blob/master/scalene-triangle/coordinate-to-surface/scalene_step4.png" width="400">
