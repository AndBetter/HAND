# HAND
Create a large scale (possibly) global HAND raster. 
The repository contains a script to create the HAND 
index. 

The script requires as an data inputs: 
1. a DEM (possibly a DTM).
2. the "occurrence" dataset of the Global Surface Water Explorer - GSWE (https://global-surface-water.appspot.com/download)

The script allowes the user to cahnge the following input parameters: 
1. The upstream area required to initialize a drainage cell
2. The flow direction method used to identify drainage directions
3. A threshold frequency used for thresholding the  "occurrence" layer of the GSWE in order to identify water bodies
4. A threshold area for filtering small contiguous water bodies (i.e. contiguous water bodies smaller than the threshold will be removed and not considered as drainage cells)

