# !!!! usare pgrid script in this folder (pgrid.py deve essere nella stessa cartella
# dello script principale)


'''
questa versione "PROVA" rispetto all' altra ritaglia l'output alla dimensione
effettiva della tile, rimuovendo i nan

fix scripts in pysheds as described in https://github.com/mdbartos/pysheds/pull/268
'''

#Setup Environment
from matplotlib import pyplot as plt
import os
import numpy as np
#from functools import partial
import warnings #Suppress warnings on occasion
#import time     #to time HAND nan-fill duration
from pathlib import Path
from pysheds.grid import Grid
from affine import Affine
import rasterio
import pyproj
#import fiona
#import shapely
#import geopandas as gpd
#import astropy
#import astropy.convolution
#from tqdm.auto import tqdm
from pysheds.view import Raster, ViewFinder
import cv2 
import rasterio.mask
from rasterio.warp import reproject, Resampling
import multiprocessing
from rasterio.windows import Window
from rasterio.transform import rowcol, xy
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds
from scipy.ndimage import label
#from scipy.interpolate import griddata
from scipy.spatial import cKDTree

################
# FUNCTIONS ####
################

# =============================================================================
# def fiona_read_vectorfile(vectorfile, get_property=None):
#     """shapes=fiona_read_vectorfile(vectorfile, get_property=None)
#        shapes, props=fiona_read_vectorfile(vectorfile, get_property='Property_Name')
#        Returns a list of shapes (and optionally properties) using fiona.
#        
#        vectorfile: any fiona compatible vector file. 
#        get_property: String for the property to be read. 
#        shapes: List of vector "geometry"
#        props:  List of vector "properties"
#     """
#     with fiona.open(vectorfile, "r") as shpf:
#         shapes   = [ feature["geometry"] for feature in shpf ]
#         print(f"Number of shapes loaded: {len(shapes)}")
#         if get_property is not None:
#             props = [ feature["properties"][get_property] for feature in shpf ]
#             return shapes, props
#         else:
#             return shapes    
# 
# 
# 
# =============================================================================


# =============================================================================
# def intersect(shapes, polygon, properties=None):
#     """
#     polygons=intersect(shapes, polygon, properties=None)
#     Returns polygons from multiple 'geometries' read by fiona.
#     
#     shapes: shapes returned by fiona_read_vectorfile()
#     polygon: a single polygon to intersect with shapes
#     properties: If not none, returns the property value instead of polygon geometry.
#     """
#     #first loop to split multi polygons to single polygons
#     polygons=[]
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")    
#         for k,shape in enumerate(tqdm(shapes)):
#             if shape['type']=='MultiPolygon':
#                 for l,p in enumerate(shape['coordinates']):
#                     s=shapely.geometry.Polygon(p[0])
#                     if polygon.intersects(s) and properties is None:            
#                         polygons.append(s)
#                     elif polygon.intersects(s) and properties is not None:
#                         if np.isscalar(properties[k]):
#                             polygons.append(properties[k])
#                         else:
#                             polygons.append(properties[k][l])
#                     
#             elif shape['type']=='Polygon':
#                 s=shapely.geometry.Polygon(shape['coordinates'][0])
#                 if polygon.intersects(s) and properties is None:
#                     polygons.append(s)
#                 elif polygon.intersects(s) and properties is not None:
#                     polygons.append(properties[k])
#     return polygons
# =============================================================================


# =============================================================================
# def reproject_(vector_file, output_crs, output_file=None):
#     """
#     output_file=reproject_(vector_file, output_crs, output_file=None)
#     Reprojects a given vector file to another reference frame (CRS). 
#     vector_file: Any vector file that can be opened with GeoPandas
#     output_crs: A rasterio opened crs (e.g. dem.crs)
#     output_file: if not defined, defaults to vector_file[:-4]+'_warp.shp'. 
#     """
#     v=gpd.GeoDataFrame.from_file(vector_file)
#     warp=v.to_crs(output_crs)
#     if output_file is None:
#         output_file=vector_file[:-4]+'_warp.shp'
#     warp.to_file(output_file)
#     return output_file
# =============================================================================


# =============================================================================
# def transform_shape(shape, s_srs='epsg:4326', t_srs='epsg:4326'):
#     transformation=partial(
#                pyproj.transform,
#                pyproj.Proj(init=s_srs), #source coordinate system
#                pyproj.Proj(init=t_srs)) #destination coordinate system
#     return shapely.ops.transform(transformation, shape)
# =============================================================================



# =============================================================================
# def xy2coord(x,y,gT):
#     '''
#     lon,lat=xy2coord(x,y,geoTransform)
#     converts pixel index to position based on geotransform.
#     '''
#     coord_x=gT[0] + x*gT[1] + y*gT[2]
#     coord_y=gT[3] + x*gT[4] + y*gT[5]
#     return coord_x, coord_y
# =============================================================================


# =============================================================================
# def get_bounding_box_rasterio(raster_path):
#     with rasterio.open(raster_path) as src:
#         # Get the affine transform, width, and height of the raster
#         transform = src.transform
#         width = src.width
#         height = src.height
# 
#         # Calculate the coordinates of the four corners
#         # (left, bottom, right, top)
#         left, top = transform * (0, 0)
#         right, bottom = transform * (width, height)
# 
#         # Return the bounding box
#         return left, bottom, right, top
# =============================================================================


def pixel_size_WGS84_array ( array, trans ):
    '''
    computes the area of each pixel in a wgs84 reference system
    and the average pixel height and width
    '''
    
    #lon_min = transform[2]
    lat_max  = trans[5]
    #d_lon    = transform[0]
    d_lat    = trans[4]
    
    #Earth radius in m
    R = 6378000
    
    rows, cols = array.shape
    
    #lon = np.linspace(lon_min, lon_min + cols*d_lon ,  cols) 
    lat = np.linspace(lat_max, lat_max + rows*d_lat ,  rows)

    #lon_grid, lat_grid = np.meshgrid(lon, lat)
    d_y = -np.radians(d_lat) * R
    d_x  = d_y * np.cos(np.radians(lat))

    d_y_mean = np.mean(d_y)
    d_x_mean = np.mean(d_x)
    pixel_area = np.ones(shape=(rows,cols)).astype(np.float16) * d_x[:, np.newaxis] * d_y

    return pixel_area, d_y_mean, d_x_mean  #lat_grid.astype(np.float16), lon_grid.astype(np.float16)


# =============================================================================
# def fill_nan(arr):
#     """
#     filled_arr=fill_nan(arr)
#     Fills Not-a-number values in arr using astropy. 
#     """    
#     kernel = astropy.convolution.Gaussian2DKernel(x_stddev=3) #kernel x_size=8*stddev
#     arr_type=arr.dtype          
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")
#         while np.any(np.isnan(arr)):
#             arr = astropy.convolution.interpolate_replace_nans(arr.astype(float), kernel, convolve=astropy.convolution.convolve)
#     return arr.astype(arr_type) 
# 
# 
# def fill_nan_bis(arr, mask=None):
#     """
#     Fills Not-a-number values in `arr` using astropy, but only at positions specified by `mask`.
#     
#     Parameters:
#     - arr: NumPy array containing `NaN` values to be filled.
#     - mask: NumPy boolean array of the same shape as `arr`. `True` indicates cells that should be filled.
#     
#     Returns:
#     - filled_arr: NumPy array with specified `NaN` values filled.
#     """
#     kernel = astropy.convolution.Gaussian2DKernel(x_stddev=3)  # Kernel with x_size=8*stddev
#     arr_type = arr.dtype
#     
#     if mask is None:
#         mask = np.isnan(arr)
#     
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")
#         
#         # Create a copy of the array to avoid modifying the original
#         arr_to_fill = np.copy(arr).astype(float)
#         
#         # Loop to fill only masked `NaN` values
#         while np.any(np.isnan(arr_to_fill) & mask):
#             arr_to_fill = astropy.convolution.interpolate_replace_nans(
#                 arr_to_fill, kernel, convolve=astropy.convolution.convolve
#             )
#             # Reapply mask to ensure only specified `NaN` positions are filled
#             arr_to_fill[~mask] = arr[~mask]
#     
#     return arr_to_fill.astype(arr_type)
# 
# 
# def fill_nan_based_on_DEM(arr, dem):
#     """
#     filled_arr=fill_nan_based_on_DEM(arr, dem)
#     Fills Not-a-number values in arr using astropy. 
#     """    
#     hond = dem - arr; #height of nearest drainage 
#     kernel = astropy.convolution.Gaussian2DKernel(x_stddev=3) #kernel x_size=8*stddev
#     arr_type=hond.dtype          
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")
#         while np.any(np.isnan(hond)):
#             hond = astropy.convolution.interpolate_replace_nans(hond.astype(float), kernel, convolve=astropy.convolution.convolve)
#     my_mask=np.isnan(arr)
#     arr[my_mask]=dem[my_mask]-hond[my_mask]
#     return arr.astype(arr_type) 
# 
# 
# # Fill NaN values with the mean of their neighbors
# def fill_nodata_with_mean(raster, valid_pixels, valid_pixels_nodata):
#     """
#     Fills NaN values in a raster by computing the mean of their neighboring pixels.
#     
#     This function identifies NaN values in the input raster and replaces them
#     with the mean of their neighboring values. It selectively fills NaNs based
#     on the `valid_pixels` array and the specified `valid_pixels_nodata` value.
#     
#     Parameters:
#     - raster: NumPy array
#       A 2D array representing the raster data with potential NaN values to fill.
#     
#     - valid_pixels: NumPy array
#       A 2D array of the same shape as `raster` indicating valid pixels. It can
#       be used to specify positions that are considered valid for filling.
#     
#     - valid_pixels_nodata: Scalar
#       A value indicating 'no data' in the `valid_pixels` array. If NaN, it indicates
#       that NaNs in `valid_pixels` should be ignored when filling `raster`.
#     
#     Returns:
#     - raster: NumPy array
#       The modified raster with NaN values filled by the mean of their neighbors.
#       If all neighbors of a NaN are NaN themselves, the function assigns a value
#       of 0 to that position.
#     
#     Notes:
#     - The function uses a moving window approach to compute the mean of neighboring
#       values. It handles edge cases where pixels may have fewer neighbors.
#     - Ensure that `valid_pixels` and `raster` have the same shape to avoid errors.
#     """   
#     if np.isnan(valid_pixels_nodata):
#         nan_indices = np.argwhere((np.isnan(raster)) & (~np.isnan(valid_pixels)))
#     else:
#         nan_indices = np.argwhere((np.isnan(raster)) & (valid_pixels != valid_pixels_nodata))
#     
#     for y, x in nan_indices:
#         neighbors = raster[max(0, y-2):y+3, max(0, x-2):x+3]
#         if np.isnan(neighbors).all():
#             # Handle the case where all neighbors are NaN
#             raster[y, x] = 0
#         raster[y, x] = np.nanmean(neighbors)
#     return raster
# 
# # Fill NaN values with the mean of their neighbors
# def fill_nodata_with_mean_mask(raster, where_to_fill_mask):
# 
#     nan_indices = np.argwhere(where_to_fill_mask>0)
#     
#     for y, x in nan_indices:
#         neighbors = raster[max(0, y-2):y+3, max(0, x-2):x+3]
#         if np.isnan(neighbors).all():
#             # Handle the case where all neighbors are NaN
#             raster[y, x] = 0
#         raster[y, x] = np.nanmean(neighbors)
#     return raster
# 
# 
# def fill_nan_cv2(image, nodata_value=np.nan, method=cv2.INPAINT_TELEA): # method=cv2.INPAINT_NS ; method=cv2.INPAINT_TELEA
#     """
#     Inpaint a grayscale image by filling nodata values.
# 
#     Parameters:
#     - image: NumPy array of the grayscale image with nodata values.
#     - nodata_value: The value used to represent nodata in the image.
#     - method: Inpainting method, either cv2.INPAINT_TELEA or cv2.INPAINT_NS.
# 
#     Returns:
#     - inpainted_image: The image with nodata values filled.
#     """
#     # Create a mask where nodata values are marked
#     if np.isnan(nodata_value):
#         mask = (np.isnan(image)).astype('uint8')
#     else:
#         mask = (image == nodata_value).astype('uint8')
# 
#     # Use OpenCV inpainting to fill the gaps
#     inpainted_image = cv2.inpaint(image, mask, inpaintRadius=10, flags=method)
# 
#     return inpainted_image
# 
# =============================================================================

def fill_missing_with_nearest_neighbor(raster,where_to_fill_mask, k=10):
    """
    Fill missing values in a raster using the nearest neighbor approach.
    
    Parameters:
    raster (np.ndarray): A 2D numpy array with missing values indicated by np.nan.
    k (int): Number of nearest neighbors to consider.
    Returns:
    np.ndarray: The raster with missing values filled.
    """
    
    where_to_fill_mask = where_to_fill_mask.astype(np.uint8)
    
    # Find indices of missing values
    missing_indices = np.argwhere(where_to_fill_mask == 1)
    
    # Dilate the mask to find the border areas
    where_to_fill_mask_dilate = cv2.morphologyEx(where_to_fill_mask, cv2.MORPH_DILATE, np.ones((3, 3), dtype=np.uint8))
    borders = where_to_fill_mask_dilate - where_to_fill_mask

    # Find indices and values of non-missing values on the border
    non_missing_indices = np.argwhere((borders == 1) & (~np.isnan(raster)))
    non_missing_values  = raster[(borders == 1) & (~np.isnan(raster))]

    # Use a KDTree for efficient nearest-neighbor search
    tree = cKDTree(non_missing_indices)

    # Iterate over missing values and fill them
    for index in missing_indices:
        # Query the k-nearest neighbors
        distances, nearest_indices = tree.query(index, k=k)
        nearest_values = non_missing_values[nearest_indices]

        # Compute the mean of the nearest values
        mean_value = np.nanmean(nearest_values)
        
        # Fill the missing value with the mean
        raster[tuple(index)] = mean_value

    return raster




def save_tiff(data_array, crs, trans, output_folder, file_name='out.tif', compression='lzw', dtype=None):
    """
    Save a numpy array as a GeoTIFF file with specified properties.

    Parameters:
    data_array (np.ndarray): The array to save.
    crs (str or dict): Coordinate reference system.
    trans (Affine): Affine transform for the raster.
    output_folder (str): Directory where the output file will be saved.
    file_name (str): Name for the output file.
    compression (str): Compression type (e.g., 'lzw', 'deflate').
    dtype (str): Data type for the output file (e.g., 'uint8', 'int16', 'float32').
    """
    # Use the provided dtype or infer from the data array
    dtype = dtype or data_array.dtype

    # Define the metadata
    metadata = {
        'driver': 'GTiff',  # GeoTIFF format
        'height': data_array.shape[0],
        'width': data_array.shape[1],
        'count': 1,  # Number of bands
        'dtype': dtype,
        'crs': crs,
        'transform': trans,
        'compress': compression  # Compression type (e.g., 'lzw', 'deflate')
    }

    # Construct the full output file path
    output_file = f'{output_folder}/{file_name}.tif'

    # Save the array as a raster file
    with rasterio.open(output_file, 'w', **metadata) as dst:
        dst.write(data_array, 1)  # Write the array to the first band

        
     
def bring_raster_to_same_grid(source_array,source_crs,source_tf, target_crs,target_tf,target_width,target_height ):

    # Prepare the destination array
    output_array = np.empty((target_height, target_width), dtype=source_array.dtype)
    
    # Reproject and resample
    raster_transformed =reproject(
                                    source        = source_array,
                                    destination   = output_array,
                                    src_transform = source_tf,
                                    src_crs       = source_crs,
                                    dst_transform = target_tf,
                                    dst_crs       = target_crs,
                                    resampling=Resampling.nearest  
                                    )
        
    return raster_transformed
       

    
def find_closest_divisor(x, y):
    # Find all divisors of x
    divisors = [i for i in range(1, x + 1) if x % i == 0]
    
    # Find the divisor closest to y
    closest_divisor = min(divisors, key=lambda d: abs(d - y))
    
    return closest_divisor


# =============================================================================
# def create_raster_windows(dem_height, dem_width,  buffer_size_window , approx_tile_size):
#     
#     tile_width  = find_closest_divisor(dem_width,  approx_tile_size)
#     tile_height = find_closest_divisor(dem_height, approx_tile_size)
#     
#     
#     windows=[]
#     window_ids=[]
#     windows_ids_progressive = []
#     ind=0
#     row_id = 0
#     # create tiles windows
#     for top in range(0, dem_height, tile_height):
#         col_id = 0
#         
#         for left in range(0, dem_width, tile_width):
#             # Define the window (tile) boundaries
#             windows.append(Window(left, top, tile_width, tile_height))
#             window_ids.append(f'row_{row_id}_col_{col_id}')
#             col_id+=1
#             windows_ids_progressive.append(ind)
#             ind+=1
#             
#         row_id+=1    
#       
#     # Create the windows extended with buffer
#     windows_extended = []
#     
#     # Create extended windows
#     for window in windows:
#         # Calculate the extended window boundaries
#         left_extended   = max(0, window.col_off - buffer_size_window)
#         top_extended    = max(0, window.row_off - buffer_size_window)
#         right_extended  = min(dem_width, window.col_off + window.width + buffer_size_window)
#         bottom_extended = min(dem_height, window.row_off + window.height + buffer_size_window)
#     
#         # Calculate width and height for the extended window
#         width_extended  = right_extended - left_extended
#         height_extended = bottom_extended - top_extended
#     
#         # Define the extended window
#         extended_window = Window(left_extended, top_extended, width_extended, height_extended)
#         windows_extended.append(extended_window)
#         
# 
#     return windows, windows_extended, window_ids, windows_ids_progressive
# =============================================================================



def create_raster_windows_fixed(dem_height, dem_width,  buffer_size_window , tile_size):
    
    
    windows=[]
    window_ids=[]
    windows_ids_progressive = []
    ind=0
    row_id = 0
    # create tiles windows
    for top in range(0, dem_height, tile_size):
        col_id = 0
        
        for left in range(0, dem_width, tile_size):
            # Define the window (tile) boundaries
            windows.append(Window(left, top, min(tile_size, dem_width - left), min(tile_size, dem_height - top)))
            window_ids.append(f'row_{row_id}_col_{col_id}')
            col_id+=1
            windows_ids_progressive.append(ind)
            ind+=1
            
        row_id+=1    
      
    # Create the windows extended with buffer
    windows_extended = []
    
    # Create extended windows
    for window in windows:
        # Calculate the extended window boundaries
        left_extended   = max(0, window.col_off - buffer_size_window)
        top_extended    = max(0, window.row_off - buffer_size_window)
        right_extended  = min(dem_width, window.col_off + window.width + buffer_size_window)
        bottom_extended = min(dem_height, window.row_off + window.height + buffer_size_window)
    
        # Calculate width and height for the extended window
        width_extended  = right_extended - left_extended
        height_extended = bottom_extended - top_extended
    
        # Define the extended window
        extended_window = Window(left_extended, top_extended, width_extended, height_extended)
        windows_extended.append(extended_window)
        

    return windows, windows_extended, window_ids, windows_ids_progressive


def identify_tiles_to_process(list_tiles_to_process):
    window_ids_array  = np.array(window_ids)
    
    windows_id_progressive_subset=[]; window_ids_subset=[]; windows_extended_subset=[]; windows_subset=[]
    
    for tile_to_process in list_tiles_to_process:
        idx = np.where(window_ids_array == tile_to_process)
        windows_id_progressive = idx[0][0]
        window_id = window_id=window_ids[windows_id_progressive]
        window_extended=windows_extended[windows_id_progressive]
        window         =windows[windows_id_progressive]
    
        windows_id_progressive_subset.append(windows_id_progressive)
        window_ids_subset.append(window_id)
        windows_extended_subset.append(window_extended)
        windows_subset.append(window)
    
    return windows_id_progressive_subset, window_ids_subset, windows_extended_subset, windows_subset


def remove_small_patches(binary_array, weights_array, threshold_area):
    
    """
    Remove small connected components from a binary array based on their weighted size.

    This function labels connected components in the binary array, calculates their
    weighted sizes using the provided weights array, and retains only those components
    whose weighted size is greater than or equal to the specified threshold.

    Parameters:
    binary_array (numpy.ndarray): A 2D binary array (uint8 or bool) where connected
                                 components are defined. Typically, this array contains
                                 values of 0 (background) and 1 (foreground).
    weights_array (numpy.ndarray): A 2D array of the same shape as binary_array,
                                   containing weight values corresponding to each pixel
                                   in binary_array.
    threshold_area (float or int): The minimum weighted size a connected component
                                   must have to be retained. Components with a weighted
                                   size smaller than this threshold will be removed.

    Returns:
    numpy.ndarray: A filtered 2D binary array with small components removed, retaining
                   only the components whose weighted size is greater than or equal to
                   the threshold_area.
    """

    # Label connected components
    labeled_array, num_features = label(binary_array)

    # Calculate the weighted size of each component
    weighted_sizes = np.bincount(labeled_array.ravel(), weights=weights_array.ravel())
    
    # Create a mask for components that are larger than the minimum weighted size
    mask = weighted_sizes >= threshold_area
    
    # Remove smaller components
    filtered_array = mask[labeled_array]
    filtered_array[binary_array==0] = 0
    
    return filtered_array.astype(np.uint8)


# =============================================================================
# dem_array = dem.read()[0,:,:].astype(np.float32)
# dem_array[dem_array==dem.nodata]=np.nan
# dem_array = dem.offsets[0] + dem.scales[0]* dem_array
# =============================================================================

# =============================================================================
# dem_array = terreno
# dem_gT = tf
# dem_proj4=dem.crs
# catch_mask=not_mask
### dtm_aux = terreno_dtm_aux_resampled
# gswe_array  = data_gswe_window_resampled
# =============================================================================


def calculate_hand(dem_array, 
                   #dtm_aux,
                   gswe_array,
                   dem_gT, 
                   dem_proj4 , 
                   catch_mask=None, 
                   acc_thresh=10, 
                   routing = 'dinf', 
                   minimum_area_threshold_water_bodies = 10,
                   geographic_coordinates = False):
    """
    Computes the Height Above Nearest Drainage (HAND) for a given Digital Elevation Model (DEM).

    This function processes a DEM to compute the HAND, which is the vertical distance from any
    point in the landscape to the nearest drainage channel. It includes several steps such as
    identifying sea coastlines, water bodies, and conditioning the DEM to handle pits, depressions,
    and flats.

    Parameters:
    dem_array (numpy.ndarray): A 2D numpy array representing the Digital Elevation Model (DEM).
    gswe_array (numpy.ndarray): Global Surface Water Explorer (GSWE) array indicating water presence.
    dem_gT (tuple or affine.Affine): GeoTransform of the input DEM, either as a tuple or an Affine object.
    dem_proj4 (str): Proj4 string representing the coordinate reference system of the DEM.
    catch_mask (numpy.ndarray, optional): A binary mask covering the catchment (1 for catchment, 0 for no catchment).
                                     If not provided, the entire DEM is evaluated.
    acc_thresh (int, optional): Accumulation threshold in square kilometers. Default is 1.
    routing (str, optional): Type of routing used to delineate the network. Options are 'd8', 'dinf', or 'mfd'. Default is 'dinf'.
    geographic_coordinates (bool, optional): Flag denoting if the reference system is geographic (True) or projected (False). Default is False.

    Returns:
    numpy.ndarray: A 2D numpy array representing the Height Above Nearest Drainage (HAND) for the input DEM.
    """

    #establish no data
    dem_array_nodata = np.float32(np.nan)


    #assigns geotransforms
    if type(dem_gT)==Affine:
        aff=dem_gT
    else:
        aff=Affine.from_gdal(*tuple(dem_gT))
        
    # what to do if there is no catchment mask 
    if catch_mask is None:
        catch_mask=np.ones(dem_array.shape, dtype=np.bool)
      
       
    #computes the area of the pixels
    if geographic_coordinates == True:   
        pixel_areas, d_y_mean, d_x_mean = pixel_size_WGS84_array ( dem_array, trans=aff  )
    else:
        pixel_areas= np.ones_like(dem_array, dtype=np.float32) * aff[0]*np.abs(aff[4])
        d_y_mean =  np.abs(aff[4])   
        d_x_mean =  aff[0]
     
        
    # IDENTIFY SEA COASTLINES

    #identify the sea based on elevation and slope
    threshold_slope   = 0.001
    
    #elev = dtm_aux# dem_array #dtm_aux_raster[0]
    gradient_i,gradient_j=np.gradient(dem_array,d_y_mean,d_x_mean)
    slope = np.sqrt(gradient_i**2+gradient_j**2)
    flat_areas  = slope < threshold_slope

    #flat_areas = cv2.medianBlur(flat_areas.astype(np.uint8), ksize=11)
    
    flat_areas[(dem_array>1) | (dem_array<-1)] = 0
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    sea = cv2.morphologyEx(flat_areas.astype(np.uint8), cv2.MORPH_CLOSE, kernel, 	iterations = 3 )
    sea = cv2.morphologyEx(sea, cv2.MORPH_DILATE, np.ones((3,3),dtype=np.uint8) )
       
     
    #because GEDTM30 has noise and biased values on the sea   
    dem_array[sea==1]=0
        
    
    # IDENTIFY WATER BODIES WITH GSWE
    gswe_array_frequency_threshold   = 90 # percent of time of water presence in the GSWE
    gswe_array_filtered = remove_small_patches(gswe_array>gswe_array_frequency_threshold, pixel_areas, minimum_area_threshold_water_bodies *1e6) 
    gswe_array_filtered = cv2.morphologyEx(gswe_array_filtered, cv2.MORPH_DILATE, np.ones((3,3),dtype=np.uint8) )

    #merge the water bodies identified with GSWE and topography 
    water_bodies = np.maximum(sea, gswe_array_filtered)    

    #create raster objects for inputs to pysheds
    raster      = Raster(dem_array,    viewfinder=ViewFinder(shape=dem_array.shape, affine = aff, crs=dem_proj4, mask=catch_mask, nodata= dem_array_nodata)) 
    pixel_areas = Raster(pixel_areas,  viewfinder=ViewFinder(shape=dem_array.shape, affine = aff, crs=dem_proj4, mask=catch_mask, nodata= dem_array_nodata))
    gridd = Grid.from_raster(raster)


    #DEM CONDITIONING
    pit_filled_dem = gridd.fill_pits        (raster,         nodata = raster.nodata, nodata_out=raster.nodata)
    flooded_dem    = gridd.fill_depressions (pit_filled_dem, nodata = raster.nodata, nodata_out=raster.nodata)
    inflated_dem   = gridd.resolve_flats    (flooded_dem,    nodata = raster.nodata, nodata_out=raster.nodata)


    #repetes to reduce singularities
    repetitions = 3
    for i in range(repetitions):
        pit_filled_dem = gridd.fill_pits        (inflated_dem   , nodata = raster.nodata, nodata_out=raster.nodata)
        flooded_dem    = gridd.fill_depressions (pit_filled_dem , nodata = raster.nodata, nodata_out=raster.nodata)
        inflated_dem   = gridd.resolve_flats    (flooded_dem    , nodata = raster.nodata, nodata_out=raster.nodata)





    if routing == 'd8':
        fdir_d8 = gridd.flowdir(inflated_dem, routing='d8' , nodata = inflated_dem.nodata, nodata_out = np.int64(0) )
        acc_d8  = gridd.accumulation(fdir_d8, weights= pixel_areas,  nodata = np.int64(0), nodata_out = np.float64(np.nan) )
        net_d8  = ((acc_d8> acc_thresh*1e6) + water_bodies)>0
        hand    = gridd.compute_hand(fdir_d8, inflated_dem, net_d8, routing='d8',nodata=np.float64(np.nan), nodata_out=np.float64(np.nan))

    elif routing == 'dinf':
        fdir_dinf = gridd.flowdir(inflated_dem, routing='dinf', nodata = inflated_dem.nodata, nodata_out=np.float64(np.nan)).astype(np.float64)
        acc_dinf  = gridd.accumulation(fdir_dinf, weights= pixel_areas, routing='dinf', nodata=np.float64(np.nan), nodata_out=np.float64(np.nan)).astype(np.float64)
        net_dinf  = ((acc_dinf> acc_thresh*1e6) + water_bodies)>0
        hand      = gridd.compute_hand(fdir_dinf, inflated_dem, net_dinf, routing='dinf', nodata=np.float64(np.nan), nodata_out=np.float64(np.nan)).astype(np.float64)


    elif routing == 'mfd':
        fdir_mfd = gridd.flowdir(inflated_dem, routing='mfd', nodata_out=np.float64(np.nan)).astype(np.float64)
        acc_mfd  = gridd.accumulation(fdir_mfd, weights= pixel_areas, routing='mfd', nodata_out=np.float64(np.nan)).astype(np.float64)
        net_mfd  = ((acc_mfd> acc_thresh*1e6) + water_bodies)>0
        hand     = gridd.compute_hand(fdir_mfd,inflated_dem, net_mfd, routing='mfd',nodata=np.float64(np.nan), nodata_out=np.float64(np.nan)).astype(np.float64)



    maschera = ( np.isnan(hand)) & (catch_mask)
    hand     = fill_missing_with_nearest_neighbor(hand, where_to_fill_mask=maschera )


    hand[np.bitwise_not(catch_mask)] = np.nan
    

    return hand



#get geotransform of a window
def get_transform(raster_path, window):
    
    with rasterio.open(raster_path) as raster:
        transform    = raster.window_transform(window)

    return transform



def extract_window_raster(dtm_path, window):
    """
    Extracts a windowed portion of a raster file and applies necessary transformations.

    This function reads a specified window from a Digital Elevation Model (DEM) raster file,
    applies any necessary scaling and offset transformations, and handles no-data values.

    Parameters:
    dtm_path (str): The file path to the Digital Elevation Model (DEM) raster file.
    window (tuple): A tuple defining the window to extract from the raster. Typically, it is a
                    4-tuple (col_off, row_off, width, height) where:
                    - col_off is the column offset,
                    - row_off is the row offset,
                    - width is the number of columns to read,
                    - height is the number of rows to read.

    Returns:
    tuple: A tuple containing:
        - terreno (numpy.ndarray): A 2D numpy array representing the extracted window of the DEM.
        - tf (affine.Affine): The affine transformation corresponding to the window.

    The function handles integer data types by converting them to float32, replaces no-data values
    with NaN, and applies any scaling and offset transformations specified in the raster metadata.
    """

    with rasterio.open(dtm_path) as dem:
        # Read the data over the specified window
        terreno = dem.read(1, window=window, boundless=True, fill_value=dem.nodata)
        tf      = dem.window_transform(window)
        nodata  = dem.nodata
        scale   = dem.scales[0]
        offset  = dem.offsets[0]
    

    if terreno.dtype.kind == 'i':
        terreno=terreno.astype(np.float32)
    
    if not nodata==None:
        if not np.isnan(nodata):
            terreno[terreno==nodata]=np.nan
    
    if scale != 1 or offset !=0:
        terreno = offset + scale * terreno
    
    return terreno, tf






#############################################################################
#############################################################################
## START INPUTS
#############################################################################
#############################################################################

#WHICH FLOW ACCUMULATION ROUTING TO USE
routing = 'mfd' # 'd8', 'dinf', 'mfd'

#UPSTREAM AREA THRESHOLD FOR INITIALIZATION OF DRAINAGE 
acc_thresh= 1 #   # in skm # None # This sets how large of an accumulation area is used for HAND. If too small, we get a very fine river network, which can be noisy. If too high, we get a very smooth HAND...                 

#ALL CONTIGUOUS WATER BODIES SMALLER THAN THIS (IN SKM) ARE REMOVED (I.E. NOT COSIDERED AS DRAINAGE)
minimum_area_threshold_water_bodies   = 1 # in skm  

# PADDING ADDED AROUND THE BASIN POLYGON (IF USED)
#pad_width=200 # Padding applied to the hydrobasins polygons for HAND processing. 

#WHERE TO SAVE
output_folder = f'/home/bettand/Desktop/bettand/WETLANDS/friuli_prove/output_tiles_routing_{routing}_upstresh_{acc_thresh}_wbodiesfill_{minimum_area_threshold_water_bodies}'

#SIZE OF THE TILES IN OUTPUT
approx_tile_size    = 10000

#ADDITIONAL PROCESSING EXTENT AROUND THE TILE TO AVOID BORDER EFFECT IN HAND COMPUTATIONS
buffer_size_window  = 500



#Select DTM
#dtm_path = Path('/home/bettand/Desktop/eos_bettand/Data/Topography/GEDTM30/legendtm_rf_30m_m_s_20000101_20231231_go_epsg.4326_v20250130.tif')
#dtm_path = Path('/home/bettand/Desktop/eos_bettand/Data/Topography/GEDTM30/GEDTM30_tiled/output.tif')
#dtm_path = Path('/home/bettand/Desktop/eos_bettand/Data/Topography/GEDTM30/GEDTM30_tiled/output_small_proj_float.tif')
#dtm_path = Path('/home/bettand/Desktop/bettand/WETLANDS/friuli_prove/friuli_proj.tif')
#dtm_path = Path('/home/bettand/Desktop/eos_bettand/Data/Topography/GEDTM30/GEDTM30_tiled/output_small.tif')
#dtm_path = Path('/home/bettand/Desktop/bettand/WETLANDS/friuli_prove/cratia_DTM_extract.tif')
#dtm_path = Path('/home/bettand/Desktop/eos_bettand/Data/Topography/GEDTM30/GEDTM30_WGS84_native_EUROPE.tif')
#dtm_path = Path('/home/bettand/Desktop/eos_bettand/Data/Topography/GEDTM30/GEDTM30_WGS84_native_EUROPE.tif')

#dtm_path = Path('/home/bettand/Desktop/eos_bettand/Data/Topography/GEDTM30/GEDTM30_projected_30m_robinson_bilinear_gdal.tif')


dtm_path = Path('/home/bettand/Desktop/eos_bettand/Data/Topography/GEDTM30//NEW_GEDTM30_projected_30m_robbinson_bilinear_gdal.tif')




#PATH TO GSWE
glob_s_w_explorer_path = Path('/home/bettand/Desktop/eos/jeodpp/data/base/Hydrography/GLOBAL/JRC/GlobalSurfaceWater/Aggregated/VER5-0/Data/occurrence/VRT/gsw_occurrence.vrt')



#uses FABDEM to correct sealevel
#dtm_aux_path = Path('/home/bettand/Desktop/eos_bettand/Data/Fabdem/World/World_complete_FABDEM_wgs4_highcompression.tif')

#USES COPDEM TO CORRECT SEALEVEL
#tm_aux_path = Path('/home/bettand/Desktop/eos/jeodpp/data/base/Elevation/GLOBAL/COP-DEM/GLO-30-DGED/LATEST/Data/VRT/Copernicus_DSM_10.vrt')



##################################
##################################
## END INPUTS
##################################
##################################





os.makedirs(output_folder, exist_ok=True)

#get metadata GEDTM
with rasterio.open(dtm_path) as dem:
    dem_width  = dem.width
    dem_height = dem.height
    dem_gT    = dem.transform
    dem_proj4 = dem.crs
    dem_dtype = dem.dtypes[0]
    #assets if input dem is in geographic coords or not
    geographic_coordinates = dem.crs.is_geographic



# =============================================================================
# #get metadata FABDEM
# with rasterio.open(dtm_aux_path) as dtm_aux:
#     dtm_aux_width = dtm_aux.width
#     dtm_aux_height= dtm_aux.height
#     dtm_aux_gT = dtm_aux.transform
#     dtm_aux_proj4 = dtm_aux.crs
#     dtm_aux_dtype = dtm_aux.dtypes[0]  
#     geographic_coordinates_dtm_aux = dtm_aux.crs.is_geographic
# =============================================================================
    
    
#get metadata Surface Water Explorer 
with rasterio.open(glob_s_w_explorer_path) as gswe:
    gswe_width = gswe.width
    gswe_height= gswe.height
    gswe_gT = gswe.transform
    gswe_proj4 = gswe.crs
    gswe_dtype = gswe.dtypes[0]  
    geographic_coordinates_gswe = gswe.crs.is_geographic  




#create tiles
windows, windows_extended, window_ids, windows_ids_progressive = create_raster_windows_fixed(dem_height,
                                                                                             dem_width,  
                                                                                             buffer_size_window = buffer_size_window , 
                                                                                             tile_size = approx_tile_size)

number_windows = len(windows)



#define a subset of tiles to be processed (optional)
list_tiles_to_process = ['row_1_col_88']



windows_id_progressive_subset, window_ids_subset, windows_extended_subset, windows_subset  = identify_tiles_to_process(list_tiles_to_process)


# if TRUE processes the tiles in the range, otherwise precesses the subsed of tiles listed above
if True:
    tiles_from = 0
    tiles_to   = None
    
    process_subset = False
    
    windows_ids_progressive = windows_ids_progressive[tiles_from:tiles_to]
    window_ids              = window_ids[tiles_from:tiles_to]
    windows                 = windows[tiles_from:tiles_to]
    windows_extended        = windows_extended[tiles_from:tiles_to]


else:
    process_subset = True
    windows_ids_progressive = windows_id_progressive_subset
    window_ids              = window_ids_subset
    windows                 = windows_subset
    windows_extended        = windows_extended_subset
    pass



# =============================================================================
# windows_id_progressive = windows_ids_progressive[0]
# window_id=window_ids[0]
# window_extended=windows_extended[0]
# window=windows[0]
# =============================================================================



#master function that runs in parallel
def process_tiles(windows_id_progressive, window_id, window_extended, window ):
    
    terreno, tf_extended           = extract_window_raster(dtm_path, window_extended)
    _, tf                          = extract_window_raster(dtm_path, window)
    
    
    #mask covering the overlap area between the windows
    buffer_mask = np.ones_like(terreno,dtype='float32')
    # removes the overlaps between the tiles
    buffer_mask[:buffer_size_window, :]  = np.nan  # First two rows
    buffer_mask[-buffer_size_window:, :] = np.nan # Last two rows  
    buffer_mask[:, :buffer_size_window]  = np.nan  # First two columns
    buffer_mask[:, -buffer_size_window:] = np.nan # Last two columns
    
    
    if np.all(np.isnan(terreno*buffer_mask)):
        print(f'All no-data in tile {window_id}\n')
        return
    

    try:
        
        # Calculate the bounding box of the processing window 
        min_x, max_y = xy(tf_extended, 0, 0, offset='ul')
        max_x, min_y = xy(tf_extended, window_extended.height, window_extended.width, offset='lr')
        bbox = (min_x, min_y, max_x, max_y)
        
        
        
        
# =============================================================================
#         #AUXILIARY DEM######
#         
#         #transform the bbox if ancillary data have different crs
#         bbox_transformed = transform_bounds(dem_proj4, dtm_aux_proj4, min_x, min_y, max_x, max_y)
#                
#         #get transforms from bbox
#         window_dtm_aux = from_bounds(*bbox_transformed, transform=dtm_aux_gT)
#          
#         
#         #extract windows of data from auxiliary dem
#         terreno_dtm_aux, tf_dtm_aux = extract_window_raster(dtm_aux_path, window_dtm_aux)
#     
#         if terreno_dtm_aux.size == 0:
#             terreno_dtm_aux_resampled = np.full(terreno.shape, np.nan)
#             print('dtm_aux outside domain')
#         else:
#             terreno_dtm_aux_resampled_raster = bring_raster_to_same_grid(source_array = terreno_dtm_aux,
#                                                                          source_crs   = dtm_aux_proj4,
#                                                                          source_tf    = tf_dtm_aux, 
#                                                                          target_crs   = dem_proj4,
#                                                                          target_tf    = tf_extended,
#                                                                          target_width = terreno.shape[1],
#                                                                          target_height= terreno.shape[0] )
#         
#             terreno_dtm_aux_resampled = terreno_dtm_aux_resampled_raster[0]
# =============================================================================
        
        
        
        #GSWE#########
        bbox_transformed = transform_bounds(dem_proj4, gswe_proj4, min_x, min_y, max_x, max_y)
        
        #get transforms from bbox
        gswe_window = from_bounds(*bbox_transformed, transform=gswe_gT)
         
        #extract windows of data 
        data_gswe_window, tf_gswe_window = extract_window_raster(glob_s_w_explorer_path, gswe_window)
    
        if data_gswe_window.size == 0:
            data_gswe_window_resampled = np.full(terreno.shape, 0, dtype=np.uint8)
            print('GSWE outside domain')
        else:
            data_gswe_window_resampled_raster = bring_raster_to_same_grid(source_array = data_gswe_window,
                                                                          source_crs   = gswe_proj4,
                                                                          source_tf    = tf_gswe_window, 
                                                                          target_crs   = dem_proj4,
                                                                          target_tf    = tf_extended,
                                                                          target_width = terreno.shape[1],
                                                                          target_height= terreno.shape[0] )
        
            data_gswe_window_resampled = data_gswe_window_resampled_raster[0]
        
        

 
        
        #mask escluding nans
        not_mask= np.bitwise_not(np.isnan(terreno))
        
        h=calculate_hand(terreno, 
                         #terreno_dtm_aux_resampled,
                         data_gswe_window_resampled,
                         tf_extended, 
                         pyproj.Proj(dem.crs.to_string()), 
                         catch_mask=not_mask, 
                         acc_thresh=acc_thresh, 
                         routing = routing,
                         minimum_area_threshold_water_bodies  = minimum_area_threshold_water_bodies,
                         geographic_coordinates = geographic_coordinates,
                         ).astype('float32')
        
        
        
        h = h*buffer_mask




        # Use the Window to slice the NumPy array and exclude the buffer area around
        row_start = window.row_off - window_extended.row_off
        row_end   = row_start + window.height
        col_start = window.col_off - window_extended.col_off
        col_end   = col_start + window.width

        # Extract the subarray
        h_extracted = h[row_start:row_end, col_start:col_end]

        
        save_tiff(h_extracted, dem_proj4, tf, output_folder,  file_name= f'{window_id}', compression = 'deflate', dtype='float32')
        
    
        print(f'Completed tile {window_id}. Progress: {int(windows_id_progressive)+1}/{number_windows}\n')
        

    except:
        print(f'Error in tile: {window_id}, progressive id: {windows_id_progressive}')




#RUN IN PARALLEL
def main():

    # Prepare arguments for each process 
    if process_subset == True:
        args_list = list(zip(windows_id_progressive_subset, window_ids_subset, windows_extended_subset, windows_subset))
    else:
        args_list = list(zip(windows_ids_progressive, window_ids, windows_extended, windows))
    

    # Set up multiprocessing pool with limited processes
    with multiprocessing.Pool(processes=15) as pool:  
        pool.starmap(process_tiles, args_list)

if __name__ == '__main__':
    main() 






