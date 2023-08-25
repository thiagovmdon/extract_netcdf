# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 15:21:01 2023

@author: nascimth
"""

import netCDF4 as nc
import geopandas as gpd
import numpy as np
import tqdm
import os
from concurrent.futures import ThreadPoolExecutor
import time

# Set your paths
path = r"C:\Users\nascimth\Documents\Thiago\Eawag\Python\Scripts\Data\EUobs\extract_netcdf\extract_netcdf"
os.chdir(path)
os.getcwd()

from utils import process_catchment

# Paths used in this script:
path = r"C:\Users\nascimth\Documents\Thiago\Eawag\Python\Data\Meteorological\EObs\datasets\\"
path_out = r"C:\Users\nascimth\Documents\Thiago\Eawag\Python\Data\Meteorological\EObs\timeseries\\"
shapefile_path = r"C:\Users\nascimth\Documents\Thiago\Eawag\Python\GIS\Catchments\EU\GIS\catchments_EU2caravan_v16.shp"

# ncdf file's name:
# Precipitation:
nc_file = "rr_ens_mean_0.25deg_reg_v27.0e.nc"
variable_name = "rr"

# Potential evapotranspiration:
#nc_file = "hargreaves2.nc"
#variable_name = "Hargreaves"

# Mean temperature:
#nc_file = "tg_ens_mean_0.25deg_reg_v27.0e.nc"
#variable_name = "tg"


"""
# For Europe and without Iceland:
lat_min = 34.75
lat_max = 70.6
lon_min = -10.5
lon_max = 45.0


# For Iceland:
lat_min = 63.0
lat_max = 67.0
lon_min = -24.5
lon_max = -10.5
"""

# Read NetCDF data
with nc.Dataset(path + nc_file, mode='r', format='NETCDF4') as nc_dataset:
    latitude = nc_dataset["latitude"][:]
    longitude = nc_dataset["longitude"][:]
    values = nc_dataset[variable_name][:]

# Here we do this simple print only for checking our variables' names:
print("Variables in NetCDF file:", nc_dataset.variables.keys())

# Set the number of workers for parallel processing
num_workers = 5

# Chunk size for reading NetCDF data
chunk_size = 100  # Adjust this value based on your available memory

# Calculate pixel extents using vectorized operations
lon_idx, lat_idx = np.meshgrid(range(len(longitude)), range(len(latitude)))
pixel_extents = np.stack(
    (longitude[lon_idx], latitude[lat_idx], longitude[lon_idx] + 0.25, latitude[lat_idx] + 0.25),
    axis=-1
)

# Read shapefile data
shapefile_all = gpd.read_file(shapefile_path)

# Set the CRS of the shapefile's geometry to EPSG:4326 (WGS 84)
shapefile_all["geometry"] = shapefile_all["geometry"].to_crs(epsg=4326)

# Extract catchment names
catchmentnames = shapefile_all.Code.tolist()

#%%
# Process catchment polygons in parallel
start = time.time()
with ThreadPoolExecutor(max_workers=num_workers) as executor:
    futures = [executor.submit(process_catchment, catchmentname, shapefile_all, values, latitude, longitude, path_out, variable_name = "rr")
               for catchmentname in tqdm.tqdm(catchmentnames)]
    
    # Wait for all futures to complete
    for future in futures:
        future.result()

end = time.time()
print("Total time:", end - start)


