import netCDF4 as nc
import geopandas as gpd
import numpy as np
import tqdm
import os
from concurrent.futures import ThreadPoolExecutor
import time
from utils import process_catchment

# Path where the ncdf files are stored:
path = r"../datasets//"
path_out = r"../timeseries//"

# ncdf file's name:
nc_file = "rr_ens_mean_0.25deg_reg_v27.0e.nc"

# Set the number of workers for parallel processing
num_workers = 4

# Read NetCDF data
with nc.Dataset(path+nc_file, mode='r', format='NETCDF4') as nc_dataset:
    latitude = nc_dataset["latitude"][:]
    longitude = nc_dataset["longitude"][:]
    values = nc_dataset["rr"][:]

# Chunk size for reading NetCDF data
chunk_size = 100  # Adjust this value based on your available memory

# Calculate pixel extents using vectorized operations
lon_idx, lat_idx = np.meshgrid(range(len(longitude)), range(len(latitude)))
pixel_extents = np.stack(
    (longitude[lon_idx], latitude[lat_idx], longitude[lon_idx] + 0.25, latitude[lat_idx] + 0.25),
    axis=-1
)

# Read shapefile data
shapefile_all = gpd.read_file(path+"catchments_EU2caravan_v15.shp")

# Set the CRS of the shapefile's geometry to EPSG:4326 (WGS 84)
shapefile_all["geometry"] = shapefile_all["geometry"].to_crs(epsg=4326)

# Extract catchment names
catchmentnames = shapefile_all.Code.tolist()

#%%
# Process catchment polygons in parallel
start = time.time()
with ThreadPoolExecutor(max_workers=num_workers) as executor:
    futures = [executor.submit(process_catchment, catchmentname, shapefile_all, values, latitude, longitude, path_out)
               for catchmentname in tqdm.tqdm(catchmentnames[35:37])]
    
    # Wait for all futures to complete
    for future in futures:
        future.result()

end = time.time()
print("Total time:", end - start)









