import netCDF4 as nc
import geopandas as gpd
import numpy as np
from shapely.geometry import Point
import tqdm
import os
from concurrent.futures import ThreadPoolExecutor
import time

# Set your paths
path = r"C:\Users\nascimth\Documents\Thiago\Eawag\Python\Scripts\Data\EUobs\extract_netcdf"
os.chdir(path)
os.getcwd()

from utils import get_pixel_indices_and_coords, process_catchment

# Path where the ncdf files are stored:
path = r"C:\Users\nascimth\Documents\Thiago\Eawag\Python\Data\Meteorological\EObs\datasets\\"
path_out = r"C:\Users\nascimth\Documents\Thiago\Eawag\Python\Data\Meteorological\EObs\timeseries\\"

# ncdf file's name:
nc_file = "rr_ens_mean_0.25deg_reg_v27.0e.nc"

# Set the number of workers for parallel processing
num_workers = 6

# Read NetCDF data
with nc.Dataset(path + nc_file, mode='r', format='NETCDF4') as nc_dataset:
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
shapefile_path = r"C:\Users\nascimth\Documents\Thiago\Eawag\Python\GIS\Catchments\EU\GIS\catchments_EU2caravan_v15.shp"
shapefile_all = gpd.read_file(shapefile_path)

# Set the CRS of the shapefile's geometry to EPSG:4326 (WGS 84)
shapefile_all["geometry"] = shapefile_all["geometry"].to_crs(epsg=4326)

# Extract catchment names
catchmentnames = shapefile_all.Code.tolist()

# Process catchment polygons in parallel
start = time.time()
with ThreadPoolExecutor(max_workers=num_workers) as executor:
    futures = [executor.submit(process_catchment, catchmentname, shapefile_all, path_out)
               for catchmentname in tqdm.tqdm(catchmentnames[0:10])]
    
    # Wait for all futures to complete
    for future in futures:
        future.result()

end = time.time()
print("Total time:", end - start)


















#%% Now we retrieve the data for each catchment:

# Create a function to process a single catchment polygon
def process_catchment(catchmentname, shapefile_all, values, latitude, longitude, path_out):
    shapefile = shapefile_all[shapefile_all.Code == catchmentname]
    polygon = shapefile.geometry.unary_union
    pixel_indices, pixel_coords = get_pixel_indices_and_coords(latitude, longitude, polygon)

    pixel_points = [Point(lon, lat) for lat, lon in zip(latitude[pixel_indices[:, 0]], longitude[pixel_indices[:, 1]])]
    pixel_polygon = gpd.GeoSeries(pixel_points).buffer(0.125)

    intersection_areas = np.array(pixel_polygon.intersection(shapefile.geometry.unary_union).area)
    weights_time_series = intersection_areas / intersection_areas.sum()

    num_time_steps = values.shape[0]
    num_pixels = len(pixel_indices)
    weighted_time_series = np.zeros((num_time_steps, num_pixels))

    for start_idx in range(0, num_time_steps, chunk_size):
        end_idx = min(start_idx + chunk_size, num_time_steps)
        chunk = values[start_idx:end_idx]

        chunk_weighted = chunk[:, pixel_indices[:, 0], pixel_indices[:, 1]]
        weighted_time_series[start_idx:end_idx] = chunk_weighted

    weighted_sum = np.sum(weighted_time_series * weights_time_series, axis=1, keepdims=True)
    pixels_array = np.hstack((pixel_coords, intersection_areas.reshape((len(intersection_areas), 1)),
                              weights_time_series.reshape((len(weights_time_series), 1))))

    np.savetxt(path_out + "pixels_" + catchmentname + ".csv", pixels_array, delimiter=',')
    np.savetxt(path_out + "rain_" + catchmentname + ".csv", weighted_sum, delimiter=',')


def process_catchments(catchmentnames, shapefile_all, values, latitude, longitude, path_out):
    with ProcessPoolExecutor() as executor:  # Use ProcessPoolExecutor for parallelism
        futures = [executor.submit(process_catchment, catchmentname, shapefile_all, values, latitude, longitude, path_out)
                   for catchmentname in catchmentnames]

        # Wait for all futures to complete
        for future in futures:
            future.result()

start = time.time()

# Process catchment polygons in parallel using processes
process_catchments(catchmentnames[0:2], shapefile_all, values, latitude, longitude, path_out)

end = time.time()
print("Total time:", end - start)


start = time.time()
# Process catchment polygons in parallel
with ThreadPoolExecutor(max_workers=4) as executor:  # Adjust the number of workers as needed
    futures = [executor.submit(process_catchment, catchmentname, shapefile_all, values, latitude, longitude, path_out)
               for catchmentname in tqdm.tqdm(catchmentnames[0:20])]

    # Wait for all futures to complete
    for future in futures:
        future.result()

end = time.time()
print(end - start)


# Get the number of available CPU cores
num_cores = multiprocessing.cpu_count()





