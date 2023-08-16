# extract_netcdf
This repository stores an optimized framework to extract the average time-series (weighted or not) encompassing different shapefiles from the E-Obs datasets 

In the script example, we are dealing with the rainfall time-series file from E-obs called: rr_ens_mean_0.25deg_reg_v27.0e.nc as the netcdf file. But the code may be adapted to other netcdfs. As shapefiles we are using the individual catchment boundaries. The task is beinjg obtimized in a way to provide a fast processing since the idea is to deal with more than 10000 different catchments and a non-parallel processing would take too much time. 

References: Cornes, R., G. van der Schrier, E.J.M. van den Besselaar, and P.D. Jones. 2018: An Ensemble Version of the E-OBS Temperature and Precipitation Datasets, J. Geophys. Res. Atmos., 123. doi:10.1029/2017JD028200
