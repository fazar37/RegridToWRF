import pyproj
import xarray as xr
import rasterio
from shapely.geometry import box
from netCDF4 import Dataset
import os
from subprocess import check_call
import numpy as np
from shapely.strtree import STRtree
# def calculate_population_in_grid(gdf, wrf_info, Lgdf, cols):
#     # Initialize an empty array for population data in WRF grid
#     o1 = np.zeros((wrf_info['ny'], wrf_info['nx']))
#     population_grid = {col:o1.copy() for col in cols}
    
#     # Create a grid of polygons representing each WRF grid cell
#     grid_polygons = []
#     for i in range(wrf_info['ny']):
#         for j in range(wrf_info['nx']):
#             x, y = wrf_info['transform'] * (j, i)
#             cell = box(x, y, x + wrf_info['transform'].a, y + wrf_info['transform'].e)
#             grid_polygons.append((i, j, cell))
    
#     # Loop over each polygon in the shapefile
#     for _, row in gdf.iterrows():
#         print(_)
#         if _ == 200: break
#         geom = row.geometry
#         # pop = row[population_column]
        
#         # Skip if the population polygon is empty or invalid
#         if geom.is_empty or not geom.is_valid: continue
        
#         intersection_with_land = geom.intersection(Lgdf)
#         if intersection_with_land.is_empty or not intersection_with_land.is_valid: continue
#         intersection_with_land_area = intersection_with_land.area
#         # Calculate the intersection area for each grid cell
#         for i, j, cell in grid_polygons:
#             intersection = geom.intersection(cell)
#             land_intersection = intersection.intersection(intersection_with_land)
#             if not land_intersection.is_empty and land_intersection.is_valid:
#                 # Calculate the area weight based on the land portion
#                 area_weight = land_intersection.area / intersection_with_land_area
                
#                 # Distribute the population based on the intersection area
#                 for ii, col in enumerate(cols):
#                     population_grid[col][i, j] += row[col] * area_weight
    
#     for col in cols:
#         for i, j, cell in grid_polygons:
#             population_grid[col][i, j] = population_grid[col][i, j] / cell.area

#     return population_grid

def calculate_population_in_grid(gdf, wrf_info, Lgdf, cols):
    # Initialize output arrays
    population_grid = {col: np.zeros((wrf_info['ny'], wrf_info['nx'])) for col in cols}

    # Create grid polygons and index
    grid_cells = []
    ij_index = {}

    for i in range(wrf_info['ny']):
        for j in range(wrf_info['nx']):
            x, y = wrf_info['transform'] * (j, i)
            cell = box(x, y, x + wrf_info['transform'].a, y + wrf_info['transform'].e)
            grid_cells.append(cell)
            ij_index[cell] = (i, j)

    # Build spatial index for the grid
    tree = STRtree(grid_cells)

    for idx, row in gdf.iterrows():
        # print(idx)
        # if idx == 200: break
        geom = row.geometry
        if geom.is_empty or not geom.is_valid: continue

        # Clip the feature to land only
        geom_land = geom.intersection(Lgdf)
        if geom_land.is_empty or not geom_land.is_valid: continue

        total_land_area = geom_land.area
        if total_land_area == 0: continue

        # Query spatial index for overlapping grid cells
        overlapping_cells = tree.query(geom_land)

        for cell_id in overlapping_cells:
            cell = tree.geometries.take(cell_id)
            intersection = geom_land.intersection(cell)
            if intersection.is_empty or not intersection.is_valid: continue

            area_weight = intersection.area / total_land_area
            i, j = ij_index[cell]
            for col in cols:
                population_grid[col][i, j] += row[col] * area_weight

    # Normalize by cell area
    cell_area = abs(wrf_info['transform'].a * wrf_info['transform'].e)
    for col in cols:
        population_grid[col] /= cell_area

    return population_grid

def get_wrf_grid_info(wrfinp):
    # Initialize WGS84 projection
    wgs_proj = pyproj.Proj(proj='latlong', datum='WGS84')
    
    # Read gridded WRF data
    dst_data = xr.open_dataset(wrfinp)
    
    # Lambert Conformal Conic
    wrf_proj = pyproj.Proj(
        proj='lcc',
        units='m',
        a=6370000,
        b=6370000,
        lat_1=dst_data.TRUELAT1,
        lat_2=dst_data.TRUELAT2,
        lat_0=dst_data.MOAD_CEN_LAT,
        lon_0=dst_data.STAND_LON
    )
    
    # Make transform
    transformer_wrf = pyproj.Transformer.from_proj(wgs_proj, wrf_proj)
    e, n = transformer_wrf.transform(dst_data.CEN_LON, dst_data.CEN_LAT)
    
    # Grid parameters
    dx, dy = dst_data.DX, dst_data.DY
    nx, ny = dst_data.dims['west_east'], dst_data.dims['south_north']
    
    # Down left corner of the domain
    x0 = -(nx - 1) / 2.0 * dx + e
    y0 = -(ny - 1) / 2.0 * dy + n
    
    wrf_transform = rasterio.transform.Affine.translation(x0 - dx / 2, y0 - dy / 2) * rasterio.transform.Affine.scale(dx, dy)
    wrf_grid_info = {'crs': wrf_proj.to_proj4(), 'transform': wrf_transform, 'nx': nx, 'ny': ny}
    
    return wrf_grid_info

def saveNC(output_netcdf:str, wrf_info, values, cols, src):
    if os.path.exists(output_netcdf): check_call(f'rm -f {output_netcdf}', shell=True)
    with Dataset(output_netcdf, "w", format="NETCDF4") as nc_file:
        # Create dimensions
        nc_file.createDimension("lat", wrf_info['ny'])
        nc_file.createDimension("lon", wrf_info['nx'])

        # Set the global attributes
        for attr_name in src.ncattrs():
            nc_file.setncattr(attr_name, src.getncattr(attr_name))

        # Create lat and lon variables
        lat_var = nc_file.createVariable("XLAT", "f4", ("lat", 'lon',))
        lon_var = nc_file.createVariable("XLONG", "f4", ('lat', "lon",))
        lat_var.units = "degrees_north"
        lon_var.units = "degrees_east"
        lat_var[:] = src['XLAT'][:,:]
        lon_var[:] = src['XLONG'][:,:]

        # Create other variables
        for col in cols:
            other_var = nc_file.createVariable(col, "f4", ( "lat", "lon"), zlib=True)
            other_var.units = "%"        
            other_var[:, :] = values[col]

def preprocess_gdf(gdf, cols):
    for col in cols:
        gdf[col] = gdf[col] * gdf.geometry.area / 100
    return gdf
