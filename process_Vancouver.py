import geopandas as gpd
from netCDF4 import Dataset
from regrid import calculate_population_in_grid, get_wrf_grid_info, saveNC, preprocess_gdf
from warnings import filterwarnings
filterwarnings('ignore')

def population(wrf_file, shapefile, output, landshapefile, cols, vanIf=True):
    print(wrf_file)
    # Load WRF grid info
    wrf_info = get_wrf_grid_info(wrf_file)

    # Load shapefile
    gdf = gpd.read_file(shapefile)[cols + ['geometry']]
    if not vanIf: gdf = preprocess_gdf(gdf, cols)
    print(gdf)
    # Transform shapefile to the WRF projection
    gdf = gdf.to_crs(wrf_info['crs'])
    Lgdf = read_Lgdf(landshapefile, wrf_info)

    # Regriding data
    population_grid = calculate_population_in_grid(gdf, wrf_info, Lgdf, cols)
    # Load the WRF data,
    ds = Dataset(wrf_file)
    saveNC(output, wrf_info, population_grid, cols, src=ds)
    ds.close()

    
def read_Lgdf(landshapefile, wrf_info):
    return gpd.read_file(landshapefile).to_crs(wrf_info['crs']).unary_union

def main():
    landshapefile = '/mnt/c/Users/foroo/Documents/Old_PC/PhD/Term6/Population/Output/2018.shp'
    
    wrfinp, output = '../../../../Term6/Population2/Victoria_domain.nc', './Victoria.nc'
    inp = './Victoria_Impervious_Tree/Victoria_imp_tree.shp'
    cols = ['PctFor86', 'PctImp86', 'PctFor05', 'PctImp05', 'PctImp11', 'PctFor11', 'PctImp19Or', 'PctFor19Or']
    population(wrfinp, inp, output, landshapefile, cols, vanIf=False)
    exit()
    wrfinp, output = '../../../../Term6/Population2/Vancouver_domain.nc', './Vancouver.nc'
    inp = './Vancouver_Impervious_Tree/2014_and_2020_Tree_Canopy_Cover_and_Impervious_Surface.shp'
    cols = ['CC2020_m2', 'CC2014_m2', 'IMP2020_m2', 'IMP2014_m2']
    population(wrfinp, inp, output, landshapefile, cols, vanIf=True)


if __name__ == '__main__':
    main()
