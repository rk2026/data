import streamlit as st
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
import pandas as pd
import leafmap.folium as leafmap

def load_dem_files():
    """
    Load predefined DEM files from GitHub or local sources
    """
    dem_files = {
        'DHADING_Thakre.tif': 'https://github.com/your-repo/DHADING_Thakre.tif',
        # Add more predefined DEM files here
    }
    return dem_files

def calculate_zonal_statistics(dem_path, vector_path):
    """
    Calculate zonal statistics for each polygon in the vector layer
    """
    # Read the raster and vector data
    with rasterio.open(dem_path) as dem_src:
        dem_data = dem_src.read(1)
        dem_transform = dem_src.transform
        
        # Read vector layer
        gdf = gpd.read_file(vector_path)
        
        # Initialize lists to store results
        results = []
        
        # Iterate through each polygon
        for index, row in gdf.iterrows():
            # Create a mask for the current polygon
            out_image, out_transform = mask(
                dem_src, 
                [row.geometry.__geo_interface__], 
                crop=True, 
                nodata=np.nan
            )
            
            # Calculate statistics
            masked_data = out_image[0]
            valid_data = masked_data[~np.isnan(masked_data)]
            
            if len(valid_data) > 0:
                min_val = np.min(valid_data)
                max_val = np.max(valid_data)
                
                # Find pixel coordinates for min and max
                min_pixel = np.where(masked_data == min_val)
                max_pixel = np.where(masked_data == max_val)
                
                # Convert pixel coordinates to geographic coordinates
                min_lon, min_lat = rasterio.transform.xy(
                    out_transform, min_pixel[0][0], min_pixel[1][0]
                )
                max_lon, max_lat = rasterio.transform.xy(
                    out_transform, max_pixel[0][0], max_pixel[1][0]
                )
                
                results.append({
                    'DISTRICT': row['DISTRICT'],
                    'GaPa_NaPa': row['GaPa_NaPa'],
                    'NEW_WARD_N': row['NEW_WARD_N'],
                    'Type_GN': row['Type_GN'],
                    'min_elevation': min_val,
                    'max_elevation': max_val,
                    'min_lon': min_lon,
                    'min_lat': min_lat,
                    'max_lon': max_lon,
                    'max_lat': max_lat
                })
        
        # Convert results to GeoDataFrame
        results_df = gpd.GeoDataFrame(results)
        return results_df

def create_interactive_map(dem_results, vector_path):
    """
    Create an interactive map with DEM and administrative boundaries
    """
    m = leafmap.Map(center=[27.7172, 85.3240], zoom=10)
    
    # Add base layers
    m.add_basemap('ESRI.WorldImagery')
    m.add_basemap('OpenStreetMap')
    
    # Add DEM layer (if possible)
    # m.add_raster(dem_path, colormap='terrain', layer_name='DEM')
    
    # Add minimum and maximum height points
    min_points = gpd.GeoDataFrame(
        dem_results, 
        geometry=gpd.points_from_xy(
            dem_results['min_lon'], 
            dem_results['min_lat']
        )
    )
    max_points = gpd.GeoDataFrame(
        dem_results, 
        geometry=gpd.points_from_xy(
            dem_results['max_lon'], 
            dem_results['max_lat']
        )
    )
    
    # Add vector layer
    vector_gdf = gpd.read_file(vector_path)
    
    m.add_gdf(vector_gdf, layer_name='Administrative Boundaries')
    m.add_gdf(min_points, layer_name='Minimum Height Points')
    m.add_gdf(max_points, layer_name='Maximum Height Points')
    
    return m

def main():
    st.title('DEM Zonal Statistics Analysis')
    
    # Sidebar for file selection
    st.sidebar.header('Data Selection')
    
    # DEM File Selection
    dem_files = load_dem_files()
    selected_dem = st.sidebar.selectbox(
        'Select DEM File', 
        list(dem_files.keys())
    )
    
    # Vector Layer Selection (hardcoded for now)
    vector_path = 'Bagmati_ward.gpkg'
    
    # Process data
    if st.sidebar.button('Analyze DEM'):
        with st.spinner('Processing DEM and calculating zonal statistics...'):
            # Calculate zonal statistics
            dem_results = calculate_zonal_statistics(
                dem_files[selected_dem], 
                vector_path
            )
            
            # Display results table
            st.dataframe(dem_results)
            
            # Create interactive map
            map_obj = create_interactive_map(dem_results, vector_path)
            map_obj.to_streamlit(height=600)

if __name__ == '__main__':
    main()

# Requirements (requirements.txt):
# streamlit
# geopandas
# rasterio
# numpy
# pandas
# leafmap
