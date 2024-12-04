import streamlit as st
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
import pandas as pd
import folium
from streamlit_folium import folium_static

def load_dem_files():
    """
    Load predefined DEM files from GitHub or local sources'''
    
    """
    vector_path = 'https://raw.githubusercontent.com/rk2026/data/main/Bagmati_ward.gpkg'
    dem_files = {
        'DHADING_Thakre.tif': 'https://raw.githubusercontent.com/rk2026/data/main/DHADING_Thakre.tif',
        'DHADING_Gajuri.tif': 'https://raw.githubusercontent.com/rk2026/data/main/DHADING_Gajuri.tif',
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
        
        # Convert results to DataFrame
        results_df = pd.DataFrame(results)
        return results_df

def create_interactive_map(dem_results, vector_path):
    """
    Create an interactive map with administrative boundaries and points
    """
    # Read vector layer
    vector_gdf = gpd.read_file(vector_path)
    
    # Create a map centered on the mean coordinates of the vector layer
    center_lat = vector_gdf.geometry.centroid.y.mean()
    center_lon = vector_gdf.geometry.centroid.x.mean()
    
    m = folium.Map(location=[center_lat, center_lon], zoom_start=10)
    
    # Add base layers (using folium tile layers)
    folium.TileLayer('OpenStreetMap').add_to(m)
    folium.TileLayer('Stamen Terrain').add_to(m)
    
    # Add vector layer with popup
    def style_function(feature):
        return {
            'fillColor': 'blue',
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.1
        }
    
    def popup_function(feature):
        props = feature['properties']
        return folium.Popup(f"""
            District: {props['DISTRICT']}
            Ward: {props['NEW_WARD_N']}
            Type: {props['Type_GN']}
        """)
    
    folium.GeoJson(
        vector_gdf, 
        style_function=style_function,
        popup=popup_function
    ).add_to(m)
    
    # Add minimum height points
    for idx, row in dem_results.iterrows():
        # Minimum height point
        folium.CircleMarker(
            location=[row['min_lat'], row['min_lon']],
            radius=5,
            popup=f"Minimum Elevation: {row['min_elevation']:.2f}",
            color='green',
            fill=True,
            fillColor='green'
        ).add_to(m)
        
        # Maximum height point
        folium.CircleMarker(
            location=[row['max_lat'], row['max_lon']],
            radius=5,
            popup=f"Maximum Elevation: {row['max_elevation']:.2f}",
            color='red',
            fill=True,
            fillColor='red'
        ).add_to(m)
    
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
            try:
                # Calculate zonal statistics
                dem_results = calculate_zonal_statistics(
                    selected_dem, 
                    vector_path
                )
                
                # Display results table
                st.dataframe(dem_results)
                
                # Create interactive map
                map_obj = create_interactive_map(dem_results, vector_path)
                folium_static(map_obj)
                
            except Exception as e:
                st.error(f"An error occurred: {str(e)}")

if __name__ == '__main__':
    main()

# requirements.txt:
# streamlit
# geopandas
# rasterio
# numpy
# pandas
# folium
# streamlit-folium
