import streamlit as st
import geopandas as gpd
import rasterio
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from rasterio.mask import mask
import requests
from io import BytesIO
import shapely
import folium
from streamlit_folium import folium_static

# List of available DEM files
DEM_FILES = [
    'Dhading_Netrawati.tif',
    'Dhading_Khaniyabash.tif', 
    'Dhading_Jwalamukhi.tif', 
    'Dhading_galchhi.tif', 
    'Dhading_Gajuri.tif', 
    'Dhading_dhunibesi.tif', 
    'Dhading_Benighat_Rorang.tif'
]

def fetch_github_file(url):
    """Fetch file from GitHub URL"""
    try:
        response = requests.get(url)
        response.raise_for_status()
        return BytesIO(response.content)
    except Exception as e:
        st.error(f"Error fetching file: {e}")
        return None

def create_2d_map(zonal_results):
    """
    Create a 2D map using Folium with elevation points
    
    Args:
        zonal_results (GeoDataFrame): Zonal statistics results
    
    Returns:
        folium.Map: Interactive map with elevation points
    """
    # Calculate center point from the data
    center_lat = zonal_results['min_lat'].mean()
    center_lon = zonal_results['min_lon'].mean()
    
    # Create base map
    m = folium.Map(location=[center_lat, center_lon], zoom_start=10)
    
    # Add markers for min and max elevation points
    for _, row in zonal_results.iterrows():
        # Minimum elevation marker
        folium.CircleMarker(
            location=[row['min_lat'], row['min_lon']],
            radius=5,
            popup=f"Minimum Elevation: {row.get('min_elevation', 'N/A')}m<br>"
                  f"Ward: {row.get('NEW_WARD_N', 'N/A')}",
            color='blue',
            fill=True,
            fillColor='blue'
        ).add_to(m)
        
        # Maximum elevation marker
        folium.CircleMarker(
            location=[row['max_lat'], row['max_lon']],
            radius=5,
            popup=f"Maximum Elevation: {row.get('max_elevation', 'N/A')}m<br>"
                  f"Ward: {row.get('NEW_WARD_N', 'N/A')}",
            color='red',
            fill=True,
            fillColor='red'
        ).add_to(m)
    
    return m

# [Rest of the previous code remains the same, just add the create_2d_map function]

def main():
    st.title("DEM Analysis and Visualization")
    
    # Dropdown for selecting DEM file
    selected_dem = st.selectbox("Select DEM File", DEM_FILES)
    
    # Construct GitHub URL for the selected DEM
    base_url = "https://github.com/rk2026/data/raw/main/"
    dem_url = base_url + selected_dem
    
    # Fixed Vector Layer URL (assuming it remains the same)
    vector_url = base_url + "Bagmati_ward.gpkg"
    
    if st.button("Process Data"):
        with st.spinner(f"Processing {selected_dem}..."):
            # Fetch files
            dem_file = fetch_github_file(dem_url)
            vector_file = fetch_github_file(vector_url)
            
            if dem_file and vector_file:
                # Process zonal statistics
                zonal_results, intersecting_gdf = process_dem_zonal_stats(dem_file, vector_file)
                
                # Display results table
                st.subheader(f"Zonal Statistics Results for {selected_dem}")
                # Select columns to display
                display_columns = [
                    'DISTRICT', 'GaPa_NaPa', 'Type_GN', 'NEW_WARD_N', 
                    'min_elevation', 'max_elevation', 
                    'min_lon', 'min_lat', 'max_lon', 'max_lat'
                ]
                
                # Filter columns that exist in the results
                available_columns = [col for col in display_columns if col in zonal_results.columns]
                
                st.dataframe(zonal_results[available_columns])
                
                # Create 3D visualization
                fig = create_3d_visualization(dem_file, intersecting_gdf, zonal_results)
                st.plotly_chart(fig, use_container_width=True)
                
                # Create and display 2D map
                st.subheader("2D Elevation Map")
                m = create_2d_map(zonal_results)
                folium_static(m)

if __name__ == "__main__":
    main()
