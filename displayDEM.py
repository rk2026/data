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
    'DHADING_Netrawati.tif',
    'DHADING_Khaniyabash.tif', 
    'DHADING_Jwalamukhi.tif', 
    'DHADING_Galchi.tif', 
    'DHADING_Gajuri.tif', 
    'DHADING_Benighat_Rorang.tif', 
    'DHADING_Nilkantha.tif',
    'DHADING_Gangajamuna.tif',
    'DHADING_Dhunibesi.tif',
    'DHADING_Siddhalek.tif',
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

def process_dem_zonal_stats(dem_path, vector_path):
    """
    Perform zonal statistics on DEM raster within vector polygons
    
    Returns:
        GeoDataFrame with zonal statistics
    """
    # Read the vector layer
    gdf = gpd.read_file(vector_path)
    
    # Read the DEM raster
    with rasterio.open(dem_path) as dem:
        # Get raster bounds
        raster_bounds = shapely.geometry.box(*dem.bounds)
        
        # Reproject vector to match raster CRS if needed
        if gdf.crs != dem.crs:
            gdf = gdf.to_crs(dem.crs)
        
        # Filter polygons that intersect with raster bounds
        intersecting_gdf = gdf[gdf.intersects(raster_bounds)]
        
        # Initialize results list
        zonal_results = []
        
        # Process each polygon
        for idx, row in intersecting_gdf.iterrows():
            try:
                # Create a mask for the current polygon
                mask_geometry = [row.geometry]
                
                # Clip raster to polygon
                out_image, out_transform = mask(dem, mask_geometry, crop=True)
                
                # Flatten and remove nodata values
                valid_pixels = out_image[out_image != dem.nodata]
                
                if len(valid_pixels) > 0:
                    # Calculate zonal statistics
                    min_val = float(np.min(valid_pixels))
                    max_val = float(np.max(valid_pixels))
                    
                    # Only include if at least one value is > 0
                    if max_val > 0:
                        # Find pixel locations for min and max
                        min_pixel_idx = np.unravel_index(np.argmin(out_image), out_image.shape)
                        max_pixel_idx = np.unravel_index(np.argmax(out_image), out_image.shape)
                        
                        # Convert pixel indices to geospatial coordinates
                        min_lon, min_lat = rasterio.transform.xy(out_transform, min_pixel_idx[0], min_pixel_idx[1])
                        max_lon, max_lat = rasterio.transform.xy(out_transform, max_pixel_idx[0], max_pixel_idx[1])
                        
                        # Combine zonal results with original polygon attributes
                        result_row = {
                            'min_elevation': min_val,
                            'max_elevation': max_val,
                            'min_lon': min_lon,
                            'min_lat': min_lat,
                            'max_lon': max_lon,
                            'max_lat': max_lat,
                            'geometry': row.geometry
                        }
                        
                        # Add additional attributes from original vector data
                        additional_attrs = ['DISTRICT', 'GaPa_NaPa', 'Type_GN', 'NEW_WARD_N']
                        for attr in additional_attrs:
                            if attr in row.index:
                                result_row[attr] = row[attr]
                        
                        zonal_results.append(result_row)
            
            except Exception as e:
                st.warning(f"Could not process polygon: {e}")
        
        # Create GeoDataFrame from results
        results_gdf = gpd.GeoDataFrame(zonal_results, crs=gdf.crs)
        return results_gdf, intersecting_gdf

def create_3d_visualization(dem_path, intersecting_gdf, zonal_results):
    # [Existing 3D visualization function remains the same]
    # ... [Copy the entire existing function from the previous code]
    pass

def create_2d_map(intersecting_gdf, zonal_results):
    """
    Create 2D map with multiple layers and toggle capability
    """
    # Compute the center of the vector layer
    center = intersecting_gdf.dissolve().centroid.iloc[0]
    
    # Create base map
    m = folium.Map(location=[center.y, center.x], zoom_start=10)
    
    # Add base map layer toggle
    folium.TileLayer('openstreetmap', name='OpenStreetMap').add_to(m)
    folium.TileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', 
                     attr='Esri World Imagery', 
                     name='Esri Satellite').add_to(m)
    
    # Add vector layer
    vector_layer = folium.GeoJson(
        intersecting_gdf, 
        name='Ward Boundaries',
        style_function=lambda x: {
            'fillColor': 'transparent',
            'color': 'red',
            'weight': 2
        }
    )
    vector_layer.add_to(m)
    
    # Add min and max elevation points
    for _, row in zonal_results.iterrows():
        # Minimum elevation point
        folium.CircleMarker(
            location=[row['min_lat'], row['min_lon']],
            radius=5,
            popup=f"Minimum Elevation: {row['min_elevation']:.2f}m\nWard: {row.get('NEW_WARD_N', 'N/A')}",
            color='blue',
            fill=True,
            fillColor='blue'
        ).add_to(m)
        
        # Maximum elevation point
        folium.CircleMarker(
            location=[row['max_lat'], row['max_lon']],
            radius=5,
            popup=f"Maximum Elevation: {row['max_elevation']:.2f}m\nWard: {row.get('NEW_WARD_N', 'N/A')}",
            color='red',
            fill=True,
            fillColor='red'
        ).add_to(m)
    
    # Add layer control
    folium.LayerControl().add_to(m)
    
    return m

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
                
                # Create visualization tabs
                tab1, tab2 = st.tabs(["3D Visualization", "2D Map"])
                
                with tab1:
                    # Create 3D visualization
                    fig = create_3d_visualization(dem_file, intersecting_gdf, zonal_results)
                    st.plotly_chart(fig, use_container_width=True)
                
                with tab2:
                    # Create 2D map
                    m = create_2d_map(intersecting_gdf, zonal_results)
                    folium_static(m)

if __name__ == "__main__":
    main()
