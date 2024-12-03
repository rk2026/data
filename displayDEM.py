import streamlit as st
import geopandas as gpd
import rasterio
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from rasterio.mask import mask
from rasterio.features import geometry_mask
import requests
from io import BytesIO

def fetch_github_file(url):
    """
    Fetch file from GitHub URL
    
    Args:
        url (str): Direct download URL for the file
    
    Returns:
        bytes: File content
    """
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
    
    Args:
        dem_path (str/BytesIO): Path or file-like object of DEM raster
        vector_path (str/BytesIO): Path or file-like object of vector layer
    
    Returns:
        GeoDataFrame with zonal statistics
    """
    # Read the vector layer
    gdf = gpd.read_file(vector_path)
    
    # Read the DEM raster
    with rasterio.open(dem_path) as dem:
        # Reproject vector to match raster CRS if needed
        if gdf.crs != dem.crs:
            gdf = gdf.to_crs(dem.crs)
        
        # Initialize results list
        zonal_results = []
        
        # Process each polygon
        for idx, row in gdf.iterrows():
            try:
                # Clip raster to polygon
                out_image, out_transform = mask(dem, [row.geometry], crop=True)
                
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
                        
                        zonal_results.append({
                            'ward_id': row.get('ward_id', idx),
                            'min_elevation': min_val,
                            'max_elevation': max_val,
                            'min_lon': min_lon,
                            'min_lat': min_lat,
                            'max_lon': max_lon,
                            'max_lat': max_lat,
                            'geometry': row.geometry
                        })
            
            except Exception as e:
                st.warning(f"Could not process polygon {idx}: {e}")
        
        # Create GeoDataFrame from results
        results_gdf = gpd.GeoDataFrame(zonal_results, crs=gdf.crs)
        return results_gdf

def create_interactive_map(dem_path, vector_path, zonal_results):
    """
    Create interactive map visualization
    
    Args:
        dem_path (str/BytesIO): Path to DEM raster
        vector_path (str/BytesIO): Path to vector layer
        zonal_results (GeoDataFrame): Processed zonal statistics
    
    Returns:
        Plotly Figure
    """
    # Read vector layer for boundaries
    vector_gdf = gpd.read_file(vector_path)
    
    # Create base map
    fig = px.choropleth_mapbox(
        vector_gdf, 
        geojson=vector_gdf.geometry, 
        locations=vector_gdf.index, 
        color_continuous_scale="Viridis",
        mapbox_style="open-street-map",
        center={"lat": vector_gdf.geometry.centroid.y.mean(), 
                "lon": vector_gdf.geometry.centroid.x.mean()},
        zoom=8,
        opacity=0.5
    )
    
    # Add min and max elevation points
    for _, row in zonal_results.iterrows():
        # Min elevation point
        fig.add_trace(go.Scattermapbox(
            mode="markers",
            lon=[row['min_lon']],
            lat=[row['min_lat']],
            marker={"size": 10, "color": "blue"},
            text=f"Min Elevation: {row['min_elevation']:.2f}m",
            name="Minimum Elevation"
        ))
        
        # Max elevation point
        fig.add_trace(go.Scattermapbox(
            mode="markers",
            lon=[row['max_lon']],
            lat=[row['max_lat']],
            marker={"size": 10, "color": "red"},
            text=f"Max Elevation: {row['max_elevation']:.2f}m",
            name="Maximum Elevation"
        ))
    
    fig.update_layout(
        title="DEM Zonal Statistics Visualization",
        mapbox_style="open-street-map",
        showlegend=True,
        margin={"r":0,"t":50,"l":0,"b":0}
    )
    
    return fig

def main():
    st.title("DEM Analysis and Visualization")
    
    # GitHub URLs (replace with actual URLs)
    dem_url = st.text_input("Enter DEM File GitHub URL", 
                            "https://github.com/your_username/repo/raw/main/DHADING_Thakre.tif")
    vector_url = st.text_input("Enter Vector Layer GitHub URL", 
                               "https://github.com/your_username/repo/raw/main/Bagmati_ward.gpkg")
    
    if st.button("Process Data"):
        with st.spinner("Processing DEM and Vector Data..."):
            # Fetch files
            dem_file = fetch_github_file(dem_url)
            vector_file = fetch_github_file(vector_url)
            
            if dem_file and vector_file:
                # Process zonal statistics
                zonal_results = process_dem_zonal_stats(dem_file, vector_file)
                
                # Display results table
                st.subheader("Zonal Statistics Results")
                st.dataframe(zonal_results.drop(columns=['geometry']))
                
                # Create interactive map
                fig = create_interactive_map(dem_file, vector_file, zonal_results)
                st.plotly_chart(fig, use_container_width=True)
                
                # Optional: 3D Terrain Visualization
                st.subheader("3D Terrain Visualization")
                with rasterio.open(dem_file) as dem:
                    dem_array = dem.read(1)
                    
                    # Create 3D surface plot
                    fig_3d = go.Figure(data=[go.Surface(z=dem_array)])
                    fig_3d.update_layout(
                        title='3D Terrain Visualization',
                        scene = dict(
                            xaxis_title='X',
                            yaxis_title='Y',
                            zaxis_title='Elevation'
                        )
                    )
                    st.plotly_chart(fig_3d, use_container_width=True)

if __name__ == "__main__":
    main()
