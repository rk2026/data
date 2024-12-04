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
    """
    Create enhanced 3D visualization combining DEM and vector layers
    
    Returns:
        Plotly Figure
    """
    # Read the DEM raster
    with rasterio.open(dem_path) as dem:
        # Read raster data
        dem_array = dem.read(1)
        
        # Get raster bounds and transform
        bounds = dem.bounds
        
        # Create x and y coordinates
        x = np.linspace(bounds.left, bounds.right, dem_array.shape[1])
        y = np.linspace(bounds.bottom, bounds.top, dem_array.shape[0])
        
        # Calculate z-scale factor (typically < 1 to avoid vertical exaggeration)
        elevation_range = dem_array.max() - dem_array.min()
        lat_range = bounds.top - bounds.bottom
        lon_range = bounds.right - bounds.left
        z_scale = min(lat_range, lon_range) / elevation_range * 0.1
        
        # Create 3D surface plot of DEM with adjusted z-scale
        surface_trace = go.Surface(
            z=dem_array * z_scale, 
            x=x, 
            y=y, 
            colorscale='Viridis', 
            showscale=False,  # Remove color scale legend
            name='Terrain Elevation',
            opacity=0.7
        )
        
        # Prepare data for traces
        traces = [surface_trace]
        
        # Add vector layer boundaries
        for _, polygon in intersecting_gdf.iterrows():
            # Extract exterior coordinates of the polygon
            if polygon.geometry.type == 'Polygon':
                coords = list(polygon.geometry.exterior.coords)
            elif polygon.geometry.type == 'MultiPolygon':
                # For multipolygon, use the first polygon's exterior
                coords = list(list(polygon.geometry.geoms)[0].exterior.coords)
            
            # Get z values for the polygon boundary
            boundary_z = np.interp(
                [coord[2] if len(coord) > 2 else dem_array.mean() for coord in coords], 
                [dem_array.min(), dem_array.max()], 
                [dem_array.min() * z_scale, dem_array.max() * z_scale]
            )
            
            # Create polygon boundary trace
            boundary_trace = go.Scatter3d(
                x=[coord[0] for coord in coords],
                y=[coord[1] for coord in coords],
                z=boundary_z,
                mode='lines',
                line=dict(color='red', width=2),
                showlegend=False
            )
            traces.append(boundary_trace)
        
        # Add min and max elevation points with draped positioning
        for _, row in zonal_results.iterrows():
            # Find the closest grid point for min elevation
            min_x_idx = np.argmin(np.abs(x - row['min_lon']))
            min_y_idx = np.argmin(np.abs(y - row['min_lat']))
            min_z = dem_array[min_y_idx, min_x_idx] * z_scale
            
            # Find the closest grid point for max elevation
            max_x_idx = np.argmin(np.abs(x - row['max_lon']))
            max_y_idx = np.argmin(np.abs(y - row['max_lat']))
            max_z = dem_array[max_y_idx, max_x_idx] * z_scale
            
            # Minimum elevation point (draped)
            min_trace = go.Scatter3d(
                x=[row['min_lon']], 
                y=[row['min_lat']], 
                z=[min_z],
                mode='markers',
                marker=dict(
                    size=10,
                    color='blue',
                    opacity=0.8
                ),
                showlegend=False,
                text=f"Min Elev: {row.get('min_elevation', 'N/A')}m<br>"
                     f"Ward: {row.get('NEW_WARD_N', 'N/A')}"
            )
            traces.append(min_trace)
            
            # Maximum elevation point (draped)
            max_trace = go.Scatter3d(
                x=[row['max_lon']], 
                y=[row['max_lat']], 
                z=[max_z],
                mode='markers',
                marker=dict(
                    size=10,
                    color='red',
                    opacity=0.8
                ),
                showlegend=False,
                text=f"Max Elev: {row.get('max_elevation', 'N/A')}m<br>"
                     f"Ward: {row.get('NEW_WARD_N', 'N/A')}"
            )
            traces.append(max_trace)
        
        # Create 3D figure
        fig = go.Figure(data=traces)
        
        # Update layout for full-screen and better view
        fig.update_layout(
            title='3D Terrain Visualization',
            scene=dict(
                xaxis_title='Longitude',
                yaxis_title='Latitude',
                zaxis_title='Elevation (scaled)',
                aspectmode='manual',
                aspectratio=dict(x=1, y=1, z=0.3),  # Flatten z-axis
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=1)
                )
            ),
            height=900,  # Increased height
            width=1600,  # Significantly increased width
            margin=dict(l=0, r=0, t=30, b=0),  # Reduced margins
            showlegend=False  # Remove legend completely
        )
        
        return fig

def main():
    st.title("DEM Analysis and 3D Visualization")
    
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

                
if __name__ == "__main__":
    main()



